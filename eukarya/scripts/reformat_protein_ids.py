#! /usr/bin/env python3
'''
    Author: John van Dam
    Created: September 1th, 2017
    Version: 0.1.0

    Purpose of script: The purpose of this script is to search and replace
    the uniform eukarya IDs with more metadata like gene symbol, species, etc.
'''
import sys
import os
import argparse
import logging
import re
import eukarya
import eukarya.regex
import sqlalchemy
import numpy
# Get the SQL alchemy objects
from eukarya.database import Proteins, Genes, Species, OrthofinderDiamond, Session, engine

# TODO: Change the "More results found than expected" warning into a critical
#       error, once the orthofinder table is corrected.
# TODO: I should maybe rethink this for large files, though processing
#       eukarya_proteomes_lt.fa is possible, so no requirement.


# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
logger.setLevel("DEBUG")

# Define the queryable columns
queryable_columns = [Proteins.protein_id,
                     Proteins.original_protein_id,
                     Proteins.original_transcript_id,
                     Genes.gene_symbol,
                     Genes.original_gene_id,
                     Species.taxonomy_id,
                     Species.abbreviation,
                     Species.scientific_name,
                     Species.common_name,
                     Species.supergroup,
                     Species.relevant_taxonomy,
                     OrthofinderDiamond.og_id]

# Deriving actual column names from SQLalchemy columns defined in queryable_columns
allowed_tags = [str(column.name) for column in queryable_columns]

parser = argparse.ArgumentParser()
parser.add_argument("file", nargs='?', metavar="FILE",
                    type=argparse.FileType('r'),
                    help="The file containing Eukarya IDs, or STDIN if ommitted.",
                    default=sys.stdin)
parser.add_argument("-p","--pattern", type=str,
                    help="The specific pattern to replace things with. Tags\
                    that can be used are:\n{}.".format(", ".join(["{%s}" % x for x in allowed_tags])),
                    default='{protein_id} {gene_symbol} {original_protein_id} {scientific_name} {og_id}')
parser.add_argument("-s","--sanitize", action='store_true', help='When set, replacements will be sanitized for troublesome characters [()[]\{\};:''"",.]. ')
args = parser.parse_args()


# Check if we are using stdin or a file to read data from
if args.file.name == '<stdin>':
    # Now, this actually only makes sense when we pipe stuff, we don;t want to
    # type things ourselves. So check if stdin is connected to tty or not
    if sys.stdin.isatty():
        logger.info("<stdin> is connected to a terminal, not file or pipe. Exiting.")
        parser.print_help()
        exit(1)
    else:
        logger.info("Reading data from STDIN.")
else:
    logger.info("Reading data from file '{}'".format(args.file.name))

# Test if there are any issues with the user provided pattern
# I can do this simply by challenging the pattern and catch any exceptions
test_tags = dict(zip(allowed_tags, ["test"]*len(allowed_tags)))
try:
    args.pattern.format(**test_tags)
except Exception as e:
    # Something else is the matter with this format string
    logger.error("Error in pattern: {}.".format(str(e)))
    exit(1)
logger.info("Using formatting pattern:")
logger.info(args.pattern)


# Read in the file completely, so we can prefetch the IDs we want to look up.
# Argparse already takes care of creating a file handle for us.
to_be_processed_file = args.file.read()
original_data = to_be_processed_file
# Get all protein_ids from the to_be_processed_file
all_protein_ids = set(eukarya.regex.protein_id().findall(to_be_processed_file))
logger.info("Found {} protein_ids to replace".format(len(all_protein_ids)))
logger.debug(" ".join(all_protein_ids))
if 0 == len(all_protein_ids):
    logger.info("Nothing to do. Outputting data without modifications.")
    sys.stdout.write(to_be_processed_file)  # Output the unmodified data
    exit(0)

# Do the query
# TODO: I should rethink this for large files.
# Note on SQLite:
# Apparently some env variable SQLITE_MAX_VARIABLE_NUMBER is set to 999 as default.
# The solution is to use session.bind.dialect.name == "sqlite" to determine
# if we are using sqlite and if so do multiple queries of 999 Something
# variables and concatenate the results.
# Solved: For some reason the chunk list comprehension is slow for large number of ids.
# query_chunks = [list(all_protein_ids)[x:x+999] for x in range(0, len(all_protein_ids), 999)]

logger.info("Querying database.")
try:
    session = Session()  # Start session
    # The below query retrieves all if not modified by a filter() method.
    query = session.query(*queryable_columns).\
                join(Genes, Genes.gene_id==Proteins.gene_id).join(Species).outerjoin(OrthofinderDiamond).distinct()  # The outerjoin is necessary because Orthofinder was only performed on the core set
    query = query.filter(Proteins.protein_id.in_(all_protein_ids))
    logger.debug(query.statement)
    results = query.all()
    session.close()  # Close the DB session
except sqlalchemy.orm.exc.NoResultFound:
    # Check if we even have results, if not then return the original data
    logger.warning("No results found matching your protein_ids. Returning the data unmodified.")
    sys.stdout.write(to_be_processed_file)
    exit(0)
except Exception as e:
    # General database exceptions handled here
    logger.error("Could not perform a database query. Exiting.")
    logger.debug(e)
    exit(1)
logger.debug("Query finished with {} results.".format(len(results)))

# Warn if some identifiers were not found in the database, this suggests that
# we might be finding matches with non-protein ids.
if len(all_protein_ids) != len(results):
    logger.warning("Number of protein_ids ({0}) does not match results found in the database ({1})!".format(len(all_protein_ids),len(results)))
    if len(all_protein_ids) > len(results):
        logger.warning("Missing identifiers in database: '{}'".format("', ''".join(sorted(all_protein_ids-set([x.protein_id for x in results])))))
    else:
        logger.warning("More results found than expected. (Currently a verified cause is an incorrect mapping of genes to multiple ogs)")
        # The above needs to become an logger.error() in future when OG issue is solved.

# Search and replace
logger.info("Search and replacing {} identifiers using the supplied pattern.".format(len(results)))
# First build a dictionary of replacement strings to speed things up
replacements = dict()
for row in results:
    # For every entry search and replace using the specified format.
    row_dict = row._asdict()  # Convert the result object into a dictionary
    # If sanitize flag is set convert spaces in each dictionary value into underscores
    if args.sanitize:
        for key, value in row_dict.items():
            row_dict[key] = replacement_string = re.sub(r'[\(\)\[\]\{\}\;\:\'\"\,\.\*\\\/\s]',"_",str(value))

    # Preparing replacement string
    replacement_string = args.pattern.format(**row_dict)

    replacements[row_dict["protein_id"]] = replacement_string

# In re.sub I use lambda function so I have to do the search/replace only once
to_be_processed_file = re.sub(eukarya.regex.protein_id(),
                              (lambda m: replacements.get(m.group(0),m.group(0))),
                              to_be_processed_file)

sys.stdout.write(to_be_processed_file)  # Output the modified data
