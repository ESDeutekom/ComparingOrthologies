#! /usr/bin/env python3
'''
    Author: John van Dam
    Created: April 6th, 2017
    Version: 0.1.0

    Purpose of script: The purpose of this script is to parse the available
    metadata of the eukarya set and some axillary data files and build a sqlite3
    database (not including the sequences themselves, that's just to much).

    The script will load all metadata and auxillary data using pandas and
    produce tables to put into SQL.
'''
import sys
import os
import argparse
import logging
import re
from Bio import SeqIO
import sqlalchemy
from eukarya.database import Genes, OrthofinderDiamond, Session, engine
from eukarya.files import eukarya_dir, eukarya_file, get_file_path

# The files
eukarya_file_proteomes = get_file_path('fasta_lt_core')

# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
#logger.setLevel("DEBUG")

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("OGs", nargs='+', help="List of orthologous groups you want to retrieve.")
args = parser.parse_args()

logger.debug("OGs requested: {}".format(", ".join(args.OGs)))

ogs = set(args.OGs)

try:
    session = Session()
    query = session.query(Genes.protein_id).join(OrthofinderDiamond).filter(OrthofinderDiamond.og_id.in_(ogs)).distinct()
    logger.debug(query.statement)
    results = query.all()
    session.close()
except sqlalchemy.orm.exc.NoResultFound:
    logger.error("No orthologous groups matched {}".format(", ".join(ogs)))
    exit(1)
except Exception as e:
    logger.error("Error querying database.")
    logger.error(e)
    exit(1)
logger.debug("Query finished with {} results.".format(len(results)))

protein_ids = [result.protein_id for result in results]

logger.info("Extracting {} fasta sequences from {}.".format(len(protein_ids),eukarya_file_proteomes))
with open(eukarya_file_proteomes) as fh:
    for record in SeqIO.parse(fh, "fasta"):
        if (record.id in protein_ids):
            SeqIO.write(record, sys.stdout, "fasta")
            protein_ids.remove(record.id)
        else:
            pass  # This is not the sequence we're looking for
        if (len(protein_ids) == 0):
            break  # Where done so we can stop now.
