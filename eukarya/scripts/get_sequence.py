'''
    Author: John van Dam
    Created: April 6th, 2017
    Version: 0.1.0

'''
import sys
import os
import argparse
import logging
import re
from Bio import SeqIO
import pandas as pd
import eukarya
from eukarya.files import eukarya_dir, eukarya_file


# The eukarya directory
eukarya_dir = eukarya.files.eukarya_dir

# The files
eukarya_file_proteomes = 'data_set/eukarya_proteomes_lt.fa'
eukarya_file_sqlite3 = 'database/eukarya_db.sqlite3'

# Some general stuff
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
# logger.setLevel("DEBUG")

# Functions!
def checkArgFileExists(path):
    ''' For ArgumentParser, Checks if path is actually a file. '''
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError("File %s does not exist." % path)
    return path

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("--fasta_file", type=argparse.FileType('r'),
                    help="fasta file", default=open(eukarya_dir+eukarya_file['fasta_lt'],'r'))
parser.add_argument("protein_id", nargs='*', help="List of protein ids you want to retrieve.")
args = parser.parse_args()

query = set()

id_format = eukarya.regex.protein_id()
if (len(args.protein_id) == 0):
    # No ids as arguments given, so take from stdin
    logger.info("Reading from STDIN.")
    for line in sys.stdin:
        query |= set(id_format.findall(line))
else:
    for i in args.protein_id:
        if id_format.match(i):
            query.add(i)
        else:
            logger.warning("Protein id " + i + " is not in the correct format. Skipping.")

if (len(query) == 0):
    logger.warning("No correct identifiers provided as arguments or STDIN. No work to be done. Exiting.")
    exit(0)

#logger.info("Extracting {} fasta sequences from {}.".format(len(query),args.fasta_file.filename))
logger.info("Extracting {} fasta sequences.".format(len(query)))
for record in SeqIO.parse(args.fasta_file, "fasta"):
    if (record.id in query):
        SeqIO.write(record, sys.stdout, "fasta")
        query.remove(record.id)
    else:
        pass  # This is not the sequence we're looking for
    if (len(query) == 0):  # If our ids are done, there's no point in searching on.
        break

if not len(query) == 0:
    logger.warning("{} identifiers were not found: {}".format(len(query), ", ".join(query)))
    logger.warning("Did you search in the full proteome set or the _lt set?")
