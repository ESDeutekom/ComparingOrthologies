#!/usr/bin/env python3
'''
    Author: John van Dam
    Created: August 28th, 2017
    Version: 0.1.0

    Purpose of script: The purpose of this script is to parse the fasta file
    and the hmmsearch output file and add a * in front of every hit in the hmm
    which is already in the fasta file

'''
import sys
import os
import argparse
import logging
from eukarya.regex import protein_id as regex_protein_id

# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
# logger.setLevel("DEBUG")

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", type=argparse.FileType('r'),
                    help="Fasta file with previously found sequences.")
args = parser.parse_args()

previously_found = set()
with open(args.fasta_file) as fh:
    for record in SeqIO.parse(fh, "fasta"):
        my_id = regex_protein_id.search(record.id)
        if my_id:
            previously_found.add(my_id.group(0))

print(previously_found)

for line in sys.stdin:
    line = line.rstrip()
    match = regex_protein_id.search(line)
    if match:
        logger.debug("Matched " + match.group(0))
        if match.group(0) in previously_found:
            print("*"+line)
        else:
            print(" "+line)
    else:
        print(" "+line)
