'''
    Author: John van Dam
    Created: April 3rd, 2017
    Version: 0.1.0

    Purpose of script: The purpose of this script is to parse the Orthofinder
    ortholog output file Orthogroups.txt and reformat it so we can more easily
    grep it and is database reasy.
'''

import sys
import os
import argparse
import logging
import re

''' Main program '''

# Some general stuff
logging.basicConfig(level=logging.INFO,
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
parser.add_argument("file", metavar="FILE", type=checkArgFileExists,
                    help="The Orthogroups.txt file path")
try:
    args = parser.parse_args()
except:
    logger.error("Provided path is not a file, or is not writable.")
    exit(1)

with open(args.file) as fh:
    print("protein_id", "orthogroup", sep="\t")
    for line in fh:
        line = line.rstrip()  # Remove trailing white space and end-of-line
        # split line by space
        orthogroup_def = line.split(" ")
        # take first element (ortho group id), and remove the trailing colon.
        og_id = orthogroup_def.pop(0).replace(":","")
        # Now for each successive protein id, print line
        for prot_id in orthogroup_def:
            print(prot_id, og_id, sep="\t")
