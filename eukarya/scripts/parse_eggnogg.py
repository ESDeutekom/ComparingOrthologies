#!/usr/bin/env python3
# coding: utf-8

# Author: John van Dam
# Email: t.j.p.vandam@uu.nl
# Creation: 11-2018
# Purpose of this script: To parse the eggnogg annotation file to properly define OGs

import sys
import os
import argparse
import logging
import re
import pandas as pd

from sqlalchemy import Column, ForeignKey, Integer, String, Boolean, create_engine, MetaData
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

import eukarya
# Get the SQL alchemy objects
from eukarya.database import Proteins, Genes, engine, Session
# Get info on the files etc.
from eukarya.files import eukarya_dir, eukarya_file, annotations_file

# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
#logger.setLevel("DEBUG")

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument("-a","--annotation_file", type=argparse.FileType('r'),
                    help="the eggnogg annotaton file", default=sys.stdin)
parser.add_argument("-o","--outputfile", type=argparse.FileType('w'), help="The reformatted eggnogg annotation file. Will print to STDOUT by default.", default=sys.stdout)
args = parser.parse_args()

og_dict = dict()  # Dictionary to store stuff.

# Open the annotation file and start processing
with(args.annotation_file):
    # The columns are:
    # 1: protein_id
    # 2: best match seed protein_id
    # 3: best protein match (e-value)
    # 4: best protein match (bit-score)
    # 5: proposed gene name based on og and best hit
    # 6: GO annotations, seperated by a comma
    # 7: KEGG KO identifiers
    # 8: BiGG metabolic reactions identifiers
    # 9: eggnoggs automatic taxonomic level determination
    # 10: list of ogs at various taxonomic levels
    # 11: Best matching OG model in HMM mode. 'NA|NA|NA' when run with diamond.
    # 12: COG functional category inferred from best matching OG
    # 13: Functional description based on OG
    # In order to get the ogs on leca level we need to extract the @euNOG
    # from column 10. Note that euNOG does not mean an og is LECA, but also
    # for instance metazoa specific, but since those are eukaryotes it gets
    # a euNOG. For our purpose this is actually quite handy.
    # All proteins that are not considered bacterial should have an euNOG.
    # NOTE: It is possible for a protein to belong to multiple euNOGs.
    singleton_counter = 1
#    print("\t".join(["protein_id","eunog_id"]),file=args.outputfile)
    for line in args.annotation_file:
        line.rstrip()  # Remove newline and trailing spaces
        if line[0] == "#":
            continue
        elements = line.rsplit("\t")
        protein_id = elements[0]
        seed_protein_id = elements[1]
        e_value = elements[2]
        bit_score = elements[3]
        proposed_gene_name = elements[4]
        go_terms = elements[5].rsplit(",")
        kegg_terms = elements[6].rsplit(",")
        bigg_terms = elements[7].rsplit(",")
        taxonomic_level = elements[8]
        eggnogg_og_list = elements[9].rsplit(",")
        eunog_ids = [nog for nog in eggnogg_og_list if "@euNOG" in nog]  # Multiple euNOGs are possible per protein
        best_og_model = elements[10]
        cog_functional_category = elements[11]
        proposed_functional_description = elements[12]

        if len(eunog_ids) == 0:  # If no euNOG is determined, create a singleton identifier
            # print("\t".join([protein_id,"S{:09}@euNOG".format(singleton_counter)]),file=args.outputfile)
            og_dict["S{:09}@euNOG".format(singleton_counter)] = [protein_id]
            singleton_counter += 1
        else:
            for euNOG in eunog_ids:
            #     print("\t".join([protein_id,euNOG]),file=args.outputfile)
                og_dict[euNOG] = og_dict.get(euNOG,[]) + [protein_id]

for eunog in sorted(og_dict.keys()):
    print("{}: {}".format(eunog," ".join(og_dict[eunog])), file=args.outputfile)
#        print("\t".join(([protein_id] + eunog_id)),file=args.outputfile)
