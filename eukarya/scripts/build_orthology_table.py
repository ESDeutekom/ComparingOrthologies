#!/usr/bin/env python3
# coding: utf-8

# Author: John van Dam
# Email: t.j.p.vandam@uu.nl
# Creation: 10-2017
# Purpose of this script: To build a SQLite3 database of Eukarya dataset related
#                         orthology datasets

import sys
import os
import argparse
import logging
import re
import pandas as pd
import argparse

from sqlalchemy import Column, ForeignKey, Integer, String, Boolean, create_engine, MetaData
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

import eukarya
# Get the SQL alchemy objects
from eukarya.database import Proteins, Genes, get_orthology_leca_tables, engine, Session
# Get info on the files etc.
from eukarya.files import eukarya_dir, eukarya_file, annotations_file

# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
#logger.setLevel("DEBUG")

# Dictionary holding all method key and sqlalchemy classes
# methods = {
#     'orthofinder_diamond': OrthofinderDiamond,
#     'orthofinder_blast_e+1': OrthofinderBlast,
#     'orthofinder_blast_e-1': OrthofinderBlast_1,
#     'orthofinder_blast_e-3': OrthofinderBlast_3,
#     'eggnog': Eggnog,
#     'julian_all': Julian,
#     'julian_all_90': Julian90,
#     'julian_all_50': Julian50,
#     'panther': Panther,
#     'manual': Manual
# }
methods = get_orthology_leca_tables()  # This is a dictionary with for each orthology both table objects

# Get arguments to determine which og we need to make a DB for
parser = argparse.ArgumentParser()
parser.add_argument("ortholog_method", choices=methods.keys(), help="Orthologous group method.")
args = parser.parse_args()

logger.info("Building a database table for {}".format(args.ortholog_method))
OGSQLClass, LECASQLClass = methods[args.ortholog_method]

logger.info("Preloading Proteins information from Eukarya database.")
session = Session()

try:
    protein2gene = dict()
    for (protein_id,gene_id) in session.query(Proteins.protein_id,Proteins.gene_id).all():
        protein2gene[protein_id] = gene_id
except:
    logger.error("Undefined error while querying database.")
    exit(1)

# Parse the orthofinder data
print(eukarya_dir+annotations_file[args.ortholog_method])
logger.info("Loading Orthology information from file...")
ortho_data = list()
try:
    ortho_fh = open(eukarya_dir+annotations_file[args.ortholog_method],"r")
except IOError:
    logger.error("Unable to read in the Orthology file.")
    exit(1)
for line in ortho_fh:
    line = line.rstrip()
    elements = line.split(": ")  # and not .split("\s") because it is not .rsplit().
    ortho_id = elements[0]  # orthofinder id is first element
    protein_ids = elements[1].split(" ")
    for protein_id in protein_ids:
        # Since we only want to match orthology info on the gene level, we look up the corresponding gene id
        try:
            gene_id = protein2gene[protein_id]
        except KeyError:
            if args.ortholog_method == "manual":
                logger.warning("Protein ID %s does not have a gene id. Skipping because manual ogs contains some non-standard protein ids!" % protein_id)
            else:
                logger.error("Protein ID %s does not have a gene id. This is a critical error!" % protein_id)
                exit(1)
        ortho_data.append((gene_id,ortho_id))
ortho_fh.close()
logger.info("Done.")

logger.info("Loading LECA ogs from file...")
leca_data = list()
try:
    leca_fh = open(eukarya_dir+annotations_file[args.ortholog_method+"_LECA"],"r")
except IOError:
    logger.error("Unable to read in the LECA ogs file.")
    exit(1)

for line in leca_fh:
    ortho_id = line.rstrip()
    leca_data.append((ortho_id))
leca_fh.close()
logger.info("Done.")

# Now make it a pandas object because it is easier to write it away that way
Ortho_df = pd.DataFrame(ortho_data, columns=('gene_id','og_id'))
LECA_df = pd.DataFrame(leca_data, columns=['og_id'])

# Since the first version of orthofinder was run on a protein set containing
# multiple proteins per gene, we need to make the entries unique
Ortho_df = Ortho_df.drop_duplicates()
LECA_df = LECA_df.drop_duplicates()

# Sort by gene _id, ascending
Ortho_df = Ortho_df.sort_values(by='gene_id')

# Sort by og_id, ascending
LECA_df = LECA_df.sort_values(by='og_id')

# Drop all tables if exists
logger.info("Dropping the relevant Ortholog tables, if exists.")
OGSQLClass.__table__.drop(checkfirst=True)
LECASQLClass.__table__.drop(checkfirst=True)
#Annotations.metadata.drop_all()

# Create all tables
logger.info("Creating Ortholog tables.")
OGSQLClass.__table__.create()
LECASQLClass.__table__.create()
#Annotations.metadata.create_all()

# Use the Pandas.to_sql() function to populate the OG table
logger.info("Populating Ortholog table.")
Ortho_df[['gene_id','og_id']].to_sql(name=OGSQLClass.__table__.name, schema="Annotations", con=engine,
                                                        if_exists='append', index=False)
LECA_df[['og_id']].to_sql(name=LECASQLClass.__table__.name, schema="Annotations", con=engine,
                                                        if_exists='append', index=False)
logger.info("Stored Ortholog Table.")
