#!/usr/bin/env python3

import sys
import os
import argparse
import logging
import zipfile
import re
import pandas as pd

# Load the Eukarya species data using the eukarya database
import eukarya
from eukarya.database import Base, Proteins, Genes, Species, engine, Session
from eukarya.database import OrthofinderDiamond, OrthofinderBlast, OrthofinderBlast_1, OrthofinderBlast_3, Eggnog, Julian, Panther, Manual

# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger

methods = {
    'orthofinder_diamond': OrthofinderDiamond,
    'orthofinder_blast_e+1': OrthofinderBlast,
    'orthofinder_blast_e-1': OrthofinderBlast_1,
    'orthofinder_blast_e-3': OrthofinderBlast_3,
    'eggnog': Eggnog,
    'julian': Julian,
    'panther': Panther,
    'manual': Manual,
}

# Get arguments to determine which og we need to make a DB for
parser = argparse.ArgumentParser()
parser.add_argument("ortholog_method", choices=methods.keys(), help="Orthologous group method.")
args = parser.parse_args()

# Setting the orthology definition to use from the arguments
logger.info("Using the {} orthology.".format(args.ortholog_method))
OrthoDef = methods[args.ortholog_method]

# Start database session
session = Session()

# Let's go straight for the profiles using sqlquery
logger.info("Setting up the sql query.")
query = session.query(OrthoDef.og_id, Species.abbreviation).\
            join(Proteins, OrthoDef.gene_id==Proteins.gene_id).\
            join(Species, Species.taxonomy_id==Proteins.taxonomy_id).order_by(OrthoDef.og_id).distinct()
logger.info(query.statement)

# Let's not use pandas right now as it is not the format that will be most convienient. Lets load as list of tubles.
logger.info("Querying the database")
results = query.all()
profiles = {}  # empty dict
for (ogid,abbreviation) in results:
    if ogid not in profiles:
        profiles[ogid] = set()
    profiles[ogid].add(abbreviation)

# Now get the correct species order:
# First get the tree
logger.info("Loading the Eukarya tree (core).")
tree = eukarya.trees.get_ete_Tree('tree_abbrv_core_annotated')
# And the leaf order:
species_order = tree.get_leaf_names()

# Now build the profiles
logger.info("Building profiles")
profile_vectors = {}
for ogid in profiles:
    species_in_profile = profiles[ogid]
    if "HSAP" not in species_in_profile:  # We only need profiles with human in them anyway
        next
    profile_vectors[ogid] = [0]*len(species_order)
    for i,species in enumerate(species_order):
        if species in species_in_profile:
            profile_vectors[ogid][i] = 1  # +=1 for counting number of genes per species

logger.info("Printing")
# Now print the profiles
print ("Orthologous_group_ID:\t" + "\t".join(species_order))
for ogid in profile_vectors:
    vector = "\t".join([str(x) for x in profile_vectors[ogid]])
    print(ogid + ":\t" + vector)
