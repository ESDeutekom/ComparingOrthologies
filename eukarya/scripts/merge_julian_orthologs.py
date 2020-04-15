#!/usr/bin/env python3

'''This script specifically tries to condens Julian's orthologs by merging
overlapping groups of sequences. Julian's orthologs are defined per pfam family.
Thus in case of multi-domain proteins, this means that proteins can occur in
mutliple orthologous groups. This script attempt to "fix" the easily mergable
orthologous groups to make a more coherent Orthologous group definition.'''

__author__ = "John van Dam"
__copyright__ = "Copyright 2019, Snel Lab"
#__credits__ = ["John van Dam"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "John van Dam"
__email__ = "tjp.vandam@uu.nl"
__status__ = "Prototype"


# Import modules
import sys
import os
import argparse
import logging
import eukarya
from eukarya.database import engine, Session, Julian, Genes


# Setting up the logger
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)
#logger.setLevel("DEBUG")

overlap_pcnt_cutoff = 50.0
logger.info("OG overlap merge cutoff is >= {}.".format(overlap_pcnt_cutoff))

logger.info("Setting up database connection")
session = Session()

# Retrieve julians data
logger.info("Querying database.")
julians_og_table = dict()
results = session.query(Julian.og_id,Genes.protein_id).join(Genes).all()
for row in results:
    julians_og_table[(row.og_id)] = julians_og_table.get((row.og_id),set()).union({row.protein_id})

logger.info("Merging OGs which are subsets of each other")
# Find which OGS are a subset of a larger og
# Get an ordered list of og ids sorted from smallest to largest set
sorted_ogs = sorted(julians_og_table.keys(),key=lambda x: len(julians_og_table[x]))

# Go through each og id and then walk through all larger og_ids to see if it is a subset. Merge ogs if so
old_set = dict()
new_set = julians_og_table.copy()
ogs_merged = 0
itterations = 0

while old_set != new_set:  # When both dicts are identical comp() returns 0, which is quite convenient in this case
    itterations += 1
    logger.info("Started itteration {}.".format(itterations))
    old_set = new_set.copy()
    sorted_old_ogs = sorted(old_set.keys(),key=lambda x: len(old_set[x]))
    try:
        for index,og_id in enumerate(sorted_old_ogs):
            for larger_og_id in sorted_old_ogs[index:]:  # Not doing index+1 to avoid self matching because at the end this would throw an exception
                if larger_og_id == og_id:
                    continue
                else:
                    if old_set[og_id].issubset(old_set[larger_og_id]):
                        # og is subset of another, so merge in the new_set
                        logger.info("{}({}) is subset of {}({})".format(og_id,len(old_set[og_id]),larger_og_id,len(old_set[larger_og_id])))
                        del new_set[og_id]  # First attempt to delete the old one, so that on KeyError the key does not already get updated.
                        new_set[(larger_og_id,og_id)] = new_set.pop(larger_og_id)  # Updating the larger set's key to include the redundant OG id.
                        ogs_merged += 1
                    elif old_set[og_id] & old_set[larger_og_id]:
                        size_x_set = len(old_set[og_id])
                        size_y_set = len(old_set[larger_og_id])
                        size_overlap = len(old_set[og_id].intersection(old_set[larger_og_id]))
                        pcnt_overlap = size_overlap/size_x_set*100
                        if pcnt_overlap >= overlap_pcnt_cutoff:
                            logger.info("{}({}) shares overlap with {}({}) of size {}({}%)".format(og_id,size_x_set,larger_og_id,size_y_set,size_overlap,pcnt_overlap))
                            # The situation now is different in that I try to add new_set[og_id] to the larger og, so here we will throw the key error before assigning anyway
                            # And we can only delete the old og after we've added it, so here the del statement comes second.
                            # NOTE: Problem is that we already 'popped' the larger_og_id before the KeyError was thrown, which now is missing... So we need to raise the exception before we do anything.
                            if og_id not in new_set.keys():
                                raise KeyError
                            new_set[(larger_og_id,og_id)] = new_set.pop(larger_og_id).union(new_set[og_id])  # Updating the larger set's key to include the redundant OG id.
                            del new_set[og_id]
                            ogs_merged += 1
    except KeyError:
        logger.info("Could not find og_id in updated orthology dictionary, likely because the key was already updated before and removed. Restarting procedure with updated orthology dictionary.")

# Once old_Set == new_Set we know we are done, so just print it!
for ogid in sorted(new_set.keys(),key=lambda x: len(new_set[x])):
    print("{}: {}".format(ogid," ".join([str(x) for x in list(new_set[ogid])])))

logger.info("Done! Merged {} OGs.".format(ogs_merged))
#
