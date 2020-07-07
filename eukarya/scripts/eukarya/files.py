#! /usr/bin/env python3
'''
    Author: John van Dam
    Created: April 6th, 2017
    Version: 0.1.0

    Purpose of module: collection of classes and functions for use in eukarya related projects.

'''

import sys
import os
import logging
import re
from Bio import SeqIO
from tempfile import TemporaryDirectory
import unittest

logger = logging.getLogger(__name__)  # Set up the logger
# logger.setLevel("DEBUG")

# The origin path of the eukarya set.
eukarya_origin_path = "/hosts/linuxhome/siluur/john/Data/eukarya/"

# The eukarya directory
eukarya_dir = os.environ.get("EUKARYA_PATH",eukarya_origin_path)
print(eukarya_dir)
# The files
# Core eukarya files
eukarya_file = dict()
eukarya_file['fasta_lt'] = 'data_set/eukarya_proteomes_lt.fa'
eukarya_file['fasta_lt_core'] = 'data_set/eukarya_proteomes_lt_core.fa'
eukarya_file['metadata'] = 'data_set/eukarya_proteomes_metadata.tsv'
eukarya_file['species'] = 'data_set/eukarya_species.tsv'
eukarya_file['sqlite3'] = 'database/eukarya_db.sqlite3'
eukarya_file['tree_abbrv'] = 'data_set/trees/eukarya_species.newick'
eukarya_file['tree_taxid'] = 'data_set/trees/eukarya_species.taxIDs.newick'
eukarya_file['tree_abbrv_annotated'] = 'data_set/trees/eukarya_species.annotated.newick'
eukarya_file['tree_abbrv_multifurcating_annotated'] = 'data_set/trees/eukarya_species.multifurcating.annotated.newick'
eukarya_file['tree_abbrv_core_annotated'] = 'data_set/trees/eukarya_species.core.annotated.newick'
eukarya_file['tree_fullnames'] = 'data_set/trees/eukarya_species.names.newick'

# Files containing annotations
annotations_file = dict()
annotations_file['orthofinder_diamond_e-3'] = 'annotations/Orthogroups_orthofinder_diamond_e-3.txt'
annotations_file['orthofinder_diamond'] = annotations_file['orthofinder_diamond_e-3']  # Legacy
annotations_file['eggnog_diamond'] = 'annotations/Orthogroups_eggnog_diamond.txt'
annotations_file['eggnog'] = annotations_file['eggnog_diamond']  # Legacy
annotations_file['orthofinder_blast_e-3'] = 'annotations/Orthogroups_orthofinder_blast_e-3.txt'
annotations_file['panther_corrected'] = 'annotations/Orthogroups_panther_different.txt'
annotations_file['manual'] = 'annotations/manual_ogs/Orthogroups_small_scale.txt'
annotations_file['sqlite3'] = 'database/eukarya_annotations_db.sqlite3'
annotations_file['eggnog_hmmer_corrected'] = 'annotations/Orthogroups_eggnog_hmmer_corrected.txt'
annotations_file['broccoli'] = 'annotations/Orthogroups_broccoli.txt'
annotations_file['swiftortho'] = 'annotations/Orthogroups_Swiftortho_c50.txt'
annotations_file['sonicparanoid_sensitive'] = 'annotations/Orthogroups_Sonicparanoid_sensitive.txt'

# Files containing leca ogs, for the above files
annotations_file['orthofinder_diamond_e-3_LECA'] = "annotations/leca_orthologous_group_list_orthofinder_diamond_e-3"
annotations_file['orthofinder_blast_e-3_LECA'] = "annotations/leca_orthologous_group_list_orthofinder_blast_e-3"
annotations_file['eggnog_diamond_LECA'] = "annotations/leca_orthologous_group_list_eggnog_diamond"
annotations_file['eggnog_hmmer_corrected_LECA'] = "annotations/leca_orthologous_group_list_eggnog_hmmer_corrected"
annotations_file['panther_corrected_LECA'] = "annotations/leca_orthologous_group_list_panther_different"
annotations_file['manual_LECA'] = 'annotations/manual_ogs/leca_orthologous_group_list_small_scale.tsv'
annotations_file['broccoli_LECA'] = 'annotations/leca_orthologous_group_list_broccoli'
annotations_file['swiftortho_LECA'] = 'annotations/leca_orthologous_group_list_Swiftortho_c50'
annotations_file['sonicparanoid_sensitive_LECA'] = 'annotations/leca_orthologous_group_list_Sonicparanoid_sensitive'


# Setup the path to a local copy. Make a temporary one if not set.
if "EUKARYA_PATH" in os.environ:
    logger.debug("EUKARYA environment variable set, using local copy.")
    eukarya_path = os.environ['EUKARYA_PATH']
    eukarya_path_type = 'local'
else:
    logger.warning("EUKARYA environment variable is not set. Local temporary copies will be created as required.")
    eukarya_tempdir = TemporaryDirectory()
    eukarya_path = eukarya_tempdir.name
    eukarya_path_type = "temp"
    logger.info("Working directory set to %s." % eukarya_path)

eukarya_database = eukarya_dir + eukarya_file['sqlite3']
annotations_database = eukarya_dir + annotations_file['sqlite3']


def get_local_eukarya_dir():
    '''Returns the local eukarya directory, or raises error.'''
    if (eukarya_dir == eukarya_origin_path):
        logger.error("The eukarya root directory is not local!")
        raise
    else:
        return eukarya_dir


def get_file_path(file_key):
    ''' Returns the full path of a registered file based on the dictionary key.
    You can specify if you a file from eukarya or annotations by setting the
    second argument to 'eukarya' or 'annotations'.
    So far this is only required to get the correct path for the annotations
    SQLite DB. If a file key matches in both eukarya and annotations, the one in
    eukarya is returned.
    '''
    try:
        return eukarya_dir + eukarya_file[file_key]
    except KeyError:
        try:
            return eukarya_dir + annotations_file[file_key]
        except KeyError:
            logger.error("Eukarya file key '{}' is unknown!".format(file_key))
            raise Exception


def get_eukarya_file_path(file_key):
    ''' See get_file_path, but this only returns file paths for the eukarya set
    '''
    try:
        return eukarya_dir + eukarya_file[file_key]
    except KeyError:
        logger.error("Eukarya file key '{}' is unknown!".format(file_key))
        raise


def get_annotations_file_path(file_key):
    ''' See get_file_path, but this only returns file paths for annotation data
    '''
    try:
        return eukarya_dir + annotations_file[file_key]
    except KeyError:
        logger.error("Annotations file key '{}' is unknown!".format(file_key))
        raise


class TestFunctions(unittest.TestCase):
    ''' Unit test to test individual functions.'''

    def test_get_file_path(self):
        ''' Test get_file_path() function. '''
        self.assertTrue(os.path.isfile(get_file_path('species')))


class TestFiles(unittest.TestCase):
    ''' Unit test to test if all file references are correct. '''

    def test_eukarya_file_paths(self):
        ''' Test get_eukarya_file_path() function for all entries. '''
        for key in eukarya_file:
            self.assertTrue(os.path.isfile(get_eukarya_file_path(key)))

    def test_annotation_file_paths(self):
        ''' Test get_annotations_file_path() function for all entries. '''
        for key in annotations_file:
            self.assertTrue(os.path.isfile(get_annotations_file_path(key)))

if __name__ == '__main__':
    unittest.main()
