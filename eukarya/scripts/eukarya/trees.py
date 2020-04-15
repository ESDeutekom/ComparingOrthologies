#! /usr/bin/env python3
'''
    Author: John van Dam
    Created: November 23, 2017
    Version: 0.1.0

    Purpose of module: collection of classes and functions for use in eukarya
    related projects.

'''

import logging
from eukarya.files import get_eukarya_file_path
from ete3 import Tree
import unittest

logger = logging.getLogger(__name__)  # Set up the logger
# logger.setLevel("DEBUG")

def get_ete_Tree(tree_label):
    '''
    Returns a ETE3 Tree object for the tree file specified by the file key
    as found in the file module.
    '''
    tree_file = get_eukarya_file_path(tree_label)
    logger.debug("Attempting to load an ETE Tree from {}.".format(tree_file))
    try:
        tree = Tree(tree_file,format=1)
    except:
        raise
    return tree

def tree_abbreviated():
    '''
    Returns the ETE3 object for the full eukarya tree with abbreviated species
    names. '''
    return get_ete_Tree('tree_abbrv')

def tree_taxid():
    '''
    Returns the ETE3 object for the full eukarya tree with abbreviated species
    names. '''
    return get_ete_Tree('tree_taxid')

def tree_abbreviated_annotated():
    '''
    Returns the ETE3 object for the full bifurcating eukarya tree with
    abbreviated species names and labeled internal nodes.
    '''
    return get_ete_Tree('tree_abbrv_annotated')

tree_abbreviated_bifurcating = tree_abbreviated_annotated  # function alias

def tree_abbreviated_multifurcating():
    '''
    Returns the ETE3 object for the full multifurcating eukarya tree with abbreviated species
    names and labeled internal nodes.
    '''
    return get_ete_Tree('tree_abbrv_annotated')

def tree_abbreviated_core_annotated():
    '''
    Returns the ETE3 object for the core bifurcating eukarya tree with
    abbreviated species names and labeled internal nodes.
    '''
    return get_ete_Tree('tree_abbrv_core_annotated')

def tree_fullnames():
    '''
    Returns the ETE3 object for the full bifurcating eukarya tree with
    full scientific species names.
    '''
    return get_ete_Tree('tree_fullnames')


class TestTrees(unittest.TestCase):
    ''' Unit test for eukarya trees. '''

    def test_trees(self):
        self.assertIsInstance(get_ete_Tree('tree_abbrv'),Tree)
        self.assertIsInstance(tree_taxid(),Tree)
        self.assertIsInstance(tree_abbreviated_annotated(),Tree)
        self.assertIsInstance(tree_abbreviated_bifurcating(),Tree)
        self.assertIsInstance(tree_abbreviated_multifurcating(),Tree)
        self.assertIsInstance(tree_abbreviated_core_annotated(),Tree)
        self.assertIsInstance(tree_abbreviated_core_annotated(),Tree)
        self.assertIsInstance(tree_fullnames(),Tree)

def main():
    unittest.main()


    # Some general stuff
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    try:
        logger.info(get_ete_Tree('tree_abbrv').write())
    except Exception as e:
        logger.info("Something did not work using key 'tree_abbrv'")
        logger.info(e)
    try:
        logger.info(get_ete_Tree('orthofinder').write())
    except Exception as e:
        logger.info("Something did not work using key 'orthofinder'")
        logger.info(e)


if __name__ == '__main__':
    main()
