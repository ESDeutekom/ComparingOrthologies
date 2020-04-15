#! /usr/bin/env python3
'''
    Author: John van Dam
    Created: September 4th, 2017
    Version: 0.1.0

    Purpose of script: The purpose of this script is to provide standardized
    regular expressions for use with the eukarya set.
'''

import re
import logging
import unittest

logger = logging.getLogger(__name__)  # Set up the logger

def protein_id():
    """Regular expression pre-compiled for the eukarya protein_id."""
    return re.compile(r"(?<!\w)[A-Z]{4}\d{6}(?!\w)")

def orthofinder_id():
    """Regular expression pre-compiled for the Orthofinder Orthologous group id."""
    return re.compile(r"(?<!\w)OG\d{7}(?!\w)")

class TestRegex(unittest.TestCase):
    ''' Unit test for the regexes. '''

    def test_protein_id(self):
        self.assertTrue(protein_id().search('HSAP000001'))
        self.assertTrue(protein_id().search(' HSAP000001 '))
        self.assertFalse(protein_id().search('HSAP00001'))
        self.assertFalse(protein_id().search('HAP000001'))
        self.assertFalse(protein_id().search('HS0P000001'))
        self.assertFalse(protein_id().search('ENSP00000418196'))

    def test_orthofinder_id(self):
        self.assertTrue(orthofinder_id().search('OG0000001'))
        self.assertTrue(orthofinder_id().search(' OG0000001 '))
        self.assertFalse(orthofinder_id().search('OG000001'))
        self.assertFalse(orthofinder_id().search('OG00000001'))
        self.assertFalse(orthofinder_id().search('O00000001'))
        self.assertFalse(orthofinder_id().search('_OG0000001'))



if __name__ == '__main__':
    unittest.main()
