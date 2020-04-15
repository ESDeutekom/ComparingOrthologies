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
import subprocess
from git import Repo
import unittest

from eukarya.files import eukarya_dir

logger = logging.getLogger(__name__)  # Set up the logger

# Setup the repo objects
def _get_repo_versions():
    ''' Provides dictionary with git repo names and their current version. '''
    # The base repo:
    eukarya_repo = Repo(eukarya_dir)
    assert not eukarya_repo.bare

    versions = dict()
    versions['eukarya'] = eukarya_repo.git.describe("--always")
    for submod in eukarya_repo.submodules:
        versions[submod.name] = submod.module().git.describe("--always")

    return versions

def log_repo_versions():
    repo_versions = _get_repo_versions()
    for name,version in repo_versions.items():
        logger.info("Git Repo version for {}: {}".format(name,version))


class TestVersion(unittest.TestCase):
    ''' unit tests for the version module '''

    def test_get_repo_versions(self):
        self.assertIsInstance(_get_repo_versions(),dict)

    def test_log_repo_versions(self):
        log_repo_versions()
        self.assertLogs(logger=logger, level='INFO')


if __name__ == '__main__':
    unittest.main()
