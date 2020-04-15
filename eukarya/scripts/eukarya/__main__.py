#! /usr/bin/env python3
'''
    Author: John van Dam
    Created: April 6th, 2017
    Version: 0.1.0

    Purpose of module: collection of classes and functions for use in eukarya related projects.

'''

import logging
import eukarya

if __name__ == '__main__':
    logger = logging.getLogger(__name__)  # Set up the logger
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    # Prints Git hash for current versions of data and code
    eukarya.version.log_repo_versions()
