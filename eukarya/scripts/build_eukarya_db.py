#!/usr/bin/env python
# coding: utf-8

# Author: John van Dam
# Email: t.j.p.vandam@uu.nl
# Creation: 10-2017
# Purpose of this script: To build a SQLite3 database of the Eukarya dataset

import sys
import os
import argparse
import logging
import re
from Bio import SeqIO
import pandas as pd
import numpy as np

from sqlalchemy import Column, ForeignKey, Integer, String, Boolean, create_engine, MetaData
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

import eukarya
import eukarya.regex
# Get the SQL alchemy objects
from eukarya.database import Base, Proteins, Genes, Species, engine
# Get info on the files etc.
from eukarya.files import eukarya_dir, eukarya_file

# Some general stuff
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)-15s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)  # Set up the logger
#logger.setLevel("DEBUG")

# Load files in pandas DataFrame
logger.info("Loading metadata and species information...")
eukarya_meta_data = pd.read_csv(eukarya_dir+eukarya_file['metadata'], sep="\t",
                                dtype={'taxonomy_id':int, 'gene_symbol':str, 'gene_location':str,
                                       'description':str})
eukarya_species = pd.read_csv(eukarya_dir+eukarya_file['species'], sep="\t")


# Prepare the primary tables
logger.info("Preparing primary tables.")
proteins = eukarya_meta_data[['protein_id',
                              'original_protein_id',
                              'original_transcript_id',
                              'original_gene_id',
                              'taxonomy_id',
                              'sequence_length',
                              'md5_sum',
                              'longest_transcript',
                              'original_header']]
# Note: 'original_gene_id' should be dropped from 'proteins' later on.

genes = eukarya_meta_data[['original_gene_id',
                           'gene_symbol',
                           'taxonomy_id',
                           'protein_id',
                           'longest_transcript']]

species = eukarya_species[['Taxonomy ID',
                           'Abbreviation',
                           'eukaryotic supergroup',
                           'relevant taxonomy',
                           'Scientific name',
                           'Common name',
                           'core set',
                           'Download location',
                           'Download date (yyyymmdd)',
                           'Genome version',
                           'Reference',
                           'Notes',
                           'fasta_header_type',
                           'proteome_fasta_file',
                           'genome_fasta_file']]

# Make SQL compatible headers.
species.columns = ['taxonomy_id',
                   'abbreviation',
                   'supergroup',
                   'relevant_taxonomy',
                   'scientific_name',
                   'common_name',
                   'core',
                   'download_location',
                   'download_date_yyyymmdd',
                   'assembly_version',
                   'reference',
                   'notes',
                   'fasta_header_type',
                   'proteome_fasta_file',
                   'genome_fasta_file']

# proteins column names are already suitable
#proteins.columns = ['protein_id', 'original_protein_id', 'original_transcript_id',
#       'taxonomy_id', 'sequence_length', 'md5_sum', 'longest_transcript',
#       'original_header']

# genes column names are already suitable
#genes.columns = ['original_gene_id', 'gene_symbol', 'taxonomy_id', 'protein_id',
#       'longest_transcript']

''' Now we need to reduce gene_table to unique entries per gene, keep a reference
to the protein id of the longest transcript and device a primary key/secondary key schema
to use to connect both tables.
'''
logger.info("Optimizing Genes table.")
# Take only the entries of the longest transcript. As they should be unique per gene.
genes = genes.loc[genes['longest_transcript'] == 1]
# Now drop longest_transcript column
del genes['longest_transcript']

# Create a standardized gene_id for use in the primary_key,secondary_key schema
# Decided on the standard primary_id solution of mysql, namely intergers in range of 1 to length of table
genes['gene_id']=range(1,len(genes)+1)


''' Update the 'proteins' table to include the gene_id as future secondary key '''
logger.info("Optimizing Proteins table.")
# TODO: Check if the indices still need to be created or if they are already defined
# Make 'genes' index based on original gene id and taxonomy id
genes.set_index(['taxonomy_id','original_gene_id'], drop=False, inplace=True)

# Make another 'proteins' index based on original gene id and taxonomy id to help a merge later on
proteins.set_index(['taxonomy_id','original_gene_id'], drop=False, inplace=True)

# Update 'proteins' table to include the gene_id
proteins = proteins.merge(genes[['gene_id']],left_index=True,right_index=True)
# and remove the original_gene_id as we no longer need it.
del proteins['original_gene_id']
# Sort the protein table, because it got messed up in the merge
proteins.sort_values(by='protein_id', ascending=True, inplace=True)

# Drop tables if exists
logger.info("Dropping tables if exists.")
try:
    Proteins.__table__.drop(checkfirst=True)
    Genes.__table__.drop(checkfirst=True)
    Species.__table__.drop(checkfirst=True)
#    Eukarya.metadata.drop_all()
except:
    logger.error("Unable to drop the tables.")
    exit(1)
# Create all tables
logger.info("Creating tables.")
try:
    Species.__table__.create()
    Genes.__table__.create()
    Proteins.__table__.create()
except:
    logger.error("Unable to create tables.")
    exit(1)

# Use the Pandas.to_sql() function to populate the tables
# Correct order is Species, Genes, Proteins
# This is due to foreign keys

# Store species table
logger.info("Storing Species data.")
try:
    species[['taxonomy_id', 'abbreviation', 'supergroup',
       'relevant_taxonomy', 'scientific_name', 'common_name', 'core',
       'download_location', 'download_date_yyyymmdd', 'assembly_version',
       'reference', 'notes', 'fasta_header_type',
       'proteome_fasta_file','genome_fasta_file']].to_sql(name=Species.__table__.name, schema="eukarya", con=engine,
                                                          if_exists='append', index=False)
except Exception as e:
    logger.error(e)
    exit(1)
# Store genes table
logger.info("Storing Genes Table.")
genes[['gene_id','gene_symbol',
       'taxonomy_id','original_gene_id',
       'protein_id']].to_sql(name=Genes.__table__.name, schema="eukarya", con=engine,
                                                      if_exists='append', index=False)

# Store proteins table
logger.info("Storing Proteins Table.")
proteins[['protein_id','original_protein_id',
          'original_transcript_id','gene_id',
          'taxonomy_id','sequence_length','md5_sum',
          'longest_transcript','original_header']].to_sql(name=Proteins.__table__.name, schema="eukarya", con=engine,
                                                          if_exists='append', index=False)
# Let SQLalchemy check if foreign keys in various tables match!
