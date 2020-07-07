#! /usr/bin/env python3
'''
    Author: John van Dam
    Created: Septermber 1st, 2017
    Version: 0.1.0

    Purpose of module: collection of classes and functions for use in eukarya
    database related projects.
'''

import sys
import os
import unittest
import logging
import pandas as pd
import numpy as np
from sqlalchemy import Column, ForeignKey, Integer, String, Boolean, create_engine, MetaData
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
import eukarya.files as files

logger = logging.getLogger(__name__)  # Set up the logger

# Setting up the database engines
try:
    engine  = create_engine('sqlite://',echo=False)  # generate an "in mem" database to attach multiple sqlite databases to.
    engine.execute("attach database '" + files.eukarya_dir + files.eukarya_file['sqlite3'] + "' as eukarya;")
    engine.execute("attach database '" + files.eukarya_dir + files.annotations_file['sqlite3'] + "' as annotations;")
    logger.debug("Successfully set Database engine.")
except Exception as e:
    logger.error("Failed to set Database engine. {0}".format(e))
    exit(1)

# Setting up the session managers
try:
    Session = sessionmaker(bind=engine)
except Exception:
    logger.error("Failed to create Session makers.")
    exit(1)


# Declaring SQLalchemy objects
""" This base SQLalchemy data object will hold all metadata related to the
Eukarya database and Annotations database for use in querying. All subsequent objects (tables),
need to be derived from this base class.
"""
Base = declarative_base(bind=engine)

# The base eukarya data set tables
class Eukarya():
    '''
    Base class for the Eukarya database.
    It presets the schema for all child classes and it makes it easier to
    obtain the tables by using the __subclasses__() routine.
    '''
    __table_args__ = {'schema': 'eukarya'}


class Species(Base,Eukarya):
    """SQLalchemy object representing the Species table."""
    __tablename__ = 'species'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    taxonomy_id = Column(Integer, unique=True)#, index=True)
    abbreviation = Column(String(4))
    supergroup = Column(String)
    relevant_taxonomy = Column(String)
    scientific_name = Column(String)
    common_name = Column(String)
    core = Column(Boolean)
    download_location = Column(String)
    download_date_yyyymmdd = Column(String)
    assembly_version = Column(String)
    reference = Column(String)
    notes = Column(String)
    fasta_header_type = Column(String)
    proteome_fasta_file = Column(String)
    genome_fasta_file = Column(String)


class Genes(Base,Eukarya):
    """SQLalchemy object representing the Genes table."""
    __tablename__ = 'genes'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, unique=True)#, index=True)
    gene_symbol = Column(String(16))#, index=True)
    taxonomy_id = Column(Integer, ForeignKey(Species.taxonomy_id))#, index=True)
    original_gene_id = Column(String)
    # I am dropping this column's (protein_id) foreign key as SQLalchemy is
    # having trouble dropping tables with circular foreign key constraints.
    # In fact theoretically we don't need this column but can be mimicked by
    # joining with the protein table, filtering for longest_transcript == 1.
    # I am keeping it, purely as a query short cut.
    protein_id = Column(String(10), unique=True)#, index=True)


class Proteins(Base,Eukarya):
    """SQLalchemy object representing the Proteins table."""
    __tablename__ = 'proteins'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    protein_id = Column(String(10), unique=True)#, index=True)
    original_protein_id = Column(String())
    original_transcript_id = Column(String)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)
    taxonomy_id = Column(Integer, ForeignKey(Species.taxonomy_id))#, index=True)
    sequence_length = Column(Integer)
    md5_sum = Column(String)
    longest_transcript = Column(Boolean)
    original_header = Column(String)


# Orthology database tables
class Annotations():
    '''
    Base class for annotation data tables.
    It presets the schema for all child classes and it makes it easier to
    obtain the tables by using the __subclasses__() routine.
    '''
    __table_args__ = {'schema': 'annotations'}

class Orthology(Annotations):
    '''
    Base class for orthologies. Inherits from Annotations().
    It presets the schema for all child classes and it makes it easier to
    obtain the tables by using the __subclasses__() routine.
    '''

class OrthofinderDiamond(Base,Orthology):
    """SQLalchemy object representing the Orthofinder table based on Diamond."""
    __tablename__ = 'orthofinder_diamond'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class OrthofinderBlast(Base,Orthology):
    """SQLalchemy object representing the Orthofinder table based on Blast with a cut-ff set to 10."""
    __tablename__ = 'orthofinder_blast_10'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class OrthofinderBlast_1(Base,Orthology):
    """SQLalchemy object representing the Orthofinder table based on Blast with cut-off set to e-1."""
    __tablename__ = 'orthofinder_blast_1'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class OrthofinderBlast_3(Base,Orthology):
    """SQLalchemy object representing the Orthofinder table based on Blast with cut-off set to e-3."""
    __tablename__ = 'orthofinder_blast_3'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class EggnogDiamond(Base,Orthology):
    """SQLalchemy object representing the Eggnog table."""
    __tablename__ = 'eggnog_diamond'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class EggnogHmmer(Base,Orthology):
    """SQLalchemy object representing the Eggnog table."""
    __tablename__ = 'eggnog_hmmer'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class EggnogHmmerCor(Base,Orthology):
    """SQLalchemy object representing the Eggnog corrected table."""
    __tablename__ = 'eggnog_hmmer_corrected'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class Julian(Base,Orthology):
    """SQLalchemy object representing the Julian (all) table."""
    __tablename__ = 'julian'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)


class Julian90(Base,Orthology):
    """SQLalchemy object representing the reduced Julian (all) table."""
    __tablename__ = 'julian90'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)


class Julian50(Base,Orthology):
    """SQLalchemy object representing the further reduced Julian (all) table."""
    __tablename__ = 'julian50'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)


class Panther(Base,Orthology):
    """SQLalchemy object representing the Panther table."""
    __tablename__ = 'panther'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class Manual(Base,Orthology):
    """SQLalchemy object representing the Manual OGs table."""
    __tablename__ = 'manual'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class Broccoli(Base,Orthology):
    """SQLalchemy object representing the Broccoli OGs table."""
    __tablename__ = 'broccoli'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class Panther_cor(Base,Orthology):
    """SQLalchemy object representing the panther panthercorrected OGs table."""
    __tablename__ = 'panther_corrected'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class Swiftortho(Base,Orthology):
    """SQLalchemy object representing the Swiftortho OGs table."""
    __tablename__ = 'swiftortho'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class Sonicparanoid_default(Base,Orthology):
    """SQLalchemy object representing the Sonicparanoid_default OGs table."""
    __tablename__ = 'sonicparanoid_default'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class Sonicparanoid_fast(Base,Orthology):
    """SQLalchemy object representing the Sonicparanoid_fast OGs table."""
    __tablename__ = 'sonicparanoid_fast'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

class Sonicparanoid_sensitive(Base,Orthology):
    """SQLalchemy object representing the Sonicparanoid_sensitive OGs table."""
    __tablename__ = 'sonicparanoid_sensitive'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    gene_id = Column(Integer, ForeignKey(Genes.gene_id))#, index=True)#, unique=True)
    og_id = Column(String(10))#, index=True)

# LECA OG lists per orthology
class OrthologyLECA(Annotations):
    '''
    Base class for LECA OG lists. Inherits from Annotations().
    It presets the schema for all child classes and it makes it easier to
    obtain the tables by using the __subclasses__() routine.
    '''
class OrthofinderDiamondLECA(Base,Orthology):
    """SQLalchemy object representing the Orthofinder LECA OG list based on Diamond."""
    __tablename__ = 'orthofinder_diamond_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(OrthofinderDiamond.og_id))#, index=True)

class OrthofinderBlastLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Orthofinder LECA OG list based on Blast with a cut-ff set to 10."""
    __tablename__ = 'orthofinder_blast_10_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(OrthofinderBlast.og_id))#, index=True)

class OrthofinderBlast_1LECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Orthofinder LECA OG list based on Blast with cut-off set to e-1."""
    __tablename__ = 'orthofinder_blast_1_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(OrthofinderBlast_1.og_id))#, index=True)

class OrthofinderBlast_3LECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Orthofinder LECA OG list based on Blast with cut-off set to e-3."""
    __tablename__ = 'orthofinder_blast_3_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(OrthofinderBlast_3.og_id))#, index=True)

class EggnogDiamondLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Eggnog LECA OG list."""
    __tablename__ = 'eggnog_diamond_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(EggnogDiamond.og_id))#, index=True)

class EggnogHmmerLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Eggnog LECA OG list."""
    __tablename__ = 'eggnog_hmmer_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(EggnogHmmer.og_id))#, index=True)

class EggnogHmmerLECACor(Base,OrthologyLECA):
    """SQLalchemy object representing the Eggnog corrected LECA OG list."""
    __tablename__ = 'eggnog_hmmer_corrected_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(EggnogHmmerCor.og_id))#, index=True)

class JulianLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Julian (all) LECA OG list."""
    __tablename__ = 'julian_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Julian.og_id))#, index=True)

class Julian90LECA(Base,OrthologyLECA):
    """SQLalchemy object representing the reduced Julian (all) LECA OG list."""
    __tablename__ = 'julian90_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Julian90.og_id))#, index=True)

class Julian50LECA(Base,OrthologyLECA):
    """SQLalchemy object representing the further reduced Julian (all) LECA OG list."""
    __tablename__ = 'julian50_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Julian50.og_id))#, index=True)

class PantherLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Panther LECA OG list."""
    __tablename__ = 'panther_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Panther.og_id))#, index=True)

class ManualLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Manual LECA OG list."""
    __tablename__ = 'manual_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Manual.og_id))

class BroccoliLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Broccoli LECA OG list."""
    __tablename__ = 'broccoli_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Broccoli.og_id))

class PantherLECA_Cor(Base,OrthologyLECA):
    """SQLalchemy object representing the Panther corrected LECA OG list."""
    __tablename__ = 'panther_corrected_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Panther_cor.og_id))#, index=True)

class SwiftorthoLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Swiftortho LECA OG list."""
    __tablename__ = 'swiftortho_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Swiftortho.og_id))#, index=True)

class Sonicparanoid_defaultLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Sonicparanoid_default LECA OG list."""
    __tablename__ = 'sonicparanoid_default_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Sonicparanoid_default.og_id))#, index=True)

class Sonicparanoid_fastLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Sonicparanoid_fast LECA OG list."""
    __tablename__ = 'sonicparanoid_fast_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Sonicparanoid_fast.og_id))#, index=True)

class Sonicparanoid_sensitiveLECA(Base,OrthologyLECA):
    """SQLalchemy object representing the Sonicparanoid_sensitive LECA OG list."""
    __tablename__ = 'sonicparanoid_sensitive_leca_list'
    id = Column(Integer, primary_key=True, unique=True, autoincrement=True)
    og_id = Column(String(10), ForeignKey(Sonicparanoid_sensitive.og_id))#, index=True)
# Legacy
Eggnog = EggnogDiamond  # Legacy object name

def get_orthology_leca_tables():
    '''
    Returns a dict with tubles containing the main orthology table and the
    table with LECA ogs.
    '''
    # Define the Orthology / leca_list pairs as a dict of tuples
    orthology_leca_tables = dict()
    orthology_leca_tables['orthofinder_diamond_e-3'] = (OrthofinderDiamond, OrthofinderDiamondLECA)
    orthology_leca_tables['eggnog_diamond'] = (EggnogDiamond, EggnogDiamondLECA)
    orthology_leca_tables['eggnog_hmmer_corrected'] = (EggnogHmmerCor, EggnogHmmerLECACor)
    orthology_leca_tables['orthofinder_blast_e-3'] = (OrthofinderBlast_3, OrthofinderBlast_3LECA)
    orthology_leca_tables['panther_corrected'] = (Panther_cor, PantherLECA_Cor)
    orthology_leca_tables['manual'] = (Manual, ManualLECA)
    orthology_leca_tables['broccoli'] = (Broccoli, BroccoliLECA)
    orthology_leca_tables['swiftortho'] = (Swiftortho, SwiftorthoLECA)
    orthology_leca_tables['sonicparanoid_sensitive'] = (Sonicparanoid_sensitive, Sonicparanoid_sensitiveLECA)
    return orthology_leca_tables

def get_orthology_tables():
    '''
    Returns a dict containing the orthology tables derived from
    get_orthology_leca_tables()
    '''
    orthologies_leca = get_orthology_leca_tables()
    return {k:v[0] for k,v in orthologies_leca.items()}


class TestConnection(unittest.TestCase):
    ''' Unit tests to test database connection.'''
    def test_session(self):
        ''' Tests a session and performs a query with three tables. '''
        try:
            session = Session()
            session.query(Genes).join(Proteins).join(Species).limit(10)
        except:
            self.assertTrue(False)
        else:
            self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
