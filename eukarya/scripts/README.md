# Eukarya scripts collection

This repo holds generic scripts that help in working
with the Eukarya set. It also contains the `eukarya` python package
for use in your own scripts and takes care of some standard operations,
such as database access, file locations, etc.

This repo should not be cloned directly. Instead clone the *base* repo
to pull everything in.

Below a description of each script.

## build_eukarya_db.py
Builds a SQLite database from the data in the Eukarya dataset for ease of use.
Run `python3 build_eukarya_db.py -h to get a full description of options.`

## build_orthofinder_table.py
Builds a SQLite table in the annotations database for a particular orthology

## found_before.py

## get_orthofinder_og_sequences.py

## get_sequence.py

## reformat_orthofinder_ogs.py

## reformat_protein_ids.py
