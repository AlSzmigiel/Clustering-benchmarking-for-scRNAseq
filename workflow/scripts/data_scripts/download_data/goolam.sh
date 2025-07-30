#!/bin/bash

# Define the output directory and create it if it doesn't exist
OUTPUT_DIR=$1
mkdir -p "$OUTPUT_DIR"

# Download data files
wget -O "${OUTPUT_DIR}/goolam_counts.tsv" \
'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/321/E-MTAB-3321/Files/Goolam_et_al_2015_count_table.tsv'

wget -O "${OUTPUT_DIR}/goolam_metadata.txt" \
'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/321/E-MTAB-3321/Files/E-MTAB-3321.sdrf.txt'