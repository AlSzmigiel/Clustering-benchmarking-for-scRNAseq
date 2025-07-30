#!/bin/bash

# Usage: ./script.sh <output_directory>
OUTPUT_DIR=$1

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Download data
wget -O "$OUTPUT_DIR/Biase_fpkm.txt.gz" \
'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE57249&format=file&file=GSE57249%5Ffpkm%2Etxt%2Egz'

# Unzip the downloaded file
gunzip "$OUTPUT_DIR/Biase_fpkm.txt.gz"