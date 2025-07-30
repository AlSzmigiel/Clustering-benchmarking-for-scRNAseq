#!/bin/bash
OUTPUT_DIR=$1
mkdir -p "$OUTPUT_DIR"
wget -O "$OUTPUT_DIR/kolodziejczyk.csv" \
'https://espresso.teichlab.sanger.ac.uk/static/counttable_es.csv'

