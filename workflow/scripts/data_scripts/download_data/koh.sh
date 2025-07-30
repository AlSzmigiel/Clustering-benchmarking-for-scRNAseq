#!/bin/bash

OUTPUT_DIR=$1

wget -O $OUTPUT_DIR/Koh.txt.gz \
'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2257302&format=file&file=GSM2257302%5FAll%5Fsamples%5Fsc%5Ftpm%2Etxt%2Egz'
gunzip $OUTPUT_DIR/Koh.txt.gz


