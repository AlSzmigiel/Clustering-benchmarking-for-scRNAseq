#!/bin/bash


OUTPUT_DIR=$1


mkdir -p "$OUTPUT_DIR"


# URL of the raw files on GitHub
sin_cell='https://github.com/LuyiTian/sc_mixology/raw/master/data/sincell_with_class.RData'
sin_cell_5='https://github.com/LuyiTian/sc_mixology/raw/master/data/sincell_with_class_5cl.RData'
output_dir_sin_cell="$OUTPUT_DIR/sincell_with_class.RData"
output_dir_sin_cell_5cl="$OUTPUT_DIR/sincell_with_class_5cl.RData"


# Download the first file
wget -O $output_dir_sin_cell  $sin_cell

# Download the second file
wget -O $output_dir_sin_cell_5cl  $sin_cell_5