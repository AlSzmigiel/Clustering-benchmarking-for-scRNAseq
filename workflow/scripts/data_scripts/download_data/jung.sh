#!/bin/bash


output_directory=$1

mkdir -p "$output_directory"  #create if does not exist

# Array of URLs
declare -a urls=(
  'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3102985&format=file&file=GSM3102985%5FSPDmat%2Etxt%2Egz'
  'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3102982&format=file&file=GSM3102982%5FPouf5SPGmat%2Etxt%2Egz'
  'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3102983&format=file&file=GSM3102983%5FSPCImat%2Etxt%2Egz'
  'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3102984&format=file&file=GSM3102984%5FSPCIImat%2Etxt%2Egz'
)

declare -a output_names=(
  'spermatid.txt.gz'
  'spermatogonia.txt.gz'
  'primary_spermatocyte.txt.gz'
  'secondary_spermatocyte.txt.gz'
)

for i in "${!urls[@]}"; do
  wget -O "$output_directory/${output_names[i]}" "${urls[i]}"
  gunzip -f "$output_directory/${output_names[i]}"
  done