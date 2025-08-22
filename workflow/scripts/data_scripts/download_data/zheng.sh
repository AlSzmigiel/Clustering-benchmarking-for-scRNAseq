#!/bin/bash

OUTPUT_DIR=$1
mkdir -p "$OUTPUT_DIR"
declare -a arr=("cd14_monocytes" "cd56_nk" "cd4_t_helper" "regulatory_t" "memory_t" "naive_t" "naive_cytotoxic" "cytotoxic_t" "cd34" "b_cells")

for i in "${arr[@]}"
do
    wget -P $OUTPUT_DIR/$i \
    "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/$i/${i}_filtered_gene_bc_matrices.tar.gz"
    tar -xvf $OUTPUT_DIR/$i/${i}_filtered_gene_bc_matrices.tar.gz -C $OUTPUT_DIR/$i
done



