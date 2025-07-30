import scanpy as sc
import anndata as ad
import numpy as np
import os
import sys

def main(input_dir, output_file):
    adatasets={}
    cell_types=['primary_spermatocyte','secondary_spermatocyte','spermatid','spermatogonia']
    for c_t in cell_types:
        input_file = os.path.join(input_dir, c_t +'.txt')  # Assuming the file is named `data.txt`
        if not os.path.exists(input_file):
                raise FileNotFoundError(f"Input file not found: {input_file}") 
        adata = sc.read_text(input_file).T  # transpose the data
        adatasets[c_t]=adata
    ds_all = ad.concat(adatasets,label="label",merge = 'same', index_unique='_')
    print(f"Created AnnData object: {ds_all}. Saving processed data to: {output_file}")
    ds_all.write(output_file)

if __name__ == "__main__":
    # Command-line arguments: input_dir and output_file
    if len(sys.argv) != 4:
        sys.exit(1)
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    main(input_dir, output_file)