import sys
import os
import pandas as pd
import anndata
import scanpy as sc


def main(input_dir, output_file):
    # Step 1: Read the input .txt file
    input_file = os.path.join(input_dir, "Biase_fpkm.txt")  # Assuming the file is named `data.txt`
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    adata = sc.read_text(input_file).T  # Tab-separated, with row names
    adata = adata[:49,:].copy()
    cell_type = ['zygote']*9+['two_cell']*20+['four_cell']*20 
    adata.obs['cell_type'] = cell_type
    # Step 3: Write the AnnData object to .h5ad format
    print(f"Created AnnData object: {adata}. Saving processed data to: {output_file}")
    adata.write(output_file)

if __name__ == "__main__":
    # Command-line arguments: input_dir and output_file
    if len(sys.argv) != 4:
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]

    main(input_dir, output_file)