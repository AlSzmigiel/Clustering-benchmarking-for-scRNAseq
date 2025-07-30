import sys
import os
import pandas as pd
import anndata
import scanpy as sc
import anndata as ad


def check_file_exists(file_path):
    """Helper function to check if a file exists."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")

def main(input_dir, output_file):
    cell_types=[ "cd56_nk", "cd4_t_helper" ,"regulatory_t", "memory_t", "naive_t", "naive_cytotoxic", "cytotoxic_t" , "b_cells","cd34"]
    #Check if paths to each cell_types exists
    mat_path = 'filtered_matrices_mex/hg19/'
    adatasets = {}
    for c_t in cell_types:
        counts_path = os.path.join(input_dir,c_t,mat_path ,'matrix.mtx')  # Assuming the file is named `data.txt`
        print(counts_path)
        genes_path = os.path.join(input_dir, c_t, mat_path, 'genes.tsv')
        barcodes_path = os.path.join(input_dir, c_t, mat_path, 'barcodes.tsv')   
        check_file_exists(counts_path)
        check_file_exists(genes_path)
        check_file_exists(barcodes_path)    
        adata = sc.read_mtx(counts_path).T  # transpose the data
        adata.var_names = pd.read_csv(genes_path, header=None, sep='\t')[1]
        adata.obs_names = pd.read_csv(barcodes_path, header=None)[0]
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
        adatasets[c_t]=adata
    PBMCs = ad.concat(adatasets, label="label", join='inner', index_unique='_')
    print(f"Created AnnData object: {PBMCs}. Saving processed data to: {output_file}")
    PBMCs.write(output_file)

if __name__ == "__main__":
    # Command-line arguments: input_dir and output_file
    if len(sys.argv) != 4:
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]

    main(input_dir, output_file)