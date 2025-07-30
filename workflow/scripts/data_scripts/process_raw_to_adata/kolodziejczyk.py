import sys
import os
import pandas as pd
import anndata
import scanpy as sc
import anndata as ad


def main(input_dir, output_file):
    # Step 1: Read the input .txt file
    input_file = os.path.join(input_dir, "kolodziejczyk.csv")  # Assuming the file is named `data.txt`
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")    
    # Split column names into parts
    counts = pd.read_csv(input_file , delim_whitespace=True)
    parts = [col.split("_") for col in counts.columns]
    # Extract the cell type (third element)
    cell_type = [part[2] for part in parts]
    # Extract the batch (concatenate third and fourth elements with an underscore)
    batch = ["_".join(part[2:4]) for part in parts]
    obs = pd.DataFrame()
    var = pd.DataFrame()
    var['genes']= counts.index.values.astype('str')[:-5]
    obs['sample'] = counts.columns.values
    obs['batch'] = batch
    obs['cell_labels']= cell_type

    adata = ad.AnnData(X = counts.values[:-5,:].T,
                                obs=obs ,
                                var=var)
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