import sys
import os
import pandas as pd
import anndata
import scanpy as sc
import anndata as ad

def main(input_dir, output_file):
    input_file = os.path.join(input_dir, "41594_2013_BFnsmb2660_MOESM31_ESM.xlsx")  # Assuming the file is named `data.txt`
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    counts = pd.read_excel(input_file)
    # chose all without hESC passage (hemberglab)
    expr = counts.iloc[:,2:92]
    cell_labels = [name.split()[0] for name in expr.columns.values]
    obs = pd.DataFrame()
    var = pd.DataFrame()
    var['genes'] = counts.Gene_ID.values.astype('str')
    obs['sample'] = expr.columns.values
    obs['cell_labels'] = cell_labels
    adata = ad.AnnData(X=expr.values.T,
                                obs=obs,
                                var=var)
    print(f"Created AnnData object: {adata}. Saving processed data to: {output_file}")
    adata.write(output_file)

if __name__ == "__main__":
    # Command-line arguments: input_dir and output_file
    if len(sys.argv) != 4:
        sys.exit(1)
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    main(input_dir, output_file)