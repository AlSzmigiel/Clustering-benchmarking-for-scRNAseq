import os 
import pandas as pd
import anndata as ad
import sys

def main(input_dir, output_file):
    # Step 1: Read the input .txt file
    input_file = os.path.join(input_dir, "Koh.txt")  # Assuming the file is named `data.txt`
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")    
    # Split column names into parts
    counts = pd.read_csv(input_file, sep='\t')
    cell_types = [element.split('.')[0] for element in list(counts.columns.values)[2:]]
    obs = pd.DataFrame()
    var = pd.DataFrame()
    var['genes'] = counts.geneID.values.astype('str')
    obs['sample'] = counts.columns.values[2:]
    obs['cell_labels'] = cell_types
    adata = ad.AnnData(X = counts.iloc[:,2:].values.T,
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