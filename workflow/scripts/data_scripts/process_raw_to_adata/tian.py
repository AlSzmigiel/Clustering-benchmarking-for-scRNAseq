import sys
import os
import pandas as pd
import anndata
import scanpy as sc
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
import anndata2ri
import conversion
#import snakemake

def main(input_dir, output_file, subdata):
    # Step 1: Read the input .txt file
    input_file = os.path.join(input_dir)  # Assuming the file is named `data.txt`
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    adatas = conversion.conversion(input_file)
    dataset = adatas[subdata]
    dataset.uns["name"] = subdata
    dataset.obs.index=dataset.obs.index.astype(str).astype("category") #because .write dont work with python type 'string' but works with 'str'
    dataset.var.index=dataset.var.index.astype(str).astype("category")
    dataset.uns['scPipe']=None #because that was an R object (after conversion) and couldnt be used with .write 
    dataset.uns['Biomart']=None
    dataset.write(os.path.join(output_file))
    print(f"Created AnnData object: {dataset}. Saving processed data to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    subdata = sys.argv[3]

    main(input_dir, output_file, subdata)