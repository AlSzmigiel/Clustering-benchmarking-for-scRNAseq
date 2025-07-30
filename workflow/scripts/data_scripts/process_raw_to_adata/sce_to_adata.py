import sys
import os
import pandas as pd
import anndata
import scanpy as sc
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
import anndata2ri
import conversion

def main(input_file, output_file):
    adata = conversion.conversion(input_file)
    # dataset = adatas[subdata]
    # dataset.uns["name"] = subdata
    adata.obs.index = adata.obs.index.astype(str).astype("category") #because .write dont work with python type 'string' but works with 'str'
    adata.var.index = adata.var.index.astype(str).astype("category")
    # dataset.uns['scPipe']=None #because that was an R object (after conversion) and couldnt be used with .write 
    # dataset.uns['Biomart']=None
    adata.write(os.path.join(output_file))
    print(f"Created AnnData object: {adata}. Saving processed data to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    input_file = sys.argv[3]

    main(input_file, output_file)