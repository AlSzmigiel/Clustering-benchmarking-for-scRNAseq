import pandas as pd 
import anndata2ri
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
import os

def conversion(path):   #specifies path to folder with SingleCellExperiment R objects to be converted 
    anndata = []
    exp_names = []
    
    # Check if path is a file or a directory
    if os.path.isfile(path):
        file_list = [path] if path.endswith(".RData") else []
    elif os.path.isdir(path):
        file_list = [os.path.join(path, filename) for filename in os.listdir(path) if filename.endswith(".RData")]
    else:
        raise ValueError("Provided path is neither a valid file nor a directory.")
    for Rdat_path in file_list:
        sceexp=robjects.r['load'](Rdat_path)
        for i in sceexp:  #can be one or more objects
            if (tuple(robjects.r[i].rclass)[0]=='SingleCellExperiment'):
                print('{name} is a SingleCellExperiment object, converting to anndata...'.format(name=i))
                exp_names.append(i)
                dataset=robjects.r[i]
                with localconverter(anndata2ri.converter):
                    anndata.append(robjects.r('as')(dataset, 'SingleCellExperiment'))
            else:
                print(print('{name} is a {obj_class} object'.format(name=i,obj_class=tuple(robjects.r[i].rclass)[0])))
    if len(anndata) == 1:
        return anndata[0]
    datasets = dict(map(lambda i,j : (i,j) , exp_names, anndata))
    return(datasets)



