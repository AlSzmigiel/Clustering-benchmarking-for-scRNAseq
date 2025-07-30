import scanpy as sc
from scipy.sparse import coo_matrix, csr_matrix, issparse
import pandas as pd 
#import compositional as coda
import rpy2.robjects as ro
from rpy2.robjects import r
import rpy2.robjects.numpy2ri
import numpy as np
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr, data


def calculate_graph(data, nn):
    distances = pd.DataFrame(data.obsp['distances_full'])
    with (ro.default_converter + pandas2ri.converter).context():
        r_from_pd_df = ro.conversion.get_conversion().py2rpy(distances)
    seurat = importr("Seurat")
    mat = seurat.FindNeighbors_dist(r_from_pd_df, k_param = nn)               #Seurat- phenograph 
    data.obsp['connectivities'] = csr_matrix(np.array(r['as.matrix'](mat[1])))
    data.obsp['distances'] = csr_matrix(np.array(r['as.matrix'](mat[0])))        
    data.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"method": None}}
    return data 

graph_mat = calculate_graph(sc.read_h5ad(str(snakemake.input)), int(snakemake.wildcards.neighbors)) 
graph_mat.write(str(snakemake.output))

#Calculating umap weights





# else:
#     if snakemake.wildcards.dist_metric=='spearman':
#         mat = sp.stats.spearmanr(adat.X, axis=1)
#         distance= 1 - mat.correlation


#     if (snakemake.wildcards.graph_method=='jaccard'):
#         #sc.pp.log1p(adat)
#         distance=pairwise_distances(adat.X,metric=snakemake.wildcards.dist_metric)
#         conn=kneighbors_graph(distance, n_neighbors=int(snakemake.wildcards.neighbors)  ,mode='distance', metric='precomputed', metric_params=None, include_self=False)
#         graph=sc.external.tl.phenograph(conn,directed=False,clustering_algo=None)
#         adat.obsp["distances"]=distance
#         adat.obsp["connectivities"]=graph[1].tocsr()
#         adat.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"method": None}}
#         graph_ds=adat
        
#     else:
#         #sc.pp.log1p(adat)
#         graph_ds=sc.pp.neighbors(adat,knn=('True'==str(snakemake.wildcards.nns)),method=snakemake.wildcards.graph_method,n_neighbors=int(snakemake.wildcards.neighbors),metric=snakemake.wildcards.dist_metric,use_rep='X',copy=True)




#boolean is read in some different ways by snakemake so the comparison of string has to be made 
