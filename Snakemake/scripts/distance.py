import scanpy as sc
from scipy.sparse import coo_matrix, csr_matrix, issparse
from scanpy.neighbors._common import _get_indices_distances_from_dense_matrix
from scanpy.neighbors._common import _get_indices_distances_from_sparse_matrix
import pandas as pd 
import scipy as sp
#import compositional as coda
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
from rpy2.robjects.packages import importr, data
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import kneighbors_graph
import scipy.sparse as spa

rpy2.robjects.numpy2ri.activate()




def calculate_r_measure(data, measure):
    propr = importr('propr')
    base = importr('base')
    at = base.__dict__["@"]
    count_matrices = ro.r.matrix(data, nrow = data.shape[0], ncol = data.shape[1])
    ro.r.assign("mat", count_matrices) 
    dist_obj= propr.propr(ro.r.mat.transpose(), metric = measure)
    distance = at(dist_obj, "matrix")
    if measure == 'rho':
        return 1 - distance
    return distance
    


def calculate_distance_matrix(data, measure):
    if (measure == 'phs' or measure == 'rho'):
        data.obsp['distances_full'] = calculate_r_measure(data.X, measure)
    elif (measure == 'spearman'):
        corr = sp.stats.spearmanr(data.X, axis=1)
        data.obsp['distances_full'] = 1 - corr.correlation
    else:
        data.obsp['distances_full'] = pairwise_distances(data.X, metric = measure)
    return data


distance_matrix = calculate_distance_matrix(sc.read_h5ad(str(snakemake.input)), snakemake.wildcards.sim_measure)

distance_matrix.write(str(snakemake.output))