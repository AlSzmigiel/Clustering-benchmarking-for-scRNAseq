import scanpy as sc
import leidenalg as la 
import louvain as lv
from scanpy._utils import get_igraph_from_adjacency
import numpy as np
import pandas as pd
from natsort import natsorted


# LOUVAIN IN SCANPY DOES NOT IMPLEMENT POSSIBILITY TO RUN FOR MULTIPLE ITERATIONS
# You can't pass n_iterations argument to louvain 



def cluster_data(data, algorithm, partition_list, resolution_list):
    G = get_igraph_from_adjacency(data.obsp['connectivities'],directed=False)
    weights =np.array(G.es["weight"]).astype(np.float64)
    for part in partition_list:
        for res in resolution_list:
            if (algorithm=='leiden'):
                sc.tl.leiden(data,resolution=res,key_added=f"{part}_{res}",partition_type=getattr(la,part),use_weights=True,directed=None, n_iterations = 10, random_state =301)
            else:
                optimiser = lv.Optimiser()
                optimiser.set_rng_seed(0)
                get_part = getattr(lv,part)
                partition = get_part(G, resolution_parameter = res, weights = weights)
                for i in range (0,10):
                     optimiser.optimise_partition(partition)
                groups = np.array(partition.membership)
                data.obs[f"{part}_{res}"] = pd.Categorical(
                values=groups.astype("U"),
                categories=natsorted(map(str, np.unique(groups))),
                )
            print(algorithm, part, res)
    del data.obsp['distances']
    del data.obsp['distances_full']
    del data.obsp['connectivities']
    return data


clustered = cluster_data(sc.read(str(snakemake.input)),snakemake.wildcards.clust_method,snakemake.params.partition,snakemake.params.res)
clustered.write(str(snakemake.output))

