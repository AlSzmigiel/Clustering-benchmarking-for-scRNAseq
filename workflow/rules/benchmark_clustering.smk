rule clustering:
    input: "results/benchmark_analysis/graphs/{dataset}/{graph_method}/{sim_measure}/{neighbors}k.h5ad"
    output: "results/benchmark_analysis/clustering/{dataset}/{clust_method}/{graph_method}/{sim_measure}/{neighbors}k.h5ad"
#     output: cluster_data="analysis/clustering/{dataset}/{rep}/{method}/{k}n_{nn}/{clust_method}_{res}.h5ad"
    params:
        res = lambda wildcards: cfg.get_clustering_params(wildcards.clust_method, key='res'),
        partition= lambda wildcards: cfg.get_clustering_params(wildcards.clust_method, key='partition_type')
    conda:
            "../envs/benchmark_clustering.yaml"
    script:"../scripts/clustering.py"