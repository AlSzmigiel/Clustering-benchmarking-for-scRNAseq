rule generate_graph_scanpy:
    input: rules.pairwise_distances.output
    output: "results/benchmark_analysis/graphs/{dataset}/{graph_method}/{sim_measure}/{neighbors}k.h5ad"
    wildcard_constraints:
        graph_method = "gauss|umap|jaccard_phenograph"
    conda: "../envs/scanpy_graphs.yaml"
    script: "../scripts/construct_graph_scanpy.py"

rule generate_graph_seurat:
    input: rules.pairwise_distances.output
    output: "results/benchmark_analysis/graphs/{dataset}/{graph_method}/{sim_measure}/{neighbors}k.h5ad"
    wildcard_constraints:
        graph_method="jaccard_seurat"  # only runs when graph_method == "seurat"
    conda: "../envs/seurat_graph.yaml"
    script: "../scripts/construct_graph_seurat.py"
