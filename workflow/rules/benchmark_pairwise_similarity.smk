rule pairwise_distances:
    input: rules.preprocess_normalize_adata.output
    output: "results/benchmark_analysis/pairwise_similarity/{dataset}/{sim_measure}.h5ad"
    conda: "../envs/r_tools_graph.yaml"
    script:
        "../scripts/pairwise_similarity.py"