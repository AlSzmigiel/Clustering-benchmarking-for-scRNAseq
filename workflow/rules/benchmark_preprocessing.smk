rule preprocess_normalize_adata:
    input:
        processed_adata=lambda wildcards: (
            "data/processed_adata/{dataset}.h5ad"
            if cfg.get_from_dataset(wildcards.dataset, key="standard") == 'gold'
            else "data/enhanced_silver_standard/{dataset}-new_lab.h5ad"
        )
    output:
        processing="results/benchmark_analysis/processed_normalized_data/{dataset}.h5ad"
    params:
        input_rep=lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key='data_type'),
        batch_key=lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key='batch_labels'),
        cell_label=lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key='cell_labels'),
        visualization="results/benchmark_analysis/visualization/batch_removal/{dataset}.png"
    conda: "../envs/benchmark_process_normalize.yaml"
    script: "../scripts/preprocess_normalize.py"