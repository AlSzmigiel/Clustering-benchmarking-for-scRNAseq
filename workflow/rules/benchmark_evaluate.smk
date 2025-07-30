rule evaluate_partitions:
    input: rules.clustering.output 
    output: "results/benchmark_analysis/evaluation/{dataset}/{clust_method}/{graph_method}/{sim_measure}/{neighbors}k.summary.txt"  #it should be possible to do it per dataset 
    params:
            partition = lambda wildcards: cfg.get_clustering_params(wildcards.clust_method, key='partition_type'),
            label_key = lambda wildcards: cfg.get_from_dataset(wildcards.dataset,key='cell_labels')        
    conda:  "../envs/benchmark_evaluation.yaml"                
    script: "../scripts/evaluate_part.py"