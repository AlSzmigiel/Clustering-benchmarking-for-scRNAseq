rule download_cellxgene:
    output:
        processed_file="data/processed_adata/{dataset}.h5ad"
    params:
        #download_link=lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key="download_script"),
        dataset_id = lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key="dataset_id")
    wildcard_constraints:
        dataset="[a-zA-Z0-9]+_cxg" 
    conda:
        "../envs/cellxgene.yaml"
    script:
        "../scripts/get_cellxgene_data.py"

    # shell:
    #     """
    #     wget -O {output.processed_file} {params.download_link}
    #     """

rule download_bash:
    output:
         touch("data/raw_data/{dataset}/.done")
    params:
        script=lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key="download_script"),
        output_dir= lambda wildcards:  f"data/raw_data/{wildcards.dataset}"
    shell:
        """
        bash {params.script} {params.output_dir}
        """
rule download_r:
    output:
        output = "data/raw_data/{dataset}.RData"
    params:
        script = lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key="download_script"),
    conda:
        "../envs/r_download.yaml"
    shell:
        """
        Rscript {params.script} {output} {wildcards.dataset}
        """


rule process_to_adata:
    input:
        raw_file=lambda wildcards: "data/raw_data/{dataset}/.done" 
            if cfg.get_from_dataset(wildcards.dataset, key="download_script").endswith(".sh") 
            else "data/raw_data/{dataset}.RData"
    output:
        processed_file="data/processed_adata/{dataset}.h5ad"
    params:
        script=lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key="process_script"),
        input_dir=lambda wildcards: f"data/raw_data/{wildcards.dataset}/",
    wildcard_constraints:
        dataset="[a-zA-Z0-9]+"  # Match only alphanumeric names without underscores
    conda:
        "../envs/process_to_adata.yaml"    
    shell:
        """
        python {params.script} {params.input_dir} {output.processed_file} {input.raw_file}
        """



rule process_to_adata_multi:
    input:
        raw_dir="data/raw_data/{dataset}/.done"
    output: processed_file = "data/processed_adata/{dataset}-{subdata}.h5ad"
    params:
        script=lambda wildcards: cfg.get_from_dataset(wildcards.dataset, key="process_script"),
        input_dir=lambda wildcards: f"data/raw_data/{wildcards.dataset}",
        subdata=lambda wildcards: wildcards.subdata
    conda:
        "../envs/process_to_adata.yaml"   
    shell:
        """
        python {params.script} {params.input_dir} {output.processed_file} {params.subdata}
        """
