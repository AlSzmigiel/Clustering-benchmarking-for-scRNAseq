rule get_scvi_model:
    output:
        "resources/models/{version}-scvi-{organism}/scvi.model/model.pt"
    params:
        version="{version}",  # Model version
        organism="{organism}"  # Organism
    conda:
        "../envs/aws_cli.yaml"  # or whatever path you saved
    shell:
        """
        aws s3 cp --no-sign-request --no-progress --only-show-errors \
        s3://cellxgene-contrib-public/models/scvi/{params.version}/{params.organism}/model.pt \
        {output}
        """

rule get_reference_embeddings:
    output: 
        ref_embedding = "resources/ref_embeddings/{organism}/{tissue}/scvi_embedding.h5ad"
    params:
        version = cfg.get_from_tl('census_version'),
        dataset_ids = cfg.get_data_ids_cellxgene()
    conda: 
        "../envs/cellxgene.yaml"
    script:
        "../scripts/get_ref_embeddings.py"


rule project_to_embeddings:
    input:
        silver_standard = "data/processed_adata/{dataset}.h5ad",
        model = lambda wildcards: f"resources/models/{cfg.get_from_tl('census_version')}-scvi-{cfg.get_from_dataset(wildcards.dataset, key='organism')}/scvi.model"
    output:
        adata_embedding = "data/projected_embeddings/{dataset}-emb.h5ad",
    wildcard_constraints:
        dataset = "|".join(cfg.get_silver_standard())  # Only allow "silver" datasets
    conda: 
        "../envs/scvi_projection.yaml"
    script:
        "../scripts/get_latent_rep.py"

rule get_labels_from_classifier:
    input:
        orig_adata = "data/processed_adata/{dataset}.h5ad",
        adata_embedding = "data/projected_embeddings/{dataset}-emb.h5ad",
        ref_embedding = lambda wildcards: f"resources/ref_embeddings/{cfg.get_from_dataset(wildcards.dataset, key='organism')}/{cfg.get_from_dataset(wildcards.dataset, key='tissue')}/scvi_embedding.h5ad"
    output:
        new_labels = "data/enhanced_silver_standard/{dataset}-new_lab.h5ad"
    conda:
        "../envs/get_enhanced_labels.yaml"
    script:
        "../scripts/get_enhanced_labels.py"
