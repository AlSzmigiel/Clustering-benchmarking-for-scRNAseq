
#  Environment Setup

To get started, first create the environment using `conda` or `mamba`. This will install `snakemake`, `snakedeploy` and a version of `scanpy` that is compatible with the workflow. Note that `anndata` is pinned to an older version to ensure compatibility with the selected scanpy version.

```bash
mamba create -c conda-forge -c bioconda --name snakemake_base snakemake snakedeploy scanpy=1.10.2 anndata=0.10.2
```

# Use the workflow with `snakedeploy`:

The following instructions follow the guidelines from the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/index.html).

## Create and enter your project working directory:

   ```bash
   mkdir -p path/to/project-workdir
   cd path/to/project-workdir
   ```

## Deploy the workflow using `snakedeploy`:

   ```bash
   snakedeploy deploy-workflow https://github.com/Campello-Lab/Clustering-benchmarking-for-scRNAseq . --tag main
```

This will create two folders:
- `workflow`: contains the deployed Snakemake module
- `config`: contains configuration file
   
## Configure Workflow:
Configuration file `config/config.yaml` contains configuration that will enable reproduction of the benchmarking performed in the paper. Edit `config/config.yaml` to adjust settings. You can comment out datasets to run the workflow on a smaller collection.

# Running the workflow

## Option A: Run Entire Workflow (Including Data Download & Preprocessing)

```bash
snakemake --cores all --sdm conda
```

- On first run, datasets will be downloaded and processed.  
- This will generate:  
  - `data/`  
  - `results/benchmark_analysis/processed_normalized_data/`  
  - `resources/`  

⚠️ The first execution may take a while, as all required environments are created.  

After datasets are downloaded and processed, run again:

```bash
snakemake --cores all --sdm conda
```

This will launch the benchmarking analysis with the specified hyperparameters and algorithms.

---

## Option B: Skipping Data Download & Preprocessing (Recommended)

### To reproduce results from the paper without reprocessing all data, download **preprocessed datasets** from [Zenodo](https://zenodo.org/records/16739211).
----


1. Download and extract `processed_normalized_data.zip` from [Zenodo](https://zenodo.org/records/16739211).    
2. Place it into your project working directory under:

```
results/benchmark_analysis/processed_normalized_data
```

This folder contains processed and normalized h5ad files. Each file has been processed to include 1,000 highly variable genes. The datasets include:
- Enhanced silver-standard datasets  
- Gold-standard datasets 

This is the **entry point** for clustering benchmarking. Now, to start the workflow, run:

```bash
snakemake --cores all --sdm conda
```


---

### 📂 Alternative Workflow Starting Point Options

Snakemake looks for certain files/folders to decide the starting point of the workflow. Depending on which ones you provide, some steps can be skipped.

1. **Raw (unprocessed) datasets**  
- Extract `processed_to_adata.zip` → place into `data processed_to_adata/`. 
- Effect: skips downloading datasets and converting them into `.h5ad` format.
- Workflow will then begin at the stage of:
   - Normalizing gold-standard and silver-standard `.h5ad` files
   - Generating enhanced labels for silver-standard datasets

2. **Skip creation of enhanced labels**
- Extract `enhanced_silver_standard.zip` → place contents into `data/enhanced_silver_standard/`.
- Effect: workflow will use enhanced labels directly (skips their creation).

3. **Speed up enhanced label creation (optional)**

- Extract `resources.zip` → place contents into `resources/`.
- This provides:
   - Pretrained scVI models → `resources/models/`
   - Reference embeddings (for training classifiers) → `resources/ref_embeddings/`
   - Trained classifiers (used to obtain enhanced labels) → `resources/trained_classification_models/`
- Effect: skips downloading reference embeddings.
   - **Still performed**:
      - Generating embeddings with pretrained scVI models
      - Training classifiers using reference + generated embeddings



# 📫 Contact

For questions or support, please open an issue or contact: `alszm@imada.sdu.dk`.


