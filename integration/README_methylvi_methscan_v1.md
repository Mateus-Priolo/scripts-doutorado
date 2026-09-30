# MethylVI package for MethSCAn outputs (syn22257780)

This package adapts `METHYLVI` to the counts already exported by MethSCAn:
- `methylated_sites.csv.gz`
- `total_sites.csv.gz`

The official methylVI workflow models **methylated counts (`mc`)** and **coverage counts (`cov`)** per region and supports nuisance-factor correction via `batch_key` or categorical covariates. The official batch integration tutorial demonstrates this with `MuData`, `mc_layer="mc"`, `cov_layer="cov"`, and a sequencing-protocol covariate. citeturn182658search0turn182658search1

## Key adaptation used here
The official tutorial uses two modalities (`mCG`, `mCH`) aggregated with ALLCools. In your current dataset, RRBS is effectively CpG-only, so this package builds a single-modality `MuData` object with one modality named `mCG`, storing:
- `layers['mc']` = methylated site counts
- `layers['cov']` = total covered sites
- `.X` = methylation fractions
- `obs['case_barcode']` as the batch key

This follows the official `setup_mudata()` interface for MethylVI while adapting it to a single methylation context. citeturn182658search0turn182658search1

## Inputs expected
- `/home/matpg/scDNAme/results_scDNAme_multimodal_v2_nobam_safe/methscan_<MODE>/<MATRIX_DIR>/methylated_sites.csv.gz`
- `/home/matpg/scDNAme/results_scDNAme_multimodal_v2_nobam_safe/methscan_<MODE>/<MATRIX_DIR>/total_sites.csv.gz`
- `/home/matpg/scDNAme/results_scDNAme_multimodal_v2_nobam_safe/methscan_<MODE>/prep/qc_keep.tsv`
- `/home/matpg/scDNAme/tables/clinical_metadata.tsv`

## Outputs
Each mode writes to:
- `.../methscan_<MODE>/methylvi_casebatch_v1/`

Files include:
- `methylvi_latent.tsv`
- `methylvi_umap.tsv`
- `cell_cluster_assignments.tsv`
- `cluster_sizes.tsv`
- `cluster_by_case.tsv`
- `cluster_by_case_fraction.tsv`
- `cell_groups/*.csv`
- `UMAP_clusters.pdf`
- `UMAP_cases.pdf`
- `UMAP_IDH.pdf`
- `methylvi_input_and_latent.h5mu`
- `model/`
- optional `dmrs/*.bed`

## Environment
The run script assumes an environment named `methylvi` by default:
```bash
SCVI_ENV=methylvi bash check_methylvi_env_v1.sh
```

Likely required packages:
- python >= 3.10 recommended
- scvi-tools
- scanpy
- anndata
- mudata
- torch
- scipy
- pandas
- numpy
- matplotlib

## Submit examples
Run on `all`, VMR counts:
```bash
sbatch run_methylvi_methscan_syn22257780_v1.sh all VMR_matrix
```

Run all modes on promoter counts:
```bash
sbatch run_methylvi_methscan_syn22257780_v1.sh all_modes promoter_matrix
```

Tune model size and feature count:
```bash
N_LATENT=12 MAX_EPOCHS=300 MAX_FEATURES=3000 sbatch run_methylvi_methscan_syn22257780_v1.sh all VMR_matrix
```

## Practical recommendation
Start with:
- `MODE=all`
- `MATRIX_DIR=VMR_matrix`
- `MAX_FEATURES=3000` or `5000`
- `N_LATENT=10`

If CPU-only training is too slow, switch to `promoter_matrix` first. The official documentation supports `accelerator`/`devices` arguments during training, so GPU-backed training can be enabled through the environment if available. citeturn182658search1
