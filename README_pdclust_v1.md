This package contains a PDclust/MDS reproduction of Fig. 1b logic using QC-pass per-cell `.cov.gz` files already produced in the MethSCAn prep step.

Files:
- `run_pdclust__syn22257780_v1.sh` — SBATCH wrapper.
- `script_pdclust__syn22257780_v1.R` — R analysis script.

Expected input files:
- `/home/matpg/scDNAme/results_scDNAme_multimodal_v2_nobam_safe/methscan_all/prep/cov_by_cell/*.cov.gz`
- `/home/matpg/scDNAme/tables/analysis_scRRBS_sequencing_qc.tsv`
- `/home/matpg/scDNAme/tables/clinical_metadata.tsv`

Outputs:
- `pdclust_pairwise.rds`
- `pdclust_dissimilarity_matrix.rds`
- `pdclust_cluster_results.rds`
- `pdclust_cluster_assignments.tsv`
- `PDclust_heatmap.pdf`
- `PDclust_MDS_by_case.pdf`
- `PDclust_MDS_by_IDH.pdf`
- `PDclust_MDS_by_cluster.pdf`
- `pdclust_cluster_by_case.tsv`
- `pdclust_cluster_by_case_fraction.tsv`
