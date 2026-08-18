#!/bin/bash
set -euo pipefail

cd /home/renanomete/projetos/matdata/met_mat_data/singlecell_novo/scripts

JOB1=$(sbatch --parsable 01b_load_gse173278_primary.sh)
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} 02b_prepare_gse173278_primary.sh)
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} 03b_harmony_integration_3datasets.sh)
JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} 04b_downstream_analysis_3datasets.sh)

echo "Pipeline V2 submetido:"
echo "${JOB1} -> ${JOB2} -> ${JOB3} -> ${JOB4}"
