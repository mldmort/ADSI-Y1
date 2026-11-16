#!/bin/bash
# Fit the Figure 1 (SPARK, hybrid) and Figure 2 (SSC probands) models as SLURM jobs.
#   bash run_models.sh            # cont_eur (primary) + all_ancestry (sensitivity) + SSC
# Prepared input tables come from the GxE_analysis pipeline (01_prepare_data.py, run from that repo with
# `--subdir cont_eur` and `--ancestry all --subdir all_ancestry`); a copy of that script is kept here.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GXE=${GXE:-/expanse/projects/sebat1/miladm/UCSD/projects/GxE_analysis}
source "$GXE/config/slurm.env"
mkdir -p "$HERE/logs"
sub() {  # name command
  sbatch --parsable --job-name="pwf_$1" --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=48G --time=04:00:00 --exclude=exp-15-32 \
    --output="$HERE/logs/$1_%j.out" --error="$HERE/logs/$1_%j.err" --export=ALL --wrap="cd '$HERE'; export OPENBLAS_NUM_THREADS=1; echo Start \$(date); $2; echo End \$(date)"
}
echo "fig1_eur=$(sub fig1_eur "'$GXE_RSCRIPT' 02_fig1_models.R '$GXE' cont_eur results")"
echo "fig1_all=$(sub fig1_all "'$GXE_RSCRIPT' 02_fig1_models.R '$GXE' all_ancestry results")"
echo "fig2_ssc=$(sub fig2_ssc "'$GXE_RSCRIPT' 03_fig2_ssc.R '$GXE' cont_eur results")"
