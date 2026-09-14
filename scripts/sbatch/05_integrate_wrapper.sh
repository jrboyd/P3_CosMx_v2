#!/usr/bin/env bash
#SBATCH --job-name=cosmx-05
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=logs/05_%j.out
#SBATCH --error=logs/05_%j.err
set -euo pipefail
R_SCRIPT="/gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/R/4.5.1/bin/Rscript"
export R_LIBS_USER="/gpfs1/home/j/r/jrboyd/R/x86_64-pc-linux-gnu-library/4.5"
NCPU="${SLURM_CPUS_PER_TASK:-1}"
export RCPP_PARALLEL_NUM_THREADS="$NCPU" OMP_NUM_THREADS="$NCPU" \
       OPENBLAS_NUM_THREADS="$NCPU" MKL_NUM_THREADS="$NCPU"
mkdir -p logs output/05_integrate
"$R_SCRIPT" --vanilla scripts/05_script_integrate.R "$@"
