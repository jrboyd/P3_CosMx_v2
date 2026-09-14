#!/usr/bin/env bash
#SBATCH --job-name=cosmx-04
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/04_%j.out
#SBATCH --error=logs/04_%j.err

# Step 4: merge the step-03 cell-typed objects and run the joint scPearsonPCA
# workflow (process_sample from step 3) on ~1.8M cells.
#
#   sbatch scripts/sbatch/04_merge_wrapper.sh                 # IB1 + RB1
#   sbatch scripts/sbatch/04_merge_wrapper.sh IB1 RB1         # explicit
#
# Pre-flight on a small subset of the same code path (minutes, not hours):
#   MERGE_TEST_NCELLS=2000 sbatch --mem=96G --time=1:00:00 \
#     --export=ALL,MERGE_TEST_NCELLS=2000 scripts/sbatch/04_merge_wrapper.sh

set -euo pipefail

R_SCRIPT="/gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/R/4.5.1/bin/Rscript"
R_LIBS_USER="/gpfs1/home/j/r/jrboyd/R/x86_64-pc-linux-gnu-library/4.5"
export R_LIBS_USER

# uwot/RcppParallel and OpenBLAS auto-detect the node's core count and ignore the
# cgroup, so UMAP ran 32 threads inside a 4-CPU allocation during the pre-flight.
# Pin every thread pool to what SLURM actually granted.
NCPU="${SLURM_CPUS_PER_TASK:-1}"
export RCPP_PARALLEL_NUM_THREADS="$NCPU"
export OMP_NUM_THREADS="$NCPU"
export OPENBLAS_NUM_THREADS="$NCPU"
export MKL_NUM_THREADS="$NCPU"
echo "pinned thread pools to ${NCPU} cpus"

mkdir -p logs output/04_merge

"$R_SCRIPT" --vanilla scripts/04_script_merge_and_process.R "$@"
