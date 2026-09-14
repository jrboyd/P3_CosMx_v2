#!/usr/bin/env bash
#SBATCH --job-name=cosmx-19a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=logs/19a_%j.out
#SBATCH --error=logs/19a_%j.err
#
# Step 19a: build the high-confidence object (symmetric tier, defect dropped) and
# cluster it jointly.  ~2.05M cells -- comparable to the 1.8M step 04 ran at 200G/24h.
#
#   sbatch scripts/sbatch/19_hc_object_wrapper.sh
#
# Two .done markers in output/19_celltypes/ make this resumable: a failure in the PCA
# does not redo the merge.
set -euo pipefail
R_SCRIPT="/gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/R/4.5.1/bin/Rscript"
R_LIBS_USER="/gpfs1/home/j/r/jrboyd/R/x86_64-pc-linux-gnu-library/4.5"
export R_LIBS_USER

# uwot/RcppParallel and OpenBLAS auto-detect the node's core count and ignore the cgroup,
# so UMAP will run 32 threads inside an 8-CPU allocation unless every pool is pinned.
NCPU="${SLURM_CPUS_PER_TASK:-1}"
export RCPP_PARALLEL_NUM_THREADS="$NCPU" OMP_NUM_THREADS="$NCPU"
export OPENBLAS_NUM_THREADS="$NCPU" MKL_NUM_THREADS="$NCPU"
echo "pinned thread pools to ${NCPU} cpus"

mkdir -p logs output/19_celltypes
"$R_SCRIPT" --vanilla scripts/19_script_hc_object.R "$@"
