#!/usr/bin/env bash
#SBATCH --job-name=cosmx-03
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/03_%j.out
#SBATCH --error=logs/03_%j.err

set -euo pipefail

# Example sbatch submission for running one sample through
# scripts/03_script_filter_and_process.R in parallel.
#
# Usage:
#   sbatch scripts/03_submit.sh /path/to/input_seurat.Rds sample_id
#
# Example:
#   sbatch scripts/03_submit.sh output/01_import_v1/seuratObject_IB1_BPCells_full.Rds IB1

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Got $# arguments: $*" >&2
  echo "Usage: sbatch $0 <input_seurat_rds> <sample_id> (<output_prefix>)" >&2
  
  exit 1
fi

INPUT_RDS="$1"
SAMPLE_ID="$2"
OUTPUT_PREFIX="${3:-$SAMPLE_ID}"

R_SCRIPT="/gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/R/4.5.1/bin/Rscript"
R_LIBS_USER="/gpfs1/home/j/r/jrboyd/R/x86_64-pc-linux-gnu-library/4.5"
export R_LIBS_USER

mkdir -p logs output/03_filter_and_process/figures

"$R_SCRIPT" --vanilla scripts/03_script_filter_and_process.R "$INPUT_RDS" "$SAMPLE_ID" "$OUTPUT_PREFIX"
