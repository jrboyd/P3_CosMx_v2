#!/usr/bin/env bash
#SBATCH --job-name=r-adhoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=logs/adhoc_%j.out
#SBATCH --error=logs/adhoc_%j.err

# Generic sbatch wrapper for ad-hoc R work that is too heavy for the login node
# (anything over 1 CPU or 4 GB RAM -- which includes reading any of this
# project's Seurat .Rds files).
#
# Usage, from the project root:
#   sbatch scripts/sbatch/run_rscript.sh path/to/script.R [args...]
#   sbatch scripts/sbatch/run_rscript.sh -e 'library(Seurat); ...'
#
# Override the defaults above on the command line as needed, e.g.
#   sbatch --mem=150G --cpus-per-task=4 --time=6:00:00 scripts/sbatch/run_rscript.sh foo.R
#
# NOTE: compute nodes cannot see the login node's node-local /tmp, so the script
# passed here (and anything it reads) must live on shared GPFS.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch $0 <script.R> [args...] | sbatch $0 -e '<R code>'" >&2
  exit 1
fi

R_SCRIPT="/gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/R/4.5.1/bin/Rscript"
R_LIBS_USER="/gpfs1/home/j/r/jrboyd/R/x86_64-pc-linux-gnu-library/4.5"
export R_LIBS_USER

mkdir -p logs

if [[ "$1" == "-e" ]]; then
  shift
  "$R_SCRIPT" --vanilla -e "$@"
else
  "$R_SCRIPT" --vanilla "$@"
fi
