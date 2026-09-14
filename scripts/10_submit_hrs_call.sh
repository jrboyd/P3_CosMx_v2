#!/usr/bin/env bash
## ============================================================================
## Submit scripts/10_script_hrs_call.R -- HRS call for IB1, RB1 as negative control.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   scripts/10_submit_hrs_call.sh
##
## Reads output/04_merge/merged_processed.Rds (1.98 GB, both samples post-QC). The old
## HRS script ran this at 200G, but step 08 peaked at 7 GB against a comparable object
## and this densifies only ~50 marker rows, so 96G is generous -- check MaxRSS.
##
## Run scripts/10_submit_hrs_views.sh afterwards. The call is not a diagnosis until the
## images have been looked at (scripts/01_import_cosmx_seurat.R:25-26).
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }
[[ -f output/04_merge/merged_processed.Rds ]] || { echo "missing merged_processed.Rds" >&2; exit 1; }

sbatch --mem=96G --cpus-per-task=2 --time=3:00:00 --job-name=hrscall \
  scripts/sbatch/run_rscript.sh scripts/10_script_hrs_call.R

echo
echo "watch:  squeue -u \$USER | grep hrscall"
echo "output: output/10_hrs/"
