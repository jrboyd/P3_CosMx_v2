#!/usr/bin/env bash
## ============================================================================
## Submit scripts/11_script_hodgkin_niche.R.
##
## Locates the Hodgkin niche in IB1 from the two canonical secreted HRS chemokines
## (CCL17, CCL22), with RB1 -- the FL recurrence, no Hodgkin component -- as the null.
## Regional, not per-cell: step 10 established this panel cannot call HRS cells.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   scripts/11_submit_hodgkin_niche.sh
##
## Reads output/04_merge/merged_processed.Rds (1.98 GB). Step 10 ran the same object
## in under 4 min; 96G is generous -- check MaxRSS and trim.
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }
[[ -f output/04_merge/merged_processed.Rds ]] || { echo "missing merged_processed.Rds" >&2; exit 1; }

sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 --job-name=hniche \
  scripts/sbatch/run_rscript.sh scripts/11_script_hodgkin_niche.R

echo
echo "watch:  squeue -u \$USER | grep hniche"
echo "output: output/11_hodgkin_niche/"
