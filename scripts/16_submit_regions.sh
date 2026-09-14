#!/usr/bin/env bash
## ============================================================================
## Submit the step-16 region builder for both QC tiers.
##
## Builds the canonical per-cell region table for IB1 and RB1: disease regions
## (cHL_hotspot / FL_focus / buffer / background in IB1) and the shared architecture
## axis (B_dense / B_mid / B_sparse, within-sample f_B quantiles) in both.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   scripts/16_submit_regions.sh [strict|matched|preqc|both|all]
##
## Prerequisites (each ~2 min, ~7 GB):
##   sbatch --mem=140G scripts/sbatch/run_rscript.sh scripts/16_script_cell_meta.R
##   sbatch --mem=160G scripts/sbatch/run_rscript.sh scripts/16_script_panel_counts.R
##
## Pure CSV + one sparse matrix, so this is far lighter than the object reads above.
## The matched tier carries 1.88M IB1 cells, hence the headroom.
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }
for f in output/16_regions/cells_meta_all.csv.gz output/16_regions/panel_counts_all.Rds; do
  [[ -f "$f" ]] || { echo "missing $f -- run the prerequisites above" >&2; exit 1; }
done

WHICH="${1:-both}"
case "$WHICH" in
  strict|matched|preqc) TIERS=("$WHICH") ;;
  both)           TIERS=(strict matched) ;;
  all)            TIERS=(strict matched preqc) ;;
  *) echo "usage: $0 [strict|matched|preqc|both|all]" >&2; exit 1 ;;
esac

for T in "${TIERS[@]}"; do
  sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 --job-name="16reg-$T" \
    --output="logs/16reg_${T}_%j.out" --error="logs/16reg_${T}_%j.err" \
    scripts/sbatch/run_rscript.sh scripts/16_script_regions.R "$T"
done

echo
echo "watch:  squeue -u \$USER | grep 16reg"
echo "output: output/16_regions/{strict,matched}/"
