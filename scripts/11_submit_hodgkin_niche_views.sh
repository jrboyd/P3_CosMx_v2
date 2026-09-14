#!/usr/bin/env bash
## ============================================================================
## Submit scripts/11_script_hodgkin_niche_views.R -- image views of the CCL17/CCL22
## hotspots, with IB1-background and RB1-top controls rendered identically.
##
##   scripts/11_submit_hodgkin_niche_views.sh            # 12 hotspots + 4 of each control
##   scripts/11_submit_hodgkin_niche_views.sh 20 6
##
## Requires scripts/11_submit_hodgkin_niche.sh to have completed.
## 96G matches the established views jobs; the OME-TIFF fetch is why this cannot run
## on the login node. TIFFs are on SCRATCH.
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }
[[ -f output/11_hodgkin_niche/niche_bins_called.csv ]] || {
  echo "missing output/11_hodgkin_niche/ -- run scripts/11_submit_hodgkin_niche.sh first" >&2
  exit 1; }
for T in IB1 RB1; do
  [[ -f "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/$T.ome.tiff" ]] || {
    echo "$T.ome.tiff not found on scratch" >&2; exit 1; }
done

sbatch --mem=96G --cpus-per-task=2 --time=4:00:00 --job-name=nview \
  scripts/sbatch/run_rscript.sh scripts/11_script_hodgkin_niche_views.R "${1:-12}" "${2:-4}"

echo
echo "watch:  squeue -u \$USER | grep nview"
echo "output: output/11_hodgkin_niche_views/"
