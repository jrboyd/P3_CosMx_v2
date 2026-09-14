#!/usr/bin/env bash
## Submit scripts/11_script_hodgkin_niche_global.R -- whole-tissue Hodgkin-chemokine
## summary, the counterpart to scripts/07_submit_lightchain_global.sh.
##   scripts/11_submit_hodgkin_niche_global.sh [sample] [bin_um]
## 96G matches the step-07 global job; the whole-slide image fetch is the cost.
set -euo pipefail
cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }
[[ -f output/11_hodgkin_niche/niche_bins_all_scales.csv.gz ]] || {
  echo "run scripts/11_submit_hodgkin_niche.sh first" >&2; exit 1; }
sbatch --mem=96G --cpus-per-task=2 --time=3:00:00 --job-name=nglob \
  scripts/sbatch/run_rscript.sh scripts/11_script_hodgkin_niche_global.R "${1:-IB1}" "${2:-100}"
echo; echo "output: output/11_hodgkin_niche_global/"
