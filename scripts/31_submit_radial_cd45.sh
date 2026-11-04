#!/usr/bin/env bash
## Step 31: radial CD45 profiles around candidate HRS cells. Three chained jobs:
##   select (R) -> profile (Python, reads the OME-TIFF on scratch) -> analyze (R)
##
##   scripts/31_submit_radial_cd45.sh            # full run, 500 + 500 cells
##   scripts/31_submit_radial_cd45.sh smoke      # 12 + 12 cells, separate output dir
set -euo pipefail
cd "$(dirname "$0")/.."
SMOKE="${1:-}"
OUT="output/31_radial_cd45${SMOKE:+_smoke}"
TIFF=/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff
[[ -f "$TIFF" ]] || { echo "$TIFF not found on scratch" >&2; exit 1; }
[[ -f output/29_bright_cells/source_cells.csv.gz ]] || { echo "run step 29 first" >&2; exit 1; }
mkdir -p logs
J1=$(sbatch --parsable --mem=16G --time=0:30:00 --job-name=s31sel \
  --output=logs/31_select_%j.out --error=logs/31_select_%j.err \
  scripts/sbatch/run_rscript.sh scripts/31_script_radial_cd45.R select $SMOKE)
J2=$(sbatch --parsable --dependency=afterok:$J1 --mem=8G --time=3:00:00 --job-name=s31prof \
  --output=logs/31_profile_%j.out --error=logs/31_profile_%j.err \
  scripts/sbatch/run_python.sh scripts/31_radial_cd45_profiles.py "$OUT")
J3=$(sbatch --parsable --dependency=afterok:$J2 --mem=16G --time=0:30:00 --job-name=s31ana \
  --output=logs/31_analyze_%j.out --error=logs/31_analyze_%j.err \
  scripts/sbatch/run_rscript.sh scripts/31_script_radial_cd45.R analyze $SMOKE)
echo "select $J1 -> profile $J2 -> analyze $J3   (output: $OUT)"
