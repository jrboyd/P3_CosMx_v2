#!/usr/bin/env bash
## Step 32: CD45 to 500 um around step 31's cells, defect burden as covariate.
##   select (R, defect counts per band) -> profile (Python, level-3 image) -> analyze (R)
##   scripts/32_submit_wide_cd45.sh [smoke]      # smoke uses step 31's smoke cells
set -euo pipefail
cd "$(dirname "$0")/.."
SMOKE="${1:-}"
OUT="output/32_wide_cd45${SMOKE:+_smoke}"
S31="output/31_radial_cd45${SMOKE:+_smoke}"
[[ -f /gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff ]] || { echo "IB1.ome.tiff not on scratch" >&2; exit 1; }
[[ -f "$S31/cells.csv" ]] || { echo "run scripts/31_submit_radial_cd45.sh $SMOKE first" >&2; exit 1; }
mkdir -p logs "$OUT"
J1=$(sbatch --parsable --mem=24G --time=0:30:00 --job-name=s32sel \
  --output=logs/32_select_%j.out --error=logs/32_select_%j.err \
  scripts/sbatch/run_rscript.sh scripts/32_script_wide_cd45.R select $SMOKE)
J2=$(sbatch --parsable --dependency=afterok:$J1 --mem=8G --time=1:00:00 --job-name=s32prof \
  --output=logs/32_profile_%j.out --error=logs/32_profile_%j.err \
  scripts/sbatch/run_python.sh scripts/32_wide_cd45_profiles.py "$OUT")
J3=$(sbatch --parsable --dependency=afterok:$J2 --mem=16G --time=0:30:00 --job-name=s32ana \
  --output=logs/32_analyze_%j.out --error=logs/32_analyze_%j.err \
  scripts/sbatch/run_rscript.sh scripts/32_script_wide_cd45.R analyze $SMOKE)
echo "select $J1 -> profile $J2 -> analyze $J3   (output: $OUT)"
