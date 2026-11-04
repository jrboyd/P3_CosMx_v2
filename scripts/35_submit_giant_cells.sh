#!/usr/bin/env bash
## Step 35: giant-cell reconstruction. select (R, waits) -> split polygons per FOV (awk)
## -> adjacency (python SLURM array, one FOV per task) -> analyze (R).
##   scripts/35_submit_giant_cells.sh [smoke]     # smoke: 6 FOVs per sample
set -euo pipefail
cd "$(dirname "$0")/.."
SMOKE="${1:-}"
OUT="output/35_giant_cells${SMOKE:+_smoke}"
for s in IB1 RB1; do
  [[ -f /gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/$s.ome.tiff ]] || { echo "$s.ome.tiff not on scratch" >&2; exit 1; }
done
[[ -f output/29_bright_cells/source_cells.csv.gz ]] || { echo "run step 29 first" >&2; exit 1; }
mkdir -p logs "$OUT/polygons" "$OUT/adj"

echo "select (waits)..."
sbatch --wait --mem=48G --time=1:00:00 --job-name=s35sel \
  --output=logs/35_select_%j.out --error=logs/35_select_%j.err \
  scripts/sbatch/run_rscript.sh scripts/35_script_giant_cells.R select $SMOKE
N=$(($(wc -l < "$OUT/fov_list.csv") - 1))
echo "FOVs to process: $N"

echo "splitting polygons per FOV..."
for s in IB1 RB1; do
  zcat input/CosMx/export1/flatFiles/$s/$s-polygons.csv.gz | \
    awk -F, -v d="$OUT/polygons" -v s="$s" 'NR==FNR{f[$1]; next} FNR>1 && ($1 in f){print $1","$2","$4","$5 > (d"/"s"_"$1".csv")}' \
    "$OUT/fovs_$s.txt" -
done
ls "$OUT/polygons" | wc -l

J2=$(sbatch --parsable --array=0-$((N-1))%16 --mem=6G --time=0:30:00 --job-name=s35adj \
  --output=logs/35_adj_%A_%a.out --error=logs/35_adj_%A_%a.err \
  scripts/sbatch/run_python.sh scripts/35_fov_adjacency.py "$OUT")
J3=$(sbatch --parsable --dependency=afterok:$J2 --mem=48G --time=4:00:00 --job-name=s35ana \
  --output=logs/35_analyze_%j.out --error=logs/35_analyze_%j.err \
  scripts/sbatch/run_rscript.sh scripts/35_script_giant_cells.R analyze $SMOKE)
echo "adjacency array $J2 -> analyze $J3   (output: $OUT)"
