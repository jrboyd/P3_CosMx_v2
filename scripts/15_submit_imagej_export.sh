#!/usr/bin/env bash
## ============================================================================
## Submit the two-stage ImageJ/Fiji export for a sample.
##
##   scripts/15_submit_imagej_export.sh                 # IB1, overview L5, 4256 px crops
##   scripts/15_submit_imagej_export.sh IB1 4 4256      # finer overview
##   scripts/15_submit_imagej_export.sh IB1 5 4256 hodgkin_hotspot,fl_nb,fl_focus
##
## Stage 1 (R)      regions -> canvas px          scripts/15_script_imagej_export.R
## Stage 2 (Python) canvas px -> TIFFs + ROIs     scripts/15_imagej_tiles.py
## Verify  (R)      ROI round-trip                scripts/15_verify_imagej_export.R
##
## Stage 2 depends on stage 1, so it is submitted with --dependency=afterok and
## the verify job likewise -- the whole chain can be fired and left alone.
##
## The OME-TIFFs live on SCRATCH, which gets cleaned. Guarded below, same as
## scripts/11_submit_hodgkin_niche_views.sh.
##
## NOTE ON THE OVERVIEW LEVEL. The ImageJ ROI format stores its bounding box as
## int16 (max 32767), so an overview at level 2 or finer cannot carry binary
## ROIs -- 15_imagej_tiles.py says so and falls back to the .ijm macro. Levels
## 3-9 are fine. See scripts/lib/imagej_roi.py.
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }

SAMPLE="${1:-IB1}"
LEVEL="${2:-5}"
CROP_PX="${3:-4256}"
GROUPS="${4:-hodgkin_hotspot,fl_nb}"

TIFF="/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/${SAMPLE}.ome.tiff"
[[ -f "$TIFF" ]] || { echo "$TIFF not found on scratch" >&2; exit 1; }
[[ -f "output/11_hodgkin_niche/niche_bins_called.csv" ]] || {
  echo "missing output/11_hodgkin_niche/ -- run scripts/11_submit_hodgkin_niche.sh first" >&2
  exit 1; }
[[ -f "output/07_lightchain/${SAMPLE}/lightchain_neighborhoods_r40um.csv" ]] || {
  echo "missing output/07_lightchain/${SAMPLE}/ -- run step 07 first" >&2; exit 1; }
[[ -f "output/10_hrs/hrs_segments.csv.gz" ]] || {
  echo "missing output/10_hrs/hrs_segments.csv.gz (needed for the mm->px fit)" >&2
  exit 1; }

mkdir -p logs

J1=$(sbatch --parsable --mem=64G --cpus-per-task=2 --time=1:00:00 --job-name=ijexp \
  --output=logs/15_ijexport_%j.out --error=logs/15_ijexport_%j.err \
  scripts/sbatch/run_rscript.sh scripts/15_script_imagej_export.R "$SAMPLE")
echo "stage 1 (regions -> canvas px):  $J1"

J2=$(sbatch --parsable --dependency=afterok:"$J1" \
  --mem=48G --cpus-per-task=2 --time=4:00:00 --job-name=ijtiles \
  --output=logs/15_ijtiles_%j.out --error=logs/15_ijtiles_%j.err \
  scripts/sbatch/run_python.sh scripts/15_imagej_tiles.py \
  "$SAMPLE" "$LEVEL" "$CROP_PX" "$GROUPS")
echo "stage 2 (TIFFs + ROIs):          $J2"

J3=$(sbatch --parsable --dependency=afterok:"$J2" \
  --mem=16G --cpus-per-task=1 --time=1:00:00 --job-name=ijverify \
  --output=logs/15_ijverify_%j.out --error=logs/15_ijverify_%j.err \
  scripts/sbatch/run_rscript.sh scripts/15_verify_imagej_export.R "$SAMPLE" "$LEVEL")
echo "verify  (ROI round-trip):        $J3"

echo
echo "watch:  squeue -u \$USER | grep -E 'ijexp|ijtiles|ijverify'"
echo "output: output/15_imagej/"
