#!/usr/bin/env bash
## ============================================================================
## Submit scripts/10_script_hrs_views.R -- image views of the called HRS objects,
## with RB1 and below-threshold controls rendered the same way.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   scripts/10_submit_hrs_views.sh            # IB1, 12 HRS + 6 of each control
##   scripts/10_submit_hrs_views.sh IB1 20 8
##
## Requires scripts/10_submit_hrs_call.sh to have completed (output/10_hrs/).
##
## Memory: 96G matches the established views job (07_submit_lightchain_views.sh).
## Windows here are far smaller than step 07's, but each cosmx_gene_points() call
## still reads ~120 MB of transcript csv per FOV, and the OME-TIFF fetch is why this
## cannot run on the login node.
##
## The TIFFs live on SCRATCH -- /gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs --
## so confirm IB1.ome.tiff is still there before a long render campaign.
## ============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }
[[ -f output/10_hrs/hrs_objects.csv.gz ]] || {
  echo "missing output/10_hrs/hrs_objects.csv.gz -- run scripts/10_submit_hrs_call.sh first" >&2
  exit 1; }
[[ -f /gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff ]] || {
  echo "IB1.ome.tiff not found on scratch -- the images may have been cleaned up" >&2
  exit 1; }

SAMPLE="${1:-IB1}"; N_HRS="${2:-12}"; N_CTRL="${3:-6}"

sbatch --mem=96G --cpus-per-task=2 --time=4:00:00 --job-name="hrsview" \
  scripts/sbatch/run_rscript.sh scripts/10_script_hrs_views.R "$SAMPLE" "$N_HRS" "$N_CTRL"

echo
echo "watch:  squeue -u \$USER | grep hrsview"
echo "output: output/10_hrs_views/$SAMPLE/"
