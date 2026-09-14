#!/usr/bin/env bash
## ============================================================================
## Submit scripts/07_script_lightchain_global.R -- one whole-tissue summary figure
## per radius, written to output/07_lightchain_global/<in_dir>/.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   COSMX_LC_DIR=output/07_lightchain scripts/07_submit_lightchain_global.sh
##   COSMX_LC_DIR=output/07_lightchain scripts/07_submit_lightchain_global.sh IB1 20
##
## The window is fixed in the script (canvas px x 15,000-100,000, y 5,000-90,000),
## so unlike the per-neighbourhood views there is no selection mode and no zoom
## series -- one figure per radius, and the radius only affects panel 2's
## neighbourhood assignment.
##
## Memory is the same order as the per-neighbourhood views: the image is fetched at a
## low resolution level because the window is ~85,000 px wide, and panels 2 and 3 draw
## cell centroids rather than polygons, so nothing scales with window area except the
## point count.
## ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }

SAMPLE="${1:-IB1}"
shift || true
RADII=("$@")
[[ ${#RADII[@]} -gt 0 ]] || RADII=(20 30 40 60)

for R in "${RADII[@]}"; do
  sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 --job-name="lcg-$R" \
    scripts/sbatch/run_rscript.sh \
    scripts/07_script_lightchain_global.R "$SAMPLE" "$R"
done

echo
echo "submitted $SAMPLE for radii: ${RADII[*]}"
echo "watch:   squeue -u \$USER | grep lcg-"
echo "figures: output/07_lightchain_global/"
