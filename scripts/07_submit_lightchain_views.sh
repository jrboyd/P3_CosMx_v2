#!/usr/bin/env bash
## ============================================================================
## Submit scripts/07_script_lightchain_views.R for every radius.
##
## This is the record of what produced everything under
## output/07_lightchain_views/ -- keep it in step with that directory. Previously
## these invocations lived only in shell history.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   COSMX_LC_DIR=output/07_lightchain scripts/07_submit_lightchain_views.sh
##   # parallel set restricted to neighbourhoods inside FOVs that pass step 02's
##   # spatial QC (only 145 of 758 IB1 FOVs do), written to ..._fovqc/:
##   COSMX_FOV_QC=1 COSMX_LC_DIR=output/07_lightchain scripts/07_submit_lightchain_views.sh
##   FOCUS_XY=7.424,12.953 COSMX_LC_DIR=output/04_merge scripts/07_submit_lightchain_views.sh
##   scripts/07_submit_lightchain_views.sh              # both modes, all radii
##   scripts/07_submit_lightchain_views.sh default      # default mode only
##   scripts/07_submit_lightchain_views.sh focus        # focus mode only
##
## Two modes, because they answer different questions:
##
##   default  per radius, the most kappa-skewed neighbourhood plus a polytypic
##            control on the expected ~2:1 ratio. NOTE the kappa pick is a
##            DIFFERENT object at each radius (1.0-3.0 mm apart, against spans of
##            0.09-0.30 mm), so this compares scales, not the same tissue.
##
##   focus    one physical focus pinned by slide-mm coordinate and followed across
##            radii. FOCUS_XY below is the r20um kappa-restricted (FL-like) focus,
##            nb 7648 at (7.424, 12.953) -- the only one in
##            output/04_merge/lightchain_restricted_foci_persistence.csv.
##            Tracking it showed the call degrade with radius:
##              r20 kappa_frac 0.895 FL-like -> r30 0.741 intermediate
##              -> r40 0.657 polytypic -> r60 0.656 polytypic
##            i.e. only the 20 um radius resolves it as clonal. See
##            output/07_lightchain_views/focus_across_radii.csv
##
## To render against a re-run of step 04 (e.g. built from a pre-QC object), export
## COSMX_LC_DIR to that output directory before submitting. It changes which cells
## carry a call/neighbourhood; it does NOT change the boundary polygons, which come
## from the raw per-FOV CellBoundaries csvs.
## ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }

SAMPLE=IB1
RADII=(20 30 40 60)
N_EACH=5        # examples per group in default mode: 5 kappa + 5 polytypic controls
## Focus mode stays at 1: the focus is a single pinned location by construction, so
## asking for 5 would return the focus plus its 4 nearest neighbours, which are
## different objects. Raise it only if you want the focus's surroundings too.
N_FOCUS=1
## Default focus: the strongest focus that PERSISTS across radii AFTER plasma cells and
## their spillover halo were dropped from the graph:
##   (8.48, 9.65)  kappa 0.967, 46 cells, found at r20 + r30
## Runner-up (1.36, -0.37) kappa 0.946, 52 cells, also 2 radii -- but it sits at the
## slide edge (y < 0), as do (0.39, 1.57) and (0.63, -0.20), so treat those with care.
##
## HISTORY, because this value has moved twice and the reason matters:
##   7.45,13.99  was the previous default, described as the best 3-radius focus. It is
##               NOT recovered at ANY radius once plasma cells are excluded (nearest
##               scored neighbourhood 0.27-0.38 mm away). It was plasma-cell spillover.
##   6.98,0.36   earlier default, superseded.
##   7.424,12.953 for the QC-filtered output/04_merge tables, which have no plasma tiering.
FOCUS_XY="${FOCUS_XY:-8.48,9.65}"

MODE="${1:-both}"
SB=(sbatch --mem=96G --cpus-per-task=2 --time=4:00:00)   # 3 zooms incl. 10x
RUNNER=scripts/sbatch/run_rscript.sh
SCRIPT=scripts/07_script_lightchain_views.R

for R in "${RADII[@]}"; do
  if [[ "$MODE" == "both" || "$MODE" == "default" ]]; then
    "${SB[@]}" --job-name="lcv-$R" "$RUNNER" "$SCRIPT" "$SAMPLE" "$R" "$N_EACH"
  fi
  if [[ "$MODE" == "both" || "$MODE" == "focus" ]]; then
    "${SB[@]}" --job-name="lcf-$R" "$RUNNER" "$SCRIPT" "$SAMPLE" "$R" "$N_FOCUS" "$FOCUS_XY"
  fi
done

echo
echo "submitted mode=$MODE for radii: ${RADII[*]}"
echo "watch:   squeue -u \$USER | grep -E 'lcv-|lcf-'"
echo "logs:    logs/adhoc_<jobid>.{out,err}"
echo "figures: output/07_lightchain_views/r<R>um/"
