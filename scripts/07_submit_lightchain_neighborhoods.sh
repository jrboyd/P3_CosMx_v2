#!/usr/bin/env bash
## ============================================================================
## Submit scripts/07_script_lightchain_neighborhoods.R.
##
## Stage 1 of step 7: calls the light-chain neighbourhoods from step 1's PRE-QC-
## FILTER object and writes output/07_lightchain/. Run this before
## scripts/07_submit_lightchain_views.sh, which renders image views of the result.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   scripts/07_submit_lightchain_neighborhoods.sh           # IB1
##   scripts/07_submit_lightchain_neighborhoods.sh IB1 RB1   # both samples
##
## Then point the views at this output instead of step 04's:
##   COSMX_LC_DIR=output/07_lightchain scripts/07_submit_lightchain_views.sh
##
## Memory: the pre-filter IB1 object is 2.54M cells (1.9 GB on disk, and the
## BPCells counts are read for two genes only). 200G is comfortable; the earlier
## step-04 version ran on 128G against 400k cells.
## ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }

SAMPLES=("$@")
[[ ${#SAMPLES[@]} -gt 0 ]] || SAMPLES=(IB1)

for S in "${SAMPLES[@]}"; do
  sbatch --mem=200G --cpus-per-task=2 --time=4:00:00 --job-name="lcnb-$S" \
    scripts/sbatch/run_rscript.sh \
    scripts/07_script_lightchain_neighborhoods.R "$S"
done

echo
echo "submitted: ${SAMPLES[*]}"
echo "watch:  squeue -u \$USER | grep lcnb-"
echo "logs:   logs/adhoc_<jobid>.{out,err}"
echo "output: output/07_lightchain/"
