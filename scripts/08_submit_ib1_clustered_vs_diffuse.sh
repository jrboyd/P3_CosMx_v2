#!/usr/bin/env bash
## ============================================================================
## Submit scripts/08_script_ib1_clustered_vs_diffuse.R.
##
## IB1 only: compares the marker phenotype and the 30 um microenvironment of B cells
## in kappa-restricted foci / other scored neighbourhoods / diffuse tissue. Requires
## scripts/07_submit_lightchain_neighborhoods.sh to have been run for the sample, since
## the groups come from output/07_lightchain/<sample>/lightchain_{cells,neighborhoods}_r*um.csv.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   scripts/08_submit_ib1_clustered_vs_diffuse.sh           # IB1
##   scripts/08_submit_ib1_clustered_vs_diffuse.sh IB1 RB1   # both, if ever wanted
##
## Memory: same 1.9 GB pre-filter object step 07 reads, but this pulls ~55 marker rows
## rather than 13 and adds a kd-tree over all 2.54M cells for the microenvironment.
## 200G matches step 07 and leaves headroom; check MaxRSS afterwards and trim if it is
## nowhere near that.
##
## This is a SEPARATE analysis from scripts/08_script_lightchain_sample_compare.R
## (IB1 vs RB1, login-node cheap, writes output/08_lightchain_compare/). That one is
## untouched -- its result is what motivates this script.
## ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }

SAMPLES=("$@")
[[ ${#SAMPLES[@]} -gt 0 ]] || SAMPLES=(IB1)

for S in "${SAMPLES[@]}"; do
  [[ -d "output/07_lightchain/$S" ]] || {
    echo "missing output/07_lightchain/$S -- run scripts/07_submit_lightchain_neighborhoods.sh $S first" >&2
    exit 1
  }
  sbatch --mem=200G --cpus-per-task=2 --time=4:00:00 --job-name="ib1cd-$S" \
    scripts/sbatch/run_rscript.sh \
    scripts/08_script_ib1_clustered_vs_diffuse.R "$S"
done

echo
echo "submitted: ${SAMPLES[*]}"
echo "watch:  squeue -u \$USER | grep ib1cd-"
echo "logs:   logs/adhoc_<jobid>.{out,err}"
echo "output: output/08_ib1_clustered_vs_diffuse/"
