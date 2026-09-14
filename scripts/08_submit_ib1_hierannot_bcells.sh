#!/usr/bin/env bash
## ============================================================================
## Submit scripts/08_script_ib1_hierannot_bcells.R.
##
## Redoes the clustered-vs-diffuse comparison on a POSITIVE B-cell definition
## (step 09's HierAnnot call) instead of step 07's exclusion-defined `B_clean` tier,
## which failed the B-identity negative control in 08a.
##
## Requires, for the sample:
##   scripts/07_submit_lightchain_neighborhoods.sh  -> output/07_lightchain/<S>/
##   scripts/09_submit_hierannot_bcell_call.sh      -> output/09_hierannot/<S>/
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   scripts/08_submit_ib1_hierannot_bcells.sh
##
## Memory: 08a peaked at 7 GB against the same 1.9 GB object, so 64G is already
## generous. Supersedes 08_script_ib1_clustered_vs_diffuse.R, which is kept because
## its control failure is itself the result that motivated step 09.
## ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }

SAMPLES=("$@")
[[ ${#SAMPLES[@]} -gt 0 ]] || SAMPLES=(IB1)

for S in "${SAMPLES[@]}"; do
  for need in "output/07_lightchain/$S" "output/09_hierannot/$S/bcell_call.csv.gz"; do
    [[ -e "$need" ]] || { echo "missing $need" >&2; exit 1; }
  done
  sbatch --mem=64G --cpus-per-task=2 --time=3:00:00 --job-name="hab-$S" \
    scripts/sbatch/run_rscript.sh scripts/08_script_ib1_hierannot_bcells.R "$S"
done

echo
echo "submitted: ${SAMPLES[*]}"
echo "watch:  squeue -u \$USER | grep hab-"
echo "output: output/08_ib1_hierannot_bcells/"
