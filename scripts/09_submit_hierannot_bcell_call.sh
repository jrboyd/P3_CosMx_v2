#!/usr/bin/env bash
## ============================================================================
## Submit scripts/09_script_hierannot_bcell_call.R.
##
## Cluster-level B-cell call for IB1 via HierAnnot (Python, through reticulate
## against the `hierannot_env` virtualenv in ~/.virtualenvs). Reads the POST-QC
## step-03 object and annotates pearson_clusters_low and pearson_clusters_high
## against every builtin hierarchy that matches immune / tonsil / tme.
##
##   cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
##   scripts/09_submit_hierannot_bcell_call.sh           # IB1
##   scripts/09_submit_hierannot_bcell_call.sh IB1 RB1
##
## Exists because step 08 failed its B-identity negative control: step 07's
## `B_clean` tier is defined by exclusion and never tests for B markers.
##
## Memory: IB1_03_HCA_cell_typed.Rds is 1.3 GB on disk (400k cells x 6175 genes);
## AverageExpression over the BPCells counts is the peak. 128G is generous --
## check MaxRSS and trim. RB1's object is 3.9 GB and would need more.
## ============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
[[ -f scripts/lib/cosmx_tiff.R ]] || { echo "run from the project root" >&2; exit 1; }
[[ -d "$HOME/.virtualenvs/hierannot_env" ]] || {
  echo "missing ~/.virtualenvs/hierannot_env -- see /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx/10a_HierAnnot.R" >&2
  exit 1
}

SAMPLES=("$@")
[[ ${#SAMPLES[@]} -gt 0 ]] || SAMPLES=(IB1)

for S in "${SAMPLES[@]}"; do
  sbatch --mem=128G --cpus-per-task=2 --time=3:00:00 --job-name="hann-$S" \
    scripts/sbatch/run_rscript.sh \
    scripts/09_script_hierannot_bcell_call.R "$S"
done

echo
echo "submitted: ${SAMPLES[*]}"
echo "watch:  squeue -u \$USER | grep hann-"
echo "logs:   logs/adhoc_<jobid>.{out,err}"
echo "output: output/09_hierannot/"
