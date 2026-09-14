#!/usr/bin/env bash
## ============================================================================
## 18_script_tx_extract.sh -- ALL transcripts of the step-18 focal and donor cells.
##
## WHY ALL GENES, not just the marker genes. The leakage statistic is a WITHIN-CELL
## contrast: of the transcripts in this cell, are the argmax-marker ones displaced toward
## the donor relative to the cell's OTHER transcripts? The cell's other transcripts are the
## null distribution, so they have to be extracted too. A marker-gene-only extract (which
## is what 17_script_tx_extract.sh does) would leave the statistic with nothing to permute
## against, and would force a between-cell null -- the exact construction that failed three
## times in step 17.
##
## Filtering by CELL instead of by GENE is what keeps this affordable. ~150k cells at ~200
## transcripts each is ~30M rows out of the ~700M in the file.
##
## Unassigned transcripts (cell_ID 0, CellComp "") cannot belong to a focal cell and are
## dropped by the cell-id filter automatically. They are 4.0% of the file.
##
## Usage: sbatch --mem=8G --time=8:00:00 scripts/18_script_tx_extract.sh IB1
## ============================================================================
#SBATCH --job-name=18tx
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=logs/18tx_%j.out
#SBATCH --error=logs/18tx_%j.err
set -euo pipefail
S="${1:?usage: $0 <IB1|RB1>}"
SRC="/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/${S}/${S}_tx_file.csv.gz"
IDS="output/18_leakage/${S}_extract_ids.txt"
OUT="output/18_leakage/${S}_tx_raw.csv.gz"
[[ -f "$SRC" ]] || { echo "missing $SRC" >&2; exit 1; }
[[ -f "$IDS" ]] || { echo "missing $IDS -- run 18_script_select_cells.R first" >&2; exit 1; }
echo "sample $S | $(wc -l < "$IDS") cell ids | source $(du -h "$SRC" | cut -f1)"

zcat "$SRC" | awk -F, -v IDS="$IDS" '
  BEGIN { while ((getline line < IDS) > 0) K[line]=1
          printf("cell ids loaded: %d\n", length(K)) > "/dev/stderr" }
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i
          split("cell x_global_px y_global_px z target CellComp", W, " ")
          for (j=1; j in W; j++) if (!(W[j] in h)) { print "missing column: " W[j] > "/dev/stderr"; exit 1 }
          print "cell_id,x_px,y_px,z,target,CellComp"; next }
  ($h["cell"] in K) {
    print $h["cell"] "," $h["x_global_px"] "," $h["y_global_px"] "," $h["z"] "," \
          $h["target"] "," $h["CellComp"]
    kept++; seen[$h["cell"]]=1 }
  { tot++ }
  END { printf("scanned %d transcripts, kept %d (%.2f%%) across %d of %d cells\n",
               tot, kept, 100*kept/tot, length(seen), length(K)) > "/dev/stderr" }' \
  | gzip -c > "$OUT"
echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
