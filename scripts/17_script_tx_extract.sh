#!/usr/bin/env bash
## ============================================================================
## Stream the per-transcript files and keep only the genes step 17 needs.
##
## The tx files are 13.0 GB (IB1) and 11.7 GB (RB1) GZIPPED -- hundreds of millions of rows,
## far past loading. One awk pass per sample filters to a target list and writes a few hundred
## MB. Streaming, so memory is constant; the cost is I/O.
##
## WHY THIS EXISTS AT ALL: the per-cell PD-1 analysis threw away every transcript that
## segmentation failed to assign (CellComp == "None"), and thresholded the rest at >=2 counts
## to clear the blank-probe floor. Transcripts have neither problem -- no cell boundary, no
## threshold -- so a transcript-level readout is the one measurement of this axis that is not
## bounded by the detection floor. CellComp is kept so the unassigned fraction is visible.
##
## Usage: sbatch --mem=8G --time=8:00:00 scripts/17_script_tx_extract.sh IB1
## ============================================================================
#SBATCH --job-name=17tx
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=logs/17tx_%j.out
#SBATCH --error=logs/17tx_%j.err
set -euo pipefail
S="${1:?usage: $0 <IB1|RB1>}"
SRC="/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/${S}/${S}_tx_file.csv.gz"
OUT="output/17_pd1/tx/${S}_tx_pd1.csv.gz"
[[ -f "$SRC" ]] || { echo "missing $SRC" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

## The PD-1 axis, its ligands, and the module genes whose transcripts are worth having.
## Kept deliberately short: every extra gene multiplies the output size.
KEEP="PDCD1,CD274,PDCD1LG2,CTLA4,LAG3,TIGIT,HAVCR2,TOX,CXCL13,CXCL9,CXCL10,STAT1,IRF1,GBP1,IFNG,JAK2,CD3E,MS4A1,CD68,BATF3,IGKC"
echo "sample $S"; echo "targets: $KEEP"
zcat "$SRC" | awk -v keep="$KEEP" -F, '
  BEGIN { n=split(keep,a,","); for(i=1;i<=n;i++) K[a[i]]=1 }
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i
          print "fov,cell,x_global_px,y_global_px,z,target,CellComp"; next }
  ($h["target"] in K) {
    print $h["fov"] "," $h["cell"] "," $h["x_global_px"] "," $h["y_global_px"] "," \
          $h["z"] "," $h["target"] "," $h["CellComp"]
    c[$h["target"]]++; kept++
  }
  END { printf("kept %d transcripts of %d rows\n", kept, NR-1) > "/dev/stderr"
        for (g in c) printf("  %-10s %d\n", g, c[g]) > "/dev/stderr" }' \
  | gzip -c > "$OUT"
echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
