#!/usr/bin/env bash
## ============================================================================
## 18_script_nuclear_counts.sh -- compartment-split marker and negprobe counts.
##
## assign_lineage() only ever touches the marker genes (scripts/lib/lineage.R: it does
## M <- cts[unique(unlist(sets)), ] and never looks at the rest) plus a per-cell negprobe
## total. So recalling lineages on a compartment subset does NOT need the 102.7M-row
## transcript table in R -- it needs a ~1.3M-row aggregate. One streaming awk pass builds
## it, and the R step stays small.
##
## Negative1..Negative20 ARE present in the transcript file, so the background floor can be
## recomputed WITHIN each compartment rather than assumed or rescaled. That matters: the
## floor is what stops a dim lineage being minted from spillover, and an arm that keeps 79%
## of transcripts must be judged against 79% of the background, not all of it.
##
## Two outputs:
##   *_nuc_totals.csv.gz    cell_id, n_nuc, n_cyt, n_mem   -- all transcripts, per compartment
##   *_nuc_markers.csv.gz   cell_id, target, comp, n       -- the 34 coarse markers + negprobes
##
## Usage: sbatch --mem=8G --time=2:00:00 scripts/18_script_nuclear_counts.sh IB1
## ============================================================================
#SBATCH --job-name=18nuc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=logs/18nuc_%j.out
#SBATCH --error=logs/18nuc_%j.err
set -euo pipefail
S="${1:?usage: $0 <IB1|RB1>}"
SRC="output/18_leakage/${S}_tx_raw.csv.gz"
TOT="output/18_leakage/${S}_nuc_totals.csv.gz"
MRK="output/18_leakage/${S}_nuc_markers.csv.gz"
[[ -f "$SRC" ]] || { echo "missing $SRC" >&2; exit 1; }

## LINEAGE_SETS_COARSE, verbatim from scripts/lib/lineage.R:95-102. The R step asserts this
## list matches, so the two cannot drift apart silently.
GENES="MS4A1,CD79A,CD79B,BANK1,CD3D,CD3E,CD3G,TRAC,NKG7,GNLY,KLRD1,NCR1,CD68,CD163,LYZ,ITGAX,C1QA,C1QB,CD14,PECAM1,VWF,CLDN5,CD34,COL1A1,COL1A2,PDGFRB,DCN,ACTA2,MZB1,XBP1,PRDM1,SDC1,DERL3,TNFRSF17"

zcat "$SRC" | awk -F, -v GENES="$GENES" -v MRK="/dev/stdout" -v TOT="$TOT.tmp" '
  BEGIN { n=split(GENES,a,","); for(i=1;i<=n;i++) K[a[i]]=1
          print "cell_id,target,comp,n" }
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
  { cell=$h["cell_id"]; tgt=$h["target"]; cc=$h["CellComp"]
    if (cc == "") next                       # unassigned; cannot reach a focal cell anyway
    T[cell SUBSEP cc]++
    if (tgt in K)            M[cell SUBSEP tgt SUBSEP cc]++
    else if (tgt ~ /^Negative[0-9]+$/) { M[cell SUBSEP "NEGPROBE" SUBSEP cc]++; NEG[tgt]=1 }
    rows++ }
  END {
    for (k in M) { split(k, p, SUBSEP); print p[1] "," p[2] "," p[3] "," M[k] }
    print "cell_id,n_nuc,n_cyt,n_mem" > TOT
    for (k in T) { split(k, p, SUBSEP); C[p[1]]=1 }
    for (c in C) printf("%s,%d,%d,%d\n", c, T[c SUBSEP "Nuclear"]+0,
                        T[c SUBSEP "Cytoplasm"]+0, T[c SUBSEP "Membrane"]+0) > TOT
    nn=0; for (g in NEG) nn++
    printf("rows %d | cells %d | distinct negprobes %d\n", rows, length(C), nn) > "/dev/stderr"
  }' | gzip -c > "$MRK"
gzip -c "$TOT.tmp" > "$TOT" && rm -f "$TOT.tmp"
echo "wrote $MRK ($(du -h "$MRK" | cut -f1))  and  $TOT ($(du -h "$TOT" | cut -f1))"
