#!/usr/bin/env bash
## ============================================================================
## 18_script_centroids.sh -- per-cell centroids in CANVAS PX, from the flat metadata file.
##
## WHY NOT output/16_regions/cells_meta_all.csv.gz. That table carries x_slide_mm /
## y_slide_mm, and getting from mm to canvas px needs the FITTED affine
## (cosmx_fit_mm_affine), which is the project's most-repeated source of mirrored
## figures. The vendor's flat metadata file carries CenterX_global_px /
## CenterY_global_px directly -- the SAME space the transcript file's x_global_px /
## y_global_px are in. So step 18 does all of its geometry in canvas px, joins
## transcripts to centroids with no transform at all, and the affine never enters.
##
## Also pulled: the segmentation-shape columns (Area, AspectRatio, Circularity,
## Solidity, Width, Height) and unassignedTranscripts. A cell whose mask is
## mis-drawn is the cell most likely to steal a neighbour's transcripts, so these
## are covariates for the leakage test, not decoration.
##
## Streaming, constant memory, one CPU -- safe on the login node.
##
## Usage: bash scripts/18_script_centroids.sh IB1
## ============================================================================
set -euo pipefail
S="${1:?usage: $0 <IB1|RB1>}"
SRC="/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/${S}/${S}_metadata_file.csv.gz"
OUT="output/18_leakage/${S}_centroids.csv.gz"
[[ -f "$SRC" ]] || { echo "missing $SRC" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

zcat "$SRC" | awk -F, -v OFS=, '
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i
          for (k in h) present[k]=1
          split("cell_id fov CenterX_global_px CenterY_global_px Area Area.um2 " \
                "AspectRatio Circularity Solidity Width Height NucArea " \
                "nCount_RNA nFeature_RNA unassignedTranscripts qcCellsPassed", W, " ")
          for (j=1; j in W; j++) if (!(W[j] in h)) { print "missing column: " W[j] > "/dev/stderr"; exit 1 }
          print "cell_id,fov,cx_px,cy_px,area_px,area_um2,aspect,circularity,solidity," \
                "width_px,height_px,nucarea_px,nCount_RNA,nFeature_RNA,unassigned_tx,qc_passed"
          next }
  { print $h["cell_id"], $h["fov"], $h["CenterX_global_px"], $h["CenterY_global_px"],
          $h["Area"], $h["Area.um2"], $h["AspectRatio"], $h["Circularity"], $h["Solidity"],
          $h["Width"], $h["Height"], $h["NucArea"], $h["nCount_RNA"], $h["nFeature_RNA"],
          $h["unassignedTranscripts"], $h["qcCellsPassed"]
    n++ }
  END { printf("%s: wrote %d cells\n", "'"$S"'", n) > "/dev/stderr" }' \
  | gzip -c > "$OUT"
echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
