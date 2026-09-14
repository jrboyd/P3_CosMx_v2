## ============================================================================
## 18_script_select_cells.R -- focal cells + donor geometry for the leakage test.
##
## THE QUESTION. argmax types a cell from 4-7 marker genes. If segmentation gives cell c a
## sliver of its neighbour's cytoplasm, c inherits a handful of the NEIGHBOUR's markers and
## argmax may call c that type, while the other ~6000 genes still cluster c correctly. That
## predicts a specific, checkable geometry: in cells where argmax and the cluster DISAGREE,
## the argmax-marker transcripts should sit displaced TOWARD the neighbour that donated them.
## In cells where they agree, they should not.
##
## THE ARBITER IS pearson_clusters_low, NOT inSituType. inSituType puts 65% of non-defect
## IB1 in two clusters, one of them `ambiguous` with margins of 0.009-0.034, and has no NK
## cluster at all -- argmax=NK could never agree with it, mechanically. pearson_clusters_low
## has 9 clusters, 7 `confident` in the step-09b audit, and gives 65-85% agreement on the
## main lineages. Its cluster 0 (low-information/defect) and cluster 6 (pDC, which has no
## argmax counterpart on the coarse sets) are dropped rather than mapped.
##
## NK IS CARRIED BUT NOT SCOREABLE. pearson_clusters_low has no NK cluster, so every argmax
## NK cell counts as a disagreement by construction. They are extracted (pearson_clusters_high
## cluster 7 IS NK, so a secondary pass can score them) and flagged `scoreable = FALSE`.
##
## HOW A DONOR IS DEFINED, AND WHY NOT BY CELL TYPE. The obvious donor is "the nearest cell
## that argmax and the cluster agree is type A". That fails on density: only 183,091 of
## IB1's 2,540,021 cells carry a pearson label at all (the object is the 400,111-cell
## qc_passing_all tier, 45% of which is the defect cluster). Nearest-labelled-neighbour
## distances would run 30-60 um -- three to six cell diameters, not adjacent, and the whole
## test is about adjacency.
##
## So a donor is defined PHYSICALLY: the nearest cell, among all 2.5M non-defect cells,
## carrying >= MIN_DONOR transcripts of family A. That needs no cell-type call on the donor,
## it is the quantity that actually matters (a cell can only donate markers it has), and it
## comes with a free continuous covariate -- how loaded the donor is.
##
## THREE DIRECTIONS ARE STORED, NOT ONE. Leakage is direction-SPECIFIC; a mis-drawn mask or
## an elongated cell is not. So each focal cell gets:
##   theta_A     toward the nearest A-marker-rich cell   <- the test direction
##   theta_ctrl  toward the nearest C-marker-rich cell, C not in {A, cluster}  <- control
##   theta_rand  a uniform random direction              <- null
## A result that appears in all three is an artefact of cell shape, not leakage.
##
## GEOMETRY IS DONE IN CANVAS PX THROUGHOUT, so the fitted mm->px affine never enters --
## CenterX/Y_global_px (18_script_centroids.sh) and the transcript file's x/y_global_px are
## the same space. Angles are computed in a right-handed PHYSICAL frame (y negated, since
## canvas px runs opposite to slide mm on both axes) so that "+y is toward the donor" is a
## true rotation of the tissue and not a mirror of it. See CLAUDE.md, THE AXIS RULE.
##
## Usage: sbatch --mem=120G --time=2:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/18_script_select_cells.R IB1
## ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(RANN); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/lineage.R")

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"
OUT    <- "output/18_leakage"
CENT_F <- file.path(OUT, sprintf("%s_centroids.csv.gz", SAMPLE))
CLU_F  <- file.path(OUT, sprintf("%s_clusters.csv.gz",  SAMPLE))
CT_F   <- "output/16_celltypes/cell_celltypes.csv.gz"
META_F <- "output/16_regions/cells_meta_all.csv.gz"
CNT_F  <- "output/16_regions/panel_counts_all.Rds"
for (f in c(CENT_F, CLU_F, CT_F, META_F, CNT_F))
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
sc <- function(n, d) { fwrite(d, file.path(OUT, n)); cat("wrote ", file.path(OUT, n),
                                                         " (", nrow(d), " rows)\n", sep = "") }

UMPX      <- 120.280945e-6 * 1000   # um per canvas px; Area/Area.um2 confirm this exactly
## MIN_DONOR is 2, not 3, because the panel is sparse where it matters. At 3, only 38,706
## of 1,042,237 non-defect cells qualify as B donors and the median nearest-B-donor sits
## 21.6 um away -- two and a half cell diameters, which is not adjacency and would leave the
## B, T, endo and PC panels of the final figure with a few hundred cells each.
MIN_DONOR <- 2L    # transcripts of family F before a cell counts as an F donor
## D_MAX_UM only sets the `adjacent` FLAG. donorA_dist is stored per cell, so the adjacency
## threshold is re-choosable at analysis time; nothing is discarded by it here.
D_MAX_UM  <- 20    # mean IB1 cell is 62.6 um2 -> 8.9 um diameter, so ~2.2 diameters
## Transcripts are extracted for every focal cell regardless of adjacency, plus donors out
## to this radius. Extraction is the one irreversible step (a streaming pass over 13 GB), so
## it is deliberately more generous than any threshold the analysis will use.
D_EXTRACT <- 25
SEED      <- 1L
set.seed(SEED)

## pearson_clusters_low id -> the argmax lineage vocabulary. Labels from the step-09b audit
## (output/09b_marker_relabel/IB1/cluster_labels_final.csv), NEVER from HierAnnot.
##   0 low-information (defect_dominated)  dropped
##   6 pDC                                 dropped: no counterpart in LINEAGE_SETS_COARSE
CLU_MAP <- c("1" = "myeloid",   # Macrophage  confident
             "2" = "fibro",     # Stromal/FDC confident
             "3" = "T",         # T cell      confident
             "4" = "B",         # B cell      confident
             "5" = "endo",      # Endothelial confident
             "7" = "myeloid",   # cDC         ambiguous
             "8" = "PC")        # Plasma cell confident
CLU_LABEL <- c("1" = "Macrophage", "2" = "Stromal/FDC", "3" = "T cell", "4" = "B cell",
               "5" = "Endothelial", "7" = "cDC", "8" = "Plasma cell")
FAMS <- names(LINEAGE_SETS_COARSE)

cat("============================================================\n")
cat("step 18 -- transcript leakage into argmax, focal-cell selection: ", SAMPLE, "\n", sep = "")
cat("============================================================\n")

## ---- 1. cells: identity, arbiter, geometry ---------------------------------
md <- fread(META_F, select = c("cell_id","samp","fov_id","nCount_RNA","is_defect",
                               "tier_strict","tier_matched","qc_passing_all"))[samp == SAMPLE]
ct <- fread(CT_F,   select = c("cell_id","samp","celltype","lineage16"))[samp == SAMPLE]
cl <- fread(CLU_F,  select = c("cell_id","pearson_clusters_low","pearson_clusters_high"),
            colClasses = "character")
ce <- fread(CENT_F)
cat(sprintf("meta %d | celltypes %d | clusters %d | centroids %d\n",
            nrow(md), nrow(ct), nrow(cl), nrow(ce)))

## nCount_RNA is in BOTH cells_meta_all and the vendor metadata, and a silent merge leaves
## nCount_RNA.x / nCount_RNA.y with neither named nCount_RNA -- which fails later as a NULL
## rather than as an error. Drop the duplicate here and assert no others survive.
ce_cols <- setdiff(names(ce), "nCount_RNA")
cells <- Reduce(function(a, b) merge(a, b, by = "cell_id", all.x = TRUE, sort = FALSE),
                list(md, ct[, .(cell_id, celltype, lineage16)], cl, ce[, ..ce_cols]))
if (nrow(cells) != nrow(md)) stop("join changed the row count", call. = FALSE)
dup <- grep("\\.(x|y)$", names(cells), value = TRUE)
if (length(dup)) stop("merge produced duplicated columns: ", paste(dup, collapse = ", "), call. = FALSE)
if (anyNA(cells$cx_px)) stop(sum(is.na(cells$cx_px)), " cells have no centroid", call. = FALSE)

## Physical, right-handed, um. Canvas px runs opposite to slide mm on BOTH axes, so negating
## y here (and leaving x) turns the px frame into a true rotation of the tissue rather than a
## mirror of it. Every angle and every rotated coordinate downstream is in THIS frame.
cells[, `:=`(X_um = cx_px * UMPX, Y_um = -cy_px * UMPX)]
cells[, arb := CLU_MAP[pearson_clusters_low]]
cells[, arb_label := CLU_LABEL[pearson_clusters_low]]

## ---- 2. per-cell marker load per family ------------------------------------
## Needed twice: as k_A for the focal cell (how much evidence argmax actually had) and as the
## donor criterion (a cell can only donate markers it carries).
M <- readRDS(CNT_F)
gi <- match(cells$cell_id, colnames(M))
if (anyNA(gi)) stop(sum(is.na(gi)), " cells absent from the count matrix", call. = FALSE)
M <- M[, gi, drop = FALSE]
cat(sprintf("count matrix subset: %d genes x %d cells\n", nrow(M), ncol(M)))
for (f in FAMS) {
  g <- intersect(LINEAGE_SETS_COARSE[[f]], rownames(M))
  if (!length(g)) stop("no panel genes for family ", f, call. = FALSE)
  cells[[paste0("k_", f)]] <- as.integer(Matrix::colSums(M[g, , drop = FALSE]))
  cat(sprintf("  %-8s %d/%d genes on panel: %s\n", f, length(g),
              length(LINEAGE_SETS_COARSE[[f]]), paste(g, collapse = ",")))
}
rm(M); invisible(gc())

## ---- 3. focal set ----------------------------------------------------------
## Focal = a real argmax lineage AND a usable arbiter label. `ambiguous`, `all_markers_0`
## and `below_background` are not calls and are excluded; defect cells are dropped per cell,
## never via the FOV-level flag.
cells[, has_call := lineage16 %in% FAMS]
foc <- cells[has_call == TRUE & !is.na(arb) & is_defect == FALSE]
foc[, `:=`(A = lineage16, agree = lineage16 == arb)]
foc[, scoreable := A != "NK"]   # pearson_clusters_low has no NK cluster
cat(sprintf("\nfocal cells: %d (%d scoreable, %d NK carried but unscoreable)\n",
            nrow(foc), foc[, sum(scoreable)], foc[, sum(!scoreable)]))
cat("\n=== focal set: argmax x arbiter ===\n")
print(dcast(foc[, .N, by = .(A, arb)], A ~ arb, value.var = "N", fill = 0L))
cat("\n=== agreement by argmax call ===\n")
print(foc[, .(cells = .N, agree = sum(agree), pct = round(100*mean(agree), 1)), by = A][order(-cells)])

## ---- 4. donor geometry -----------------------------------------------------
## Donor pool for family F: all non-defect cells with k_F >= MIN_DONOR. Searched over the
## FULL 2.5M-cell sample, not the labelled subset, so "nearest" means nearest.
pool <- cells[is_defect == FALSE]
cat(sprintf("\ndonor pool (non-defect): %d cells\n", nrow(pool)))
nn_family <- function(qry, fam) {
  d <- pool[get(paste0("k_", fam)) >= MIN_DONOR]
  if (nrow(d) < 2L) return(NULL)
  cat(sprintf("  %-8s donors %7d  ", fam, nrow(d)))
  r <- RANN::nn2(d[, .(X_um, Y_um)], qry[, .(X_um, Y_um)], k = 2L)
  ## k=2 then drop self: a focal cell can be its own family's donor.
  self <- d$cell_id[r$nn.idx[, 1]] == qry$cell_id
  idx  <- ifelse(self, r$nn.idx[, 2], r$nn.idx[, 1])
  dst  <- ifelse(self, r$nn.dists[, 2], r$nn.dists[, 1])
  cat(sprintf("median nn %.1f um, within %g um: %.1f%%\n",
              median(dst), D_MAX_UM, 100*mean(dst <= D_MAX_UM)))
  list(id = d$cell_id[idx], dist = dst, X = d$X_um[idx], Y = d$Y_um[idx], k = d[[paste0("k_", fam)]][idx],
       n = d$nCount_RNA[idx])
}
cat("\n=== nearest-donor search, per family ===\n")
NNF <- setNames(lapply(FAMS, function(f) nn_family(foc, f)), FAMS)

## An absent field must be an error naming the field, not a length-zero replacement three
## frames down -- that is exactly how the nCount_RNA merge collision above presented.
for (f in FAMS) {
  if (is.null(NNF[[f]])) stop("no donor search result for family ", f, call. = FALSE)
  for (fl in c("id","dist","X","Y","k","n"))
    if (is.null(NNF[[f]][[fl]])) stop("NNF[['", f, "']][['", fl, "']] is NULL", call. = FALSE)
}
pick <- function(fam_vec, field, mode = c("num","chr")) {
  mode <- match.arg(mode)
  out <- if (mode == "num") rep(NA_real_, length(fam_vec)) else rep(NA_character_, length(fam_vec))
  for (f in FAMS) { i <- which(fam_vec == f); if (!length(i)) next
    out[i] <- NNF[[f]][[field]][i] }
  out
}
foc[, `:=`(donorA_id   = pick(A, "id", "chr"), donorA_dist = pick(A, "dist"),
           donorA_X    = pick(A, "X"),    donorA_Y    = pick(A, "Y"),
           donorA_k    = pick(A, "k"),    donorA_nCount = pick(A, "n"))]

## Control family C: not the argmax call, not the arbiter's call. Chosen per cell as the
## family whose donor is CLOSEST, so the control direction is matched to the test direction
## on adjacency -- the one covariate that could otherwise explain a difference between them.
ctrl_cands <- lapply(FAMS, function(f) {
  d <- NNF[[f]]$dist; d[foc$A == f | foc$arb == f] <- Inf; d })
ctrl_mat <- do.call(cbind, ctrl_cands); colnames(ctrl_mat) <- FAMS
ci <- max.col(-ctrl_mat, ties.method = "first")
foc[, C := FAMS[ci]]
foc[, `:=`(donorC_id = pick(C, "id", "chr"), donorC_dist = pick(C, "dist"),
           donorC_X  = pick(C, "X"),  donorC_Y    = pick(C, "Y"), donorC_k = pick(C, "k"))]

## Arbiter-direction donor: toward the nearest cell rich in the markers of the type the
## CLUSTER says this cell is. Under leakage this direction should show nothing special;
## under "the cluster is wrong and argmax is right" it is where the signal would sit.
foc[, `:=`(donorB_id = pick(arb, "id", "chr"), donorB_dist = pick(arb, "dist"),
           donorB_X  = pick(arb, "X"),  donorB_Y    = pick(arb, "Y"))]

ang <- function(dx, dy) atan2(dy, dx)
foc[, `:=`(theta_A    = ang(donorA_X - X_um, donorA_Y - Y_um),
           theta_C    = ang(donorC_X - X_um, donorC_Y - Y_um),
           theta_B    = ang(donorB_X - X_um, donorB_Y - Y_um),
           theta_rand = runif(.N, -pi, pi))]
foc[, r_eff_um := sqrt(area_um2 / pi)]
foc[, k_A := as.integer(mapply(function(f, i) foc[[paste0("k_", f)]][i], A, seq_len(.N)))]
foc[, adjacent := is.finite(donorA_dist) & donorA_dist <= D_MAX_UM]

cat(sprintf("\nfocal cells with an A-donor within %g um: %d / %d (%.1f%%)\n",
            D_MAX_UM, foc[, sum(adjacent)], nrow(foc), 100*foc[, mean(adjacent)]))
cat("\n=== adjacency and donor load, by argmax call x agreement ===\n")
print(foc[, .(cells = .N, adj_pct = round(100*mean(adjacent), 1),
              med_donorA_dist = round(median(donorA_dist), 1),
              med_k_A = as.numeric(median(k_A)), med_nCount = as.numeric(median(nCount_RNA)),
              med_area = round(median(area_um2), 1)), by = .(A, agree)][order(A, agree)])

## Angular separation between the test and control directions. If they are collinear the
## control cannot discriminate, so this has to be inspected, not assumed.
foc[, sep_AC_deg := abs(((theta_A - theta_C + pi) %% (2*pi)) - pi) * 180/pi]
cat(sprintf("\nangle between theta_A and theta_C: median %.0f deg, %% under 30 deg: %.1f\n",
            median(foc$sep_AC_deg, na.rm = TRUE), 100*mean(foc$sep_AC_deg < 30, na.rm = TRUE)))

KEEP <- c("cell_id","fov_id","A","celltype","arb","arb_label","pearson_clusters_low",
          "pearson_clusters_high","agree","scoreable","adjacent",
          "nCount_RNA","nFeature_RNA","area_um2","r_eff_um","aspect","circularity",
          "unassigned_tx","qc_passing_all","tier_matched","cx_px","cy_px","X_um","Y_um",
          "k_A","C", paste0("k_", FAMS),
          "donorA_id","donorA_dist","donorA_k","donorA_nCount","theta_A",
          "donorC_id","donorC_dist","donorC_k","theta_C",
          "donorB_id","donorB_dist","theta_B","theta_rand","sep_AC_deg")
sc(sprintf("%s_focal_cells.csv.gz", SAMPLE), foc[, ..KEEP])

## The cell-id list the transcript extraction filters on. Focal cells AND their donors --
## the donor's own transcripts are wanted so a donor-side depletion can be checked later.
ids <- unique(c(foc$cell_id, foc$donorA_id[is.finite(foc$donorA_dist) &
                                           foc$donorA_dist <= D_EXTRACT]))
ids <- ids[!is.na(ids)]
writeLines(ids, file.path(OUT, sprintf("%s_extract_ids.txt", SAMPLE)))
cat(sprintf("\nwrote %s (%d cell ids: %d focal + %d donors)\n",
            file.path(OUT, sprintf("%s_extract_ids.txt", SAMPLE)), length(ids),
            nrow(foc), length(ids) - nrow(foc)))
cat("\ndone\n")
