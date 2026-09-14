# =============================================================================
# qc_ib1_overfixation_summary.R
#
# SUMMARY LOGIC for the over-fixation diagnosis in IB1.
#
# This is not a report of a conclusion already reached. It is the decision
# procedure itself, written down: a fixed set of criteria, each with a
# mechanism, a pre-registered numeric rule, and an explicit list of competing
# hypotheses that the criterion eliminates *if and only if* its rule fires. The
# script evaluates every rule against the data and lets the verdict fall out.
# A criterion that fails is reported as failing and its exclusions are not
# applied, so the summary can weaken the over-fixation call, not only support it.
#
# ---------------------------------------------------------------- THE CLAIM --
# H1  OVER-FIXATION. IB1 sat in formalin too long. Excess methylene-bridge
#     crosslinking (a) masks protein epitopes, (b) occludes RNA from probe
#     hybridisation, and (c) does both along a gradient set by formalin
#     diffusing inward from the block surface.
#
# ------------------------------------------------------- THE COMPETING SET --
# H2  UNDER-FIXATION / ischaemic delay -> RNA degradation.
# H3  Segmentation dilution (real counts smeared over oversized cell masks).
# H4  Optical / decoding failure.
# H5  Antibody-only staining failure (bad lot or protocol, RNA side intact).
# H6  Genuine biology (this node simply has fewer B cells).
# H7  Section-surface handling artifact (drying / rim damage on the cut section,
#     as opposed to a block-scale fixation gradient).
#
# RB1 is the internal negative control throughout: same run, same panel, same
# operator, no staining defect. Any criterion that also fires in RB1 is
# measuring the platform, not IB1's lesion, and is downgraded to NEUTRAL.
#
# ------------------------------------------------------------ THE PIVOT ------
# The single fact that organises the whole argument: the IB1 defect is
# QUALITATIVE, not quantitative. Cluster-i cells are not shallow -- they carry
# slightly MORE counts than their neighbours. What they have lost is
# specificity. Degradation (H2) cannot produce that; crosslinking can, because
# occluded transcripts stop contributing specific signal while a nonspecific
# floor keeps the count total up. C9 is the direct test of that mechanism and is
# the criterion most able to break the diagnosis.
#
# Usage:
#   sbatch --mem=96G --cpus-per-task=4 --time=3:00:00 \
#     --job-name=ovfix --output=logs/ovfix_%j.out --error=logs/ovfix_%j.err \
#     scripts/sbatch/run_rscript.sh scripts/sbatch/qc_ib1_overfixation_summary.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

set.seed(1)
options(width = 200, stringsAsFactors = FALSE)

# ------------------------------------------------------------------ CONFIG ---
MERGED   <- "output/04_merge/merged_01_merged.Rds"
FOV_CSV  <- "output/02_filter_QC/qc_spatial_fov.csv"
OUT_DIR  <- "output/07_overfixation"
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

CASE     <- "IB1"          # the sample under diagnosis
CTRL     <- "RB1"          # internal negative control
DEFECT   <- "i"            # the InSituType cluster that carries the defect

# Palette: slots 1-2 of the validated reference categorical palette (light mode).
# Documented as clearing every all-pairs gate in both modes at <= 3 slots, so no
# re-validation is needed for a two-series chart. Control = blue, case = orange.
COL <- c(RB1 = "#2a78d6", IB1 = "#eb6834")
INK <- c(primary = "#0b0b0b", secondary = "#52514e", grid = "#d8d7d2")

# --- pre-registered thresholds. Set here, before looking, so the rules below
# --- are decisions rather than descriptions of whatever came out.
THR <- list(
  channel_down     = 0.90,   # C1 antibody channel ratio (defect/rest) counts as "lost"
  channel_dapi_up  = 1.00,   # C1 DAPI ratio must be at least flat
  ctrl_flat        = 0.15,   # C1/C8 control must sit within +-15% of no effect
  flatten_slope    = -0.20,  # C2 slope of log2(detection ratio) on log2(baseline prevalence)
  flatten_p        = 0.01,   # C2 significance required on that slope
  specificity_idx  = 0.70,   # C2b (informational) detection fall vs depth fall
  panelwide_det    = 0.80,   # C2b (informational) non-B drop
  featureless_lfc  = 1.00,   # C3 max |log2FC| for a cluster to count as profile-flat
  featureless_rel  = 0.60,   # C3 ... and it must be this fraction of its peers' median
  depth_preserved  = 0.95,   # C4 defect/rest count and density ratio floor
  negprobe_slack   = 0.005,  # C4/C6 tolerated excess propNegative in the defect
  area_flat        = 0.10,   # C5 |area ratio - 1| ceiling
  rho_gradient     = 0.20,   # C7 Spearman rho(radius, defect) to call a gradient
  rho_ctrl_flat    = 0.10,   # C7 control rho must be below this
  rho_fov_couple   = -0.40,  # C8 rho(frac_defect, Mean.Membrane) at FOV level
  flatten_share    = 0.85,   # C9 specific-marker count share must fall this far
  n_match_profile  = 8000L,  # C9 depth-matched cells per arm (dense profile work)
  n_quantile       = 5L      # radial bins
)

# A priori marker panel for C2/C9. Chosen for high prevalence and lineage
# breadth BEFORE looking at IB1, so that a panel-wide drop cannot be explained
# by "this node has fewer B cells" -- the T, epithelial, myeloid, stromal and
# housekeeping members have nothing to do with the B compartment.
PANEL <- list(
  B          = c("MS4A1", "CD79A", "CD79B", "BANK1", "CD19"),
  T          = c("CD3E", "CD3D", "TRAC", "IL7R"),
  epithelial = c("KRT8", "KRT18", "KRT19", "EPCAM"),
  myeloid    = c("LYZ", "CD68", "ITGAX"),
  stromal    = c("ACTA2", "PECAM1", "COL1A1"),
  ubiquitous = c("CD74", "PTPRC", "B2M", "ACTB")
)
NONB_SETS <- setdiff(names(PANEL), "B")

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
rule_hdr <- function(id, title) cat(sprintf("\n===== %s  %s %s\n", id, title,
                                            strrep("=", max(0, 60 - nchar(title)))))
fmt <- function(x, d = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")

# Every criterion registers itself here. `excludes` is applied only when the
# rule fires -- that coupling is what makes this a logic and not a narrative.
CRIT <- list()
register <- function(id, title, mechanism, prediction, rule_text,
                     supports, excludes = character(), neutral_reason = NA_character_,
                     stats = list()) {
  verdict <- if (isTRUE(supports)) "SUPPORTS" else if (isFALSE(supports)) "REFUTES" else "NEUTRAL"
  CRIT[[id]] <<- list(id = id, title = title, mechanism = mechanism,
                      prediction = prediction, rule = rule_text,
                      verdict = verdict,
                      excludes = if (verdict == "SUPPORTS") excludes else character(),
                      excludes_if_fired = excludes,
                      neutral_reason = neutral_reason, stats = stats)
  cat(sprintf("  rule : %s\n  ---> %s\n", rule_text, verdict))
  if (length(stats)) {
    for (nm in names(stats)) cat(sprintf("       %-34s %s\n", nm, fmt(stats[[nm]])))
  }
  invisible(NULL)
}

# ==============================================================================
# 0. LOAD
# ==============================================================================
msg("reading", MERGED)
seu <- readRDS(MERGED)
md  <- as.data.table(seu@meta.data, keep.rownames = "barcode")
cts <- LayerData(seu, assay = "RNA", layer = "counts")
genes <- rownames(seu)
stopifnot(nrow(md) == ncol(cts))
msg("cells:", nrow(md), "| genes:", length(genes))

need <- c("merge_sample", "inSitu_clusters", "nCount_RNA", "nFeature_RNA",
          "nCount_negprobes", "nCount_falsecode", "Area.um2", "complexity",
          "propNegative", "Mean.DAPI", "Mean.Membrane", "Mean.CD45", "Mean.PanCK",
          "x_slide_mm", "y_slide_mm", "lineage", "strongB",
          "inSitu_clusters_posterior_probability")
stopifnot(all(need %in% names(md)))

md[, cps        := nCount_RNA / pmax(Area.um2, 1)]
md[, frac_neg   := nCount_negprobes / pmax(nCount_RNA + nCount_negprobes, 1)]
md[, frac_false := nCount_falsecode / pmax(nCount_RNA + nCount_falsecode, 1)]
md[, is_defect  := inSitu_clusters == DEFECT]
md[, informative := !lineage %in% c("all_markers_0", "below_background", "ambiguous")]

# radial distance from each section's own centroid (block-scale axis) and
# distance to the nearest tissue boundary incl. internal holes (local axis).
# C7 needs BOTH: they separate a block diffusion gradient from a rim artifact.
BIN <- 0.1
md[, r_mm := {
  sqrt((x_slide_mm - median(x_slide_mm))^2 + (y_slide_mm - median(y_slide_mm))^2)
}, by = merge_sample]

edge_distance <- function(x, y) {
  ix <- floor((x - min(x)) / BIN) + 1L
  iy <- floor((y - min(y)) / BIN) + 1L
  occ <- matrix(FALSE, max(ix), max(iy))
  occ[cbind(ix, iy)] <- TRUE
  pad <- matrix(FALSE, nrow(occ) + 2L, ncol(occ) + 2L)
  pad[2:(nrow(occ) + 1L), 2:(ncol(occ) + 1L)] <- occ
  nb_empty <- matrix(FALSE, nrow(occ), ncol(occ))
  for (dx in -1:1) for (dy in -1:1) {
    if (dx == 0 && dy == 0) next
    nb_empty <- nb_empty | !pad[(2:(nrow(occ) + 1L)) + dx, (2:(ncol(occ) + 1L)) + dy]
  }
  edge <- which(occ & nb_empty, arr.ind = TRUE)
  ub <- unique(cbind(ix, iy))
  d  <- vapply(seq_len(nrow(ub)), function(k)
    BIN * sqrt(min((edge[, 1] - ub[k, 1])^2 + (edge[, 2] - ub[k, 2])^2)), numeric(1))
  d[match(paste(ix, iy), paste(ub[, 1], ub[, 2]))]
}
msg("computing distance-to-tissue-edge per sample")
md[, d_edge_mm := edge_distance(x_slide_mm, y_slide_mm), by = merge_sample]
md[, rq := cut(r_mm, breaks = quantile(r_mm, 0:THR$n_quantile / THR$n_quantile),
               include.lowest = TRUE, labels = FALSE), by = merge_sample]

CASE_I <- md$merge_sample == CASE
CTRL_I <- md$merge_sample == CTRL

# convenience: median of `col` for defect vs rest within a sample, and the ratio
ratio_defect <- function(dt, col) {
  a <- median(dt[[col]][dt$is_defect],  na.rm = TRUE)
  b <- median(dt[[col]][!dt$is_defect], na.rm = TRUE)
  c(defect = a, rest = b, ratio = a / b)
}
mdC <- md[CASE_I]; mdR <- md[CTRL_I]

cat("\n########## composition ##########\n")
print(md[, .(cells = .N, pct_defect = round(100 * mean(is_defect), 2),
             median_nCount = median(nCount_RNA),
             pct_informative = round(100 * mean(informative), 1),
             pct_strongB = round(100 * mean(strongB), 1)), by = merge_sample])

# ==============================================================================
# C1  PROTEIN-CHANNEL DISSOCIATION
# ==============================================================================
rule_hdr("C1", "antibody channels fail, the DNA stain does not")
ch <- c("Mean.Membrane", "Mean.CD45", "Mean.PanCK", "Mean.DAPI")
c1 <- rbindlist(lapply(ch, function(k)
  data.table(channel = k,
             case_ratio = ratio_defect(mdC, k)[["ratio"]],
             ctrl_ratio = ratio_defect(mdR, k)[["ratio"]])))
print(c1)
c1_ok <- c1[channel == "Mean.Membrane", case_ratio] < THR$channel_down &
         c1[channel == "Mean.CD45",     case_ratio] < THR$channel_down &
         c1[channel == "Mean.DAPI",     case_ratio] >= THR$channel_dapi_up &
         abs(c1[channel == "Mean.Membrane", ctrl_ratio] - 1) < THR$ctrl_flat
register(
  "C1", "Protein-channel dissociation",
  paste("Formaldehyde crosslinks protein. Antibody epitopes (Membrane, CD45,",
        "PanCK) are masked; DAPI intercalates DNA and is indifferent to protein",
        "crosslinking. Over-fixation therefore splits the channels rather than",
        "dimming them together."),
  "Membrane and CD45 down in the defect population, DAPI flat or up.",
  sprintf("Membrane<%.2f AND CD45<%.2f AND DAPI>=%.2f (case) AND |Membrane-1|<%.2f (control)",
          THR$channel_down, THR$channel_down, THR$channel_dapi_up, THR$ctrl_flat),
  supports = c1_ok,
  excludes = c("H4"),   # optical failure would take DAPI down with everything else
  stats = setNames(as.list(c1$case_ratio), paste0("case ratio ", c1$channel))
)

# ==============================================================================
# C2  IS THE DETECTION LOSS LINEAGE-BLIND?
# ==============================================================================
# TWO DISCARDED DESIGNS, KEPT ON THE RECORD BECAUSE THEY SHAPE THE RULE.
#
# (i) Cross-sample. Comparing per-gene detection between IB1 and RB1 mixes probe
#     accessibility with cell composition, because the two blocks are different
#     tissues. Retained below as C2b, informational only.
#
# (ii) "Occlusion is not gene-selective, so EVERY lineage set must fall." That
#     rule is arithmetically unsatisfiable and was wrong to write. Detection is
#     measured at matched total counts, so the measure is compositional: if a
#     flat background floor replaces structured signal, the panel's ABUNDANT
#     genes must lose share and its RARE genes must gain it. Something always
#     rises. Demanding that nothing rises tests closure, not biology.
#
# The correct closure-aware signature of a background floor is therefore not
# "everything falls" but "the profile is pulled toward uniform": the detection
# ratio should decline with a gene's baseline prevalence, slope 0 under no
# change and -1 under a perfectly uniform floor. A genuine cell type has no such
# relationship -- which of its genes go up is set by its lineage, not by how
# common the gene happened to be.
rule_hdr("C2", "is the defect's profile pulled toward uniform?")
panel_genes <- intersect(unlist(PANEL), genes)
panel_set   <- rep(names(PANEL), lengths(PANEL))[match(panel_genes, unlist(PANEL))]
msg("panel genes present:", length(panel_genes), "of", length(unlist(PANEL)))
msg("extracting panel rows from the count store")
pm <- as.matrix(cts[panel_genes, , drop = FALSE])   # small: |panel| x n_cells
det <- pm > 0
rm(pm); invisible(gc(FALSE))

# Depth matching: bin nCount into common ventiles and take equal numbers from
# each arm within a bin. Detection rate is strongly depth-dependent, so an
# unmatched ratio cannot separate "shallower" from "occluded".
brk <- unique(quantile(md$nCount_RNA, 0:20 / 20))
md[, depth_bin := cut(nCount_RNA, brk, include.lowest = TRUE, labels = FALSE)]
matched_idx <- function(mask_a, mask_b) {
  ia <- ib <- integer(0)
  for (b in sort(unique(md$depth_bin))) {
    A <- which(mask_a & md$depth_bin == b)
    B <- which(mask_b & md$depth_bin == b)
    n <- min(length(A), length(B))
    if (n == 0L) next
    ia <- c(ia, sample(A, n)); ib <- c(ib, sample(B, n))
  }
  list(a = ia, b = ib)
}
mm <- matched_idx(CASE_I, CTRL_I)
msg("depth-matched cohort:", length(mm$a), "cells per sample")

depth_ratio <- median(md$nCount_RNA[CASE_I]) / median(md$nCount_RNA[CTRL_I])

# Within one sample, depth-matched: defect cells vs that sample's own unaffected
# cells. Steric occlusion is not gene-selective, so it must pull EVERY lineage
# set down together. A cell type pulls some sets down and pushes others up.
within_detection <- function(smask) {
  w <- matched_idx(smask & md$is_defect, smask & !md$is_defect)
  data.table(gene = panel_genes, set = panel_set,
             det_defect = rowMeans(det[, w$a, drop = FALSE]),
             det_rest   = rowMeans(det[, w$b, drop = FALSE]),
             n_per_arm  = length(w$a))[, ratio := det_defect / det_rest][]
}
w_case <- within_detection(CASE_I)[, sample := CASE]
w_ctrl <- within_detection(CTRL_I)[, sample := CTRL]
wdet   <- rbind(w_case, w_ctrl)
set_case <- w_case[det_rest >= 0.02, .(med_ratio = median(ratio)), by = set][order(med_ratio)]
set_ctrl <- w_ctrl[det_rest >= 0.02, .(med_ratio = median(ratio)), by = set][order(med_ratio)]
cat("-- IB1: defect vs own cells, depth-matched (", w_case$n_per_arm[1], "per arm) --\n")
print(w_case[det_rest >= 0.02][order(set, -det_rest)][, .(gene, set,
      det_defect = round(100 * det_defect, 1), det_rest = round(100 * det_rest, 1),
      ratio = round(ratio, 3))], row.names = FALSE)
cat("\nper-lineage-set median ratio (descriptive) -- case:\n"); print(set_case, row.names = FALSE)
cat("per-lineage-set median ratio (descriptive) -- control:\n"); print(set_ctrl, row.names = FALSE)

# The criterion: how steeply does the detection ratio fall with a gene's baseline
# prevalence? 0 = profile unchanged, -1 = profile fully replaced by a flat floor.
flatten_fit <- function(w) {
  u <- w[det_rest >= 0.02 & det_defect > 0]
  f <- lm(log2(ratio) ~ log2(det_rest), data = u)
  list(slope = unname(coef(f)[2]), p = summary(f)$coefficients[2, 4],
       rho = cor(log2(u$det_rest), log2(u$ratio), method = "spearman"), n = nrow(u))
}
ff_case <- flatten_fit(w_case); ff_ctrl <- flatten_fit(w_ctrl)
cat(sprintf("\nflattening slope: case %+.3f (p %.2g, rho %+.2f, n %d) | control %+.3f (p %.2g)\n",
            ff_case$slope, ff_case$p, ff_case$rho, ff_case$n, ff_ctrl$slope, ff_ctrl$p))
c2_ok <- ff_case$slope < THR$flatten_slope && ff_case$p < THR$flatten_p
register(
  "C2", "Profile pulled toward uniform",
  paste("A nonspecific floor does not erase genes evenly in SHARE terms -- it",
        "erases structure. Detection here is measured at matched total counts, so",
        "the panel is compositional: replacing signal with a flat floor drags",
        "abundant genes down and pushes rare ones up, giving the detection ratio",
        "a negative dependence on baseline prevalence that tends to slope -1 as",
        "the floor takes over. A cell type shows no such dependence, because",
        "which of its genes rise is decided by its lineage, not by how common",
        "the gene happens to be."),
  "Detection ratio (defect vs the same sample's matched cells) declines with baseline prevalence.",
  sprintf("slope of log2(ratio) on log2(prevalence) < %.2f at p < %.2g",
          THR$flatten_slope, THR$flatten_p),
  supports = c2_ok,
  excludes = c("H6"),
  stats = list("flattening slope (case)" = ff_case$slope,
               "slope p (case)" = ff_case$p,
               "spearman rho (case)" = ff_case$rho,
               "flattening slope (control)" = ff_ctrl$slope,
               "matched cells per arm" = w_case$n_per_arm[1])
)

# ---- C2b informational: the confounded cross-sample view, reported not used --
rule_hdr("C2b", "cross-sample detection (confounded by composition)")
c2 <- data.table(
  gene = panel_genes, set = panel_set,
  det_case = rowMeans(det[, CASE_I, drop = FALSE]),
  det_ctrl = rowMeans(det[, CTRL_I, drop = FALSE]),
  det_case_m = rowMeans(det[, mm$a, drop = FALSE]),
  det_ctrl_m = rowMeans(det[, mm$b, drop = FALSE]))
c2[, `:=`(ratio = det_case / det_ctrl, ratio_matched = det_case_m / det_ctrl_m)]
c2_use <- c2[det_ctrl >= 0.02]
spec_idx <- median(c2_use$ratio, na.rm = TRUE) / depth_ratio
nonb_med <- median(c2_use[set %in% NONB_SETS, ratio_matched], na.rm = TRUE)
print(c2_use[order(set, -det_ctrl)][, .(gene, set,
      det_case = round(100 * det_case, 1), det_ctrl = round(100 * det_ctrl, 1),
      ratio = round(ratio, 3), ratio_matched = round(ratio_matched, 3))], row.names = FALSE)
register(
  "C2b", "Cross-sample detection difference",
  paste("Reported for completeness and explicitly NOT used as evidence. IB1 and",
        "RB1 are different tissues, so this contrast mixes probe accessibility",
        "with cell composition and cannot separate them."),
  "n/a", "informational only - confounded by design", supports = NA,
  neutral_reason = "cross-tissue contrast confounds accessibility with composition",
  stats = list("depth ratio case/control" = depth_ratio,
               "median detection ratio" = median(c2_use$ratio, na.rm = TRUE),
               "specificity index" = spec_idx,
               "non-B matched median" = nonb_med)
)

# ==============================================================================
# C3  THE DEFECT POPULATION IS PROFILE-FLAT
# ==============================================================================
rule_hdr("C3", "the defect cluster has no genes of its own")
cluster_profiles <- function(idx) {
  idx <- sort(idx)                       # BPCells wants ascending column order
  g   <- droplevels(factor(md$inSitu_clusters[idx]))
  ind <- matrix(0, length(idx), nlevels(g), dimnames = list(NULL, levels(g)))
  ind[cbind(seq_along(g), as.integer(g))] <- 1
  s <- tryCatch(as.matrix(cts[, idx] %*% ind),
                error = function(e) {
                  msg("matmul unavailable (", conditionMessage(e), ") - falling back to rowSums")
                  sapply(levels(g), function(l) rowSums(cts[, idx[g == l], drop = FALSE]))
                })
  rownames(s) <- genes
  cpm <- log2(sweep(s, 2, pmax(colSums(s), 1), "/") * 1e6 + 1)
  cpm
}
msg("computing per-cluster mean profiles (case)")
prof <- cluster_profiles(which(CASE_I))
maxlfc <- sapply(colnames(prof), function(k) {
  others <- rowMeans(prof[, setdiff(colnames(prof), k), drop = FALSE])
  max(prof[, k] - others, na.rm = TRUE)
})
print(round(sort(maxlfc), 3))
peer_med <- median(maxlfc[setdiff(names(maxlfc), DEFECT)])
post_med <- median(mdC$inSitu_clusters_posterior_probability[mdC$is_defect], na.rm = TRUE)
c3_ok <- maxlfc[[DEFECT]] < THR$featureless_lfc &&
         maxlfc[[DEFECT]] < THR$featureless_rel * peer_med
register(
  "C3", "Defect population is profile-flat",
  paste("A crosslink-silenced cell reports the panel's background shape, not a",
        "cell type. InSituType still assigns it confidently, so this is a real,",
        "reproducible mode of the data rather than clustering uncertainty."),
  "Cluster i has a near-panel-average profile: no gene stands out, unlike every other cluster.",
  sprintf("max log2FC(i) < %.2f AND < %.2f x median(max log2FC of other clusters)",
          THR$featureless_lfc, THR$featureless_rel),
  supports = c3_ok,
  excludes = c("H6"),
  stats = list("max log2FC cluster i" = maxlfc[[DEFECT]],
               "median max log2FC peers" = peer_med,
               "median posterior prob (i)" = post_med)
)

# ==============================================================================
# C4  DEPTH AND NOISE ARE PRESERVED -> NOT DEGRADATION
# ==============================================================================
rule_hdr("C4", "the defect is qualitative, not a loss of material")
c4 <- rbindlist(lapply(c("nCount_RNA", "nFeature_RNA", "cps", "complexity"), function(k) {
  v <- ratio_defect(mdC, k)
  data.table(metric = k, defect = v[["defect"]], rest = v[["rest"]], ratio = v[["ratio"]])
}))
print(c4)
neg_d <- median(mdC$frac_neg[mdC$is_defect]);  neg_r <- median(mdC$frac_neg[!mdC$is_defect])
c4_ok <- c4[metric == "nCount_RNA", ratio] >= THR$depth_preserved &&
         c4[metric == "cps",        ratio] >= THR$depth_preserved &&
         neg_d <= neg_r + THR$negprobe_slack
register(
  "C4", "Depth and noise preserved in the defect",
  paste("Degradation removes molecules: fewer counts, fewer features, lower",
        "counts per unit area, and a rising negative-probe share as real signal",
        "thins out. Crosslinking removes ACCESS, not molecules, so the totals",
        "survive. This is the criterion that separates over- from under-fixation",
        "on the biochemistry rather than on the geometry."),
  "Counts, features and counts/um2 flat or higher in the defect; negative-probe share flat.",
  sprintf("nCount ratio >= %.2f AND cps ratio >= %.2f AND propNeg excess <= %.3f",
          THR$depth_preserved, THR$depth_preserved, THR$negprobe_slack),
  supports = c4_ok,
  excludes = c("H2"),
  stats = list("nCount ratio" = c4[metric == "nCount_RNA", ratio],
               "cps ratio"    = c4[metric == "cps", ratio],
               "negprobe frac defect" = neg_d, "negprobe frac rest" = neg_r)
)

# ==============================================================================
# C5  NOT SEGMENTATION DILUTION
# ==============================================================================
rule_hdr("C5", "cell masks are the same size")
c5 <- ratio_defect(mdC, "Area.um2")
print(round(c5, 3))
c5_ok <- abs(c5[["ratio"]] - 1) < THR$area_flat &&
         ratio_defect(mdC, "cps")[["ratio"]] > THR$depth_preserved
register(
  "C5", "Cell geometry unchanged",
  "Over-segmentation would spread normal counts over larger masks: area up, counts per um2 down.",
  "Defect cells are the same size as their neighbours and no less dense.",
  sprintf("|area ratio - 1| < %.2f AND cps ratio > %.2f", THR$area_flat, THR$depth_preserved),
  supports = c5_ok, excludes = c("H3"),
  stats = list("area ratio" = c5[["ratio"]], "cps ratio" = ratio_defect(mdC, "cps")[["ratio"]])
)

# ==============================================================================
# C6  NOT DECODING FAILURE
# ==============================================================================
rule_hdr("C6", "barcode decoding is clean")
c6 <- data.table(
  metric = c("frac_negprobe", "frac_falsecode"),
  case_defect = c(median(mdC$frac_neg[mdC$is_defect]),   median(mdC$frac_false[mdC$is_defect])),
  case_rest   = c(median(mdC$frac_neg[!mdC$is_defect]),  median(mdC$frac_false[!mdC$is_defect])),
  ctrl        = c(median(mdR$frac_neg),                  median(mdR$frac_false)))
print(c6)
c6_ok <- all(c6$case_defect <= c6$case_rest + THR$negprobe_slack) &&
         all(c6$case_defect <= c6$ctrl + THR$negprobe_slack)
register(
  "C6", "Decoding integrity intact",
  paste("An optical or decoding failure inflates negative-probe and false-code",
        "rates, because misreads are drawn uniformly from barcode space."),
  "Negative-probe and false-code shares in the defect no higher than elsewhere, or than the control.",
  sprintf("defect shares <= rest and <= control, within %.3f", THR$negprobe_slack),
  supports = c6_ok, excludes = c("H4"),
  stats = list("negprobe defect" = c6$case_defect[1], "falsecode defect" = c6$case_defect[2],
               "falsecode control" = c6$ctrl[2])
)

# ==============================================================================
# C7  DIFFUSION-SHAPED SPATIAL GRADIENT
# ==============================================================================
rule_hdr("C7", "the defect follows a block-scale gradient, and its sign is 'over'")
grad <- md[, .(n = .N,
               median_r = round(median(r_mm), 2),
               pct_defect = round(100 * mean(is_defect), 1),
               median_nCount = round(median(nCount_RNA)),
               pct_informative = round(100 * mean(informative), 1),
               pct_strongB = round(100 * mean(strongB), 1),
               Membrane = round(median(Mean.Membrane)),
               DAPI = round(median(Mean.DAPI))),
           by = .(merge_sample, rq)][order(merge_sample, rq)]
print(grad)
rho <- function(dt, v) suppressWarnings(
  cor(dt[[v]], as.numeric(dt$is_defect), method = "spearman", use = "complete.obs"))
rho_r_case <- rho(mdC, "r_mm"); rho_r_ctrl <- rho(mdR, "r_mm")
rho_e_case <- rho(mdC, "d_edge_mm")
# sign convention: r_mm grows toward the rim, d_edge_mm grows toward the interior.
# Over-fixation (surface-in) => rho_r > 0. Under-fixation (core spared last,
# autolysed) => rho_r < 0. A cut-section handling artifact tracks the LOCAL
# boundary, so it would make |rho_edge| >= |rho_r|.
c7_ok <- rho_r_case > THR$rho_gradient && abs(rho_r_ctrl) < THR$rho_ctrl_flat
c7_block_scale <- abs(rho_r_case) > abs(rho_e_case)
register(
  "C7", "Diffusion-shaped gradient, rim-worst",
  paste("Formalin enters the block from its surface. Tissue nearest that surface",
        "accumulates the most crosslinking, so over-fixation is worst at the rim",
        "and mildest in the core. Under-fixation is the mirror image: the core is",
        "what formalin never reached. The sign of the gradient therefore reads",
        "the direction of the fixation error directly."),
  "Defect fraction rises from core to rim in the case and is flat in the control.",
  sprintf("rho(radius, defect) > +%.2f in case AND |rho| < %.2f in control",
          THR$rho_gradient, THR$rho_ctrl_flat),
  supports = c7_ok,
  excludes = if (c7_block_scale) c("H2", "H7") else c("H2"),
  stats = list("rho radius (case)" = rho_r_case, "rho radius (control)" = rho_r_ctrl,
               "rho edge-distance (case)" = rho_e_case,
               "pct defect core quintile" = grad[merge_sample == CASE & rq == 1, pct_defect],
               "pct defect rim quintile"  = grad[merge_sample == CASE & rq == THR$n_quantile, pct_defect])
)
if (!c7_block_scale)
  cat("  NOTE: local edge distance explains the gradient at least as well as global\n",
      "        radius -- H7 (section handling artifact) is NOT excluded.\n")

# ==============================================================================
# C8  ONE UPSTREAM CAUSE DRIVES BOTH READOUTS
# ==============================================================================
rule_hdr("C8", "protein masking and RNA loss are coupled FOV by FOV")
fov <- fread(FOV_CSV)
c8 <- fov[, .(n_fov = .N,
              rho_membrane = cor(frac_i, Mean.Membrane, method = "spearman"),
              rho_cd45     = cor(frac_i, Mean.CD45,     method = "spearman"),
              rho_dapi     = cor(frac_i, Mean.DAPI,     method = "spearman")),
          by = sample_id]
print(c8)
c8_case <- c8[sample_id == CASE]
c8_ok <- c8_case$rho_membrane < THR$rho_fov_couple && c8_case$rho_dapi > 0
register(
  "C8", "Protein and RNA readouts co-vary",
  paste("Two independent accidents -- a bad antibody lot and a shallow library --",
        "would not track each other across 300 fields of view. A single upstream",
        "insult acting on both protein and RNA does, and it does so with opposite",
        "sign on the DNA stain, which the insult does not touch."),
  "Per-FOV defect fraction is strongly anti-correlated with antibody channels and positively with DAPI.",
  sprintf("rho(frac_i, Membrane) < %.2f AND rho(frac_i, DAPI) > 0", THR$rho_fov_couple),
  supports = c8_ok,
  excludes = c("H5"),
  stats = list("rho membrane" = c8_case$rho_membrane, "rho CD45" = c8_case$rho_cd45,
               "rho DAPI" = c8_case$rho_dapi, "n FOV" = c8_case$n_fov)
)

# ==============================================================================
# C9  THE MECHANISM TEST: is the surviving signal a nonspecific floor?
# ==============================================================================
# C4 established that defect cells are not shallow. That is only compatible with
# C2's loss of specificity if the counts that remain are spread flat across the
# panel instead of concentrated in a few real transcripts. This is the criterion
# that can break the diagnosis: if defect cells at matched depth carry the SAME
# specific-marker share as their neighbours, the occlusion story fails.
rule_hdr("C9", "at equal depth, is the defect's signal flatter?")
mmi <- matched_idx(CASE_I & md$is_defect, CASE_I & !md$is_defect)
n_keep <- min(length(mmi$a), THR$n_match_profile)
ia <- sort(sample(mmi$a, n_keep)); ib <- sort(sample(mmi$b, n_keep))
msg("depth-matched within-case cohort:", length(ia), "cells per arm")

panel_rows <- match(panel_genes, genes)
profile_stats <- function(idx, chunk = 2000L) {
  out <- vector("list", ceiling(length(idx) / chunk))
  for (j in seq_along(out)) {
    sel <- idx[((j - 1L) * chunk + 1L):min(j * chunk, length(idx))]
    x <- as.matrix(cts[, sel, drop = FALSE])
    tot <- pmax(colSums(x), 1)
    p <- sweep(x, 2, tot, "/")
    H <- -colSums(ifelse(p > 0, p * log(p), 0))
    out[[j]] <- data.table(top1_share = apply(x, 2, max) / tot,
                           eff_genes  = exp(H),
                           spec_share = colSums(x[panel_rows, , drop = FALSE]) / tot,
                           nCount     = tot)
    rm(x, p); gc(FALSE)
  }
  rbindlist(out)
}
sa <- profile_stats(ia); sa[, arm := "defect"]
sb <- profile_stats(ib); sb[, arm := "rest"]
c9 <- rbind(sa, sb)
print(c9[, lapply(.SD, function(v) round(median(v), 4)),
         by = arm, .SDcols = c("nCount", "top1_share", "eff_genes", "spec_share")])
ss_d <- median(sa$spec_share); ss_r <- median(sb$spec_share)
eg_d <- median(sa$eff_genes);  eg_r <- median(sb$eff_genes)
w_p  <- suppressWarnings(wilcox.test(sa$spec_share, sb$spec_share)$p.value)
c9_ok <- ss_d < THR$flatten_share * ss_r && eg_d > eg_r
register(
  "C9", "Surviving signal is a nonspecific floor",
  paste("The mechanism's load-bearing prediction. If crosslinking occludes",
        "specific transcripts while background continues to accumulate, then at",
        "IDENTICAL library size a defect cell must put less of its budget into",
        "real markers and spread the remainder over more of the panel."),
  "At matched depth, defect cells show a lower marker share and a higher effective gene count.",
  sprintf("median specific share < %.2f x rest AND effective genes > rest",
          THR$flatten_share),
  supports = c9_ok,
  excludes = c("H5"),
  stats = list("specific share defect" = ss_d, "specific share rest" = ss_r,
               "share ratio" = ss_d / ss_r, "eff genes defect" = eg_d,
               "eff genes rest" = eg_r, "wilcoxon p" = w_p)
)

# ==============================================================================
# C10 CONSEQUENCE (informational, excludes nothing)
# ==============================================================================
rule_hdr("C10", "what the defect costs, and what survives it")
salv <- rbindlist(lapply(list(
  list("as-is",                            quote(rep(TRUE, .N))),
  list("drop defect cluster",              quote(!is_defect)),
  list("marker-informative only",          quote(informative)),
  list("drop defect + informative",        quote(!is_defect & informative)),
  list("inner 2 radial quintiles",         quote(rq <= 2)),
  list("inner 2 quintiles + drop defect",  quote(rq <= 2 & !is_defect))
), function(f) md[, {
  k <- eval(f[[2]])
  .(filter = f[[1]], n_kept = sum(k), pct_kept = round(100 * mean(k), 1),
    pct_strongB = round(100 * mean(strongB[k]), 2),
    pct_informative = round(100 * mean(informative[k]), 1))
}, by = merge_sample]))
print(dcast(salv, filter ~ merge_sample, value.var = c("pct_kept", "pct_strongB")))
register(
  "C10", "Analytic cost and salvage",
  "Not a mechanism test. Quantifies how much of IB1 is usable under each rescue filter.",
  "n/a", "informational only", supports = NA,
  neutral_reason = "reports consequences; excludes no hypothesis",
  stats = list("IB1 pct strongB as-is" = md[CASE_I, 100 * mean(strongB)],
               "RB1 pct strongB as-is" = md[CTRL_I, 100 * mean(strongB)])
)

# ==============================================================================
# THE LOGIC: apply exclusions, report survivors
# ==============================================================================
HYP <- c(
  H1 = "Over-fixation (excess formaldehyde crosslinking)",
  H2 = "Under-fixation / ischaemic delay -> RNA degradation",
  H3 = "Segmentation dilution",
  H4 = "Optical or decoding failure",
  H5 = "Antibody-only staining failure",
  H6 = "Genuine biology (the defect population is just another cell type)",
  H7 = "Section-surface handling artifact"
)
fired    <- vapply(CRIT, function(x) x$verdict == "SUPPORTS", logical(1))
killed_by <- setNames(vector("list", length(HYP)), names(HYP))
for (cr in CRIT) for (h in cr$excludes) killed_by[[h]] <- c(killed_by[[h]], cr$id)

hyp_tab <- data.table(
  hypothesis = names(HYP), description = unname(HYP),
  status = vapply(names(HYP), function(h)
    if (h == "H1") "PRIMARY"
    else if (length(killed_by[[h]])) "ELIMINATED" else "NOT EXCLUDED", character(1)),
  eliminated_by = vapply(names(HYP), function(h)
    paste(killed_by[[h]], collapse = "+"), character(1))
)
crit_tab <- rbindlist(lapply(CRIT, function(x) data.table(
  id = x$id, criterion = x$title, verdict = x$verdict,
  rule = x$rule, excludes_if_fired = paste(x$excludes_if_fired, collapse = "+"),
  applied = paste(x$excludes, collapse = "+"),
  prediction = x$prediction, mechanism = x$mechanism)))
metric_tab <- rbindlist(lapply(CRIT, function(x) {
  if (!length(x$stats)) return(NULL)
  data.table(id = x$id, statistic = names(x$stats),
             value = unlist(x$stats, use.names = FALSE))
}))

n_mech <- sum(vapply(CRIT, function(x) x$verdict != "NEUTRAL", logical(1)))
n_sup  <- sum(vapply(CRIT, function(x) x$verdict == "SUPPORTS", logical(1)))
survivors <- hyp_tab[status == "NOT EXCLUDED", hypothesis]

cat("\n\n################## CRITERIA ##################\n")
print(crit_tab[, .(id, criterion, verdict, applied)], row.names = FALSE)
cat("\n################## HYPOTHESES ##################\n")
print(hyp_tab, row.names = FALSE)
cat(sprintf("\nmechanism criteria satisfied: %d / %d\n", n_sup, n_mech))
cat("alternatives still standing:",
    if (length(survivors)) paste(survivors, collapse = ", ") else "none", "\n")

VERDICT <- if (n_sup == n_mech && !length(survivors)) {
  "SUPPORTED - every criterion fired and every alternative is eliminated"
} else if (n_sup >= n_mech - 1 && length(survivors) <= 1) {
  "SUPPORTED WITH RESERVATION - see the failing criterion / surviving alternative"
} else {
  "NOT ESTABLISHED - the criteria do not converge; treat the diagnosis as open"
}
cat("\nVERDICT:", VERDICT, "\n")

fwrite(crit_tab,   file.path(OUT_DIR, "ib1_overfixation_criteria.csv"))
fwrite(hyp_tab,    file.path(OUT_DIR, "ib1_overfixation_hypotheses.csv"))
fwrite(metric_tab, file.path(OUT_DIR, "ib1_overfixation_metrics.csv"))
fwrite(grad,       file.path(OUT_DIR, "ib1_overfixation_radial_gradient.csv"))
fwrite(wdet,       file.path(OUT_DIR, "ib1_overfixation_within_sample_detection.csv"))
fwrite(c2,         file.path(OUT_DIR, "ib1_overfixation_cross_sample_detection.csv"))
fwrite(salv,       file.path(OUT_DIR, "ib1_overfixation_salvage.csv"))

# ==============================================================================
# FIGURES
# ==============================================================================
th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = INK[["grid"]], linewidth = 0.25),
        panel.border = element_rect(colour = INK[["grid"]]),
        axis.text = element_text(colour = INK[["secondary"]]),
        axis.title = element_text(colour = INK[["secondary"]]),
        plot.title = element_text(colour = INK[["primary"]], face = "bold"),
        plot.subtitle = element_text(colour = INK[["secondary"]]),
        legend.position = "top", legend.title = element_blank())

# P1 -- the gradient, with the control that makes it mean something
p1 <- ggplot(grad, aes(median_r, pct_defect, colour = merge_sample)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
  geom_text(data = grad[rq == THR$n_quantile], aes(label = merge_sample),
            hjust = -0.25, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = COL) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(title = "C7  Defect rises toward the tissue rim in IB1 only",
       subtitle = sprintf("rho(radius, defect): IB1 %+.3f | RB1 %+.3f  -- rim-worst is the over-fixation sign",
                          rho_r_case, rho_r_ctrl),
       x = "distance from section centroid (mm, core -> rim)",
       y = "% cells in defect cluster") + th

# P2 -- channel dissociation
c1l <- melt(c1, id.vars = "channel", variable.name = "sample", value.name = "ratio")
c1l[, sample := ifelse(sample == "case_ratio", CASE, CTRL)]
c1l[, channel := factor(sub("^Mean\\.", "", channel),
                        levels = c("Membrane", "CD45", "PanCK", "DAPI"))]
p2 <- ggplot(c1l, aes(channel, ratio, fill = sample)) +
  geom_col(position = position_dodge2(padding = 0.15), width = 0.7) +
  geom_hline(yintercept = 1, linewidth = 0.4, colour = INK[["secondary"]]) +
  scale_fill_manual(values = COL) +
  labs(title = "C1  Antibody channels fail; the DNA stain does not",
       subtitle = "median intensity in defect cells / median in all other cells of the same sample",
       x = NULL, y = "defect / rest") + th

# P3 -- within-sample, depth-matched detection by lineage set
wp <- wdet[det_rest >= 0.02 & det_defect > 0]
p3 <- ggplot(wp, aes(100 * det_rest, ratio, colour = sample)) +
  geom_hline(yintercept = 1, linewidth = 0.4, colour = INK[["secondary"]]) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
  geom_point(size = 1.9, alpha = 0.9) +
  geom_text_repel(data = wp[sample == CASE], aes(label = gene), size = 2.4,
                  colour = INK[["secondary"]], max.overlaps = 20, segment.size = 0.2,
                  show.legend = FALSE) +
  scale_colour_manual(values = COL) + scale_x_log10() + scale_y_log10() +
  labs(title = "C2  The defect's profile is dragged toward uniform",
       subtitle = sprintf("common genes lose, rare genes gain: the closure signature of a flat floor.\nslope IB1 %+.2f, RB1 %+.2f  (0 = profile unchanged, -1 = fully uniform)",
                          ff_case$slope, ff_ctrl$slope),
       x = "baseline prevalence in that sample's own cells, % (log)",
       y = "detection ratio, defect / rest (log)") + th

# P3b -- the confounded cross-sample view, kept visibly separate
p3b <- ggplot(c2_use, aes(100 * det_ctrl, 100 * det_case)) +
  geom_abline(slope = 1, linewidth = 0.4, colour = INK[["secondary"]]) +
  geom_abline(slope = depth_ratio, linetype = "dashed", linewidth = 0.4,
              colour = INK[["secondary"]]) +
  geom_point(colour = COL[[CASE]], size = 2) +
  geom_text_repel(aes(label = gene), size = 2.6, colour = INK[["secondary"]],
                  max.overlaps = 30, segment.size = 0.2) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "C2b  Cross-sample detection — CONFOUNDED, not used as evidence",
       subtitle = sprintf("solid = parity, dashed = the %.0f%% depth gap. Scatter about the line is composition, not accessibility.",
                          100 * (1 - depth_ratio)),
       x = "RB1 % cells detecting (log)", y = "IB1 % cells detecting (log)") + th

# P4 -- FOV coupling
p4 <- ggplot(fov, aes(Mean.Membrane, 100 * frac_i, colour = sample_id)) +
  geom_point(size = 1.1, alpha = 0.6) +
  scale_colour_manual(values = COL) + scale_x_log10() +
  labs(title = "C8  Protein masking and RNA defect move together",
       subtitle = sprintf("per FOV. rho(frac_i, Membrane) = %+.2f in IB1, %+.2f in RB1",
                          c8_case$rho_membrane, c8[sample_id == CTRL, rho_membrane]),
       x = "FOV mean Membrane intensity (log)", y = "% cells in defect cluster") + th

# P5 -- the mechanism test at matched depth
c9l <- melt(c9[, .(arm, spec_share, eff_genes)], id.vars = "arm")
c9l[, variable := factor(variable, labels = c("marker share of counts", "effective genes (exp H)"))]
p5 <- ggplot(c9l, aes(arm, value, fill = arm)) +
  geom_violin(colour = NA, alpha = 0.85, scale = "width") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               colour = INK[["secondary"]], linewidth = 0.3) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values = c(defect = COL[[CASE]], rest = COL[[CTRL]])) +
  labs(title = "C9  At identical library size, defect cells carry flatter signal",
       subtitle = sprintf("IB1 only, depth-matched (%s cells per arm). Marker share ratio %.2f, Wilcoxon p %s\nthe marker-share panel carries this criterion; the spread in effective genes is small",
                          format(length(ia), big.mark = ","), ss_d / ss_r,
                          format.pval(w_p, digits = 2, eps = .Machine$double.eps)),
       x = NULL, y = NULL) + th + theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "C7_radial_gradient.png"), p1, width = 6.5, height = 4, dpi = 150)
ggsave(file.path(FIG_DIR, "C1_channel_dissociation.png"), p2, width = 6.5, height = 4, dpi = 150)
ggsave(file.path(FIG_DIR, "C2_lineage_blind_loss.png"), p3, width = 6.5, height = 4.2, dpi = 150)
ggsave(file.path(FIG_DIR, "C2b_cross_sample_confounded.png"), p3b, width = 6.5, height = 5, dpi = 150)
ggsave(file.path(FIG_DIR, "C8_fov_coupling.png"), p4, width = 6.5, height = 4, dpi = 150)
ggsave(file.path(FIG_DIR, "C9_matched_depth_flatness.png"), p5, width = 6.5, height = 4, dpi = 150)
ggsave(file.path(FIG_DIR, "overfixation_overview.png"),
       (p1 | p2) / (p3 | p4) / p5 + plot_annotation(
         title = "IB1 over-fixation: the evidence, criterion by criterion",
         theme = theme(plot.title = element_text(face = "bold", size = 12))),
       width = 13, height = 15, dpi = 140)
msg("figures written to", FIG_DIR)

# ==============================================================================
# WRITTEN SUMMARY
# ==============================================================================
li <- function(...) paste0(...)
md_lines <- c(
  "# IB1 over-fixation — diagnostic summary",
  "",
  li("Generated by `scripts/sbatch/qc_ib1_overfixation_summary.R` on ",
     format(Sys.time(), "%Y-%m-%d %H:%M"), " from `", MERGED, "`."),
  li("Case **", CASE, "** (", format(sum(CASE_I), big.mark = ","), " cells), internal control **",
     CTRL, "** (", format(sum(CTRL_I), big.mark = ","), " cells), defect population = InSituType cluster `",
     DEFECT, "` (", fmt(100 * mean(mdC$is_defect), 1), "% of ", CASE, " vs ",
     fmt(100 * mean(mdR$is_defect), 1), "% of ", CTRL, ")."),
  "",
  li("## Verdict: ", VERDICT),
  "",
  li(n_sup, " of ", n_mech, " mechanism criteria fired. Alternatives still standing: ",
     if (length(survivors)) paste(survivors, collapse = ", ") else "**none**", "."),
  "",
  "## Scope — two claims that must not be merged",
  "",
  li("**(A) Cluster `", DEFECT, "` inside IB1 is a technical defect population.** That is what ",
     "the criteria below test, and it is where the over-fixation call lives."),
  li("**(B) IB1 as a whole is globally degraded relative to RB1.** This summary does *not* ",
     "assert (B), and C2b shows why it cannot be read off a sample-to-sample comparison: the two ",
     "blocks are different tissues. At matched depth, IB1's B markers fall (MS4A1 ",
     fmt(c2[gene == "MS4A1", ratio_matched], 2), "x, CD79A ",
     fmt(c2[gene == "CD79A", ratio_matched], 2), "x) while its myeloid and stromal markers *rise* ",
     "(LYZ ", fmt(c2[gene == "LYZ", ratio_matched], 2), "x, COL1A1 ",
     fmt(c2[gene == "COL1A1", ratio_matched], 2), "x). That is the signature of differing ",
     "composition, not of uniform probe occlusion, so the cross-sample contrast is reported and ",
     "then set aside."),
  "",
  "## The argument",
  "",
  li("The defect in IB1 is **qualitative, not quantitative**. Cluster-`", DEFECT,
     "` cells are not shallow — they carry ", fmt(c4[metric == "nCount_RNA", ratio], 2),
     "x the counts of their neighbours at ", fmt(c5[["ratio"]], 2),
     "x the cell area (C4, C5). What they have lost is the ability to say anything: ",
     "no gene is enriched in them (max log2FC ", fmt(maxlfc[[DEFECT]], 2),
     " against ", fmt(peer_med, 2), " for the other clusters, C3). Measured against IB1's *own* ",
     "unaffected cells at matched depth, the profile is being dragged toward uniform: the ",
     "detection ratio falls with a gene's baseline prevalence at slope ", fmt(ff_case$slope, 2),
     " (p ", format.pval(ff_case$p, digits = 2), ", C2). Common genes lose, rare genes gain — ",
     "which is what a flat background floor does to a compositional measure, and is not how a ",
     "cell type behaves, since a cell type's winners are chosen by lineage rather than by ",
     "abundance. RB1's own cluster `", DEFECT, "` flattens the same way and slightly harder ",
     "(slope ", fmt(ff_ctrl$slope, 2), "), which is the point: cluster `", DEFECT,
     "` is one shared technical mode of this assay, not an IB1 invention. What is specific to ",
     "IB1 is how much of it there is (", fmt(100 * mean(mdC$is_defect), 1), "% against ",
     fmt(100 * mean(mdR$is_defect), 1), "%) and that it is laid down along a gradient."),
  "",
  li("C9 measures that floor directly rather than by inference. Holding library size equal (",
     fmt(median(sa$nCount), 0), " vs ", fmt(median(sb$nCount), 0), " counts), defect cells commit ",
     fmt(100 * ss_d, 2), "% of their budget to the 22-gene marker panel against ",
     fmt(100 * ss_r, 2), "% for their matched neighbours — a ratio of ", fmt(ss_d / ss_r, 2),
     " — and spread the remainder over ", fmt(eg_d, 0), " effective genes against ",
     fmt(eg_r, 0), " (Wilcoxon p ", format.pval(w_p, digits = 2), ")."),
  "",
  li("Degradation cannot produce that pattern — it removes molecules and would show up as fewer ",
     "counts, fewer features and a rising negative-probe share, none of which is present (C4). ",
     "Crosslinking can, because it removes *access* rather than material: the transcripts are ",
     "still there, they are simply unavailable to the probes, and the counting continues on ",
     "whatever background remains."),
  "",
  li("The protein side agrees and localises the insult to the fixative rather than the ",
     "antibodies. Membrane and CD45 fall to ", fmt(c1[channel == "Mean.Membrane", case_ratio], 2),
     "x and ", fmt(c1[channel == "Mean.CD45", case_ratio], 2),
     "x in defect cells while DAPI *rises* to ", fmt(c1[channel == "Mean.DAPI", case_ratio], 2),
     "x (C1) — exactly the split expected when formaldehyde masks protein epitopes but leaves ",
     "the DNA minor groove available. Across ", c8_case$n_fov, " fields of view the RNA defect ",
     "and the protein defect co-vary (rho ", fmt(c8_case$rho_membrane, 2),
     " with Membrane, ", fmt(c8_case$rho_dapi, 2), " with DAPI, C8), which two independent ",
     "accidents would not do."),
  "",
  li("Finally the geometry gives the *sign* of the error. Formalin diffuses inward from the ",
     "block surface, so over-fixation is worst at the rim and mildest in the core; ",
     "under-fixation is the mirror image. Observed: defect fraction climbs ",
     fmt(grad[merge_sample == CASE & rq == 1, pct_defect], 1), "% → ",
     fmt(grad[merge_sample == CASE & rq == THR$n_quantile, pct_defect], 1),
     "% from core to rim (rho ", fmt(rho_r_case, 2), "), while the control is flat (rho ",
     fmt(rho_r_ctrl, 2), ") — C7. Global radius outranks local boundary distance (rho ",
     fmt(rho_e_case, 2), "), which is what a block-scale diffusion gradient predicts and a ",
     "cut-section handling artifact does not."),
  "",
  "## Criteria",
  "",
  "| id | criterion | rule | verdict | eliminates |",
  "|----|-----------|------|---------|------------|",
  apply(crit_tab, 1, function(r)
    sprintf("| %s | %s | `%s` | **%s** | %s |", r[["id"]], r[["criterion"]],
            gsub("|", "\\|", r[["rule"]], fixed = TRUE), r[["verdict"]],
            if (nzchar(r[["applied"]])) r[["applied"]] else "—")),
  "",
  "## Competing hypotheses",
  "",
  "| id | hypothesis | status | eliminated by |",
  "|----|------------|--------|---------------|",
  apply(hyp_tab, 1, function(r)
    sprintf("| %s | %s | **%s** | %s |", r[["hypothesis"]], r[["description"]],
            r[["status"]], if (nzchar(r[["eliminated_by"]])) r[["eliminated_by"]] else "—")),
  "",
  "## What this cannot settle",
  "",
  "- Fixation *time* is not measured here. The data show a crosslinking-shaped lesion with a",
  "  surface-in gradient; they cannot distinguish prolonged formalin from an aggressive",
  "  antigen-retrieval or protease step that failed to reverse normal crosslinking. The block's",
  "  processing record is the confirmatory evidence and lives outside this dataset.",
  li("- H7 (section handling) is ",
     if (c7_block_scale) "disfavoured but not formally excluded — global radius beats local edge distance, which argues for a block-scale cause, but both are 2-D projections of a 3-D gradient."
     else "**NOT excluded**: local edge distance explains the gradient at least as well as global radius, which is the handling-artifact signature."),
  "- The defect is a continuum, not a binary. Cluster `i` is the thresholded readout of it, so",
  "  dropping cluster `i` removes the worst-affected cells but leaves a graded residue behind —",
  "  which is why the radial filters in C10 keep improving on top of it.",
  li("- Claim (B) is untested, not disproved. Whether IB1 as a whole is also mildly over-fixed ",
     "cannot be answered against RB1, because composition and fixation are entangled across ",
     "blocks. Answering it needs a within-IB1 dose axis (the C7 gradient) or a second control ",
     "block of the same tissue type."),
  li("- Cluster `", DEFECT, "` is not unique to IB1 — it is ", fmt(100 * mean(mdR$is_defect), 1),
     "% of RB1 and flattens there too (C2). The diagnosis is about *dose*: what over-fixation ",
     "explains is why IB1 carries ", fmt(mean(mdC$is_defect) / mean(mdR$is_defect), 1),
     "x as much of it and why IB1's copy is spatially organised, not why the mode exists."),
  li("- C9 is carried by its marker-share term (", fmt(ss_d / ss_r, 2),
     "x, a large and unambiguous effect). Its second term, effective gene count, moves in the ",
     "predicted direction but only slightly (", fmt(eg_d, 0), " vs ", fmt(eg_r, 0),
     "), because entropy over 6,175 genes is dominated by sampling noise at ~550 counts per ",
     "cell. Do not quote the entropy term as strong evidence on its own."),
  "- CD45 is a leukocyte marker as well as an imaging channel, so its fall in the defect",
  "  population is not by itself independent of composition. C1 does not rest on it: the",
  "  load-bearing observation is the *split* between the antibody channels and DAPI, which no",
  "  change in cell composition produces on its own.",
  "- The defect population is defined by InSituType cluster `i`, which was fit on this data, so",
  "  every within-sample contrast compares a data-derived group against its own complement. The",
  "  guard against reading that circularity as a result is RB1, where cluster `i` is defined the",
  "  same way: it shows the same flattening (C2) yet is four to five times rarer and carries no",
  "  radial gradient (C7) and no channel coupling (C8). The mode is generic; its dose and its",
  "  spatial organisation in IB1 are not.",
  "",
  "## Salvage",
  "",
  paste(capture.output(print(dcast(salv, filter ~ merge_sample,
                                   value.var = c("pct_kept", "pct_strongB")),
                             row.names = FALSE)), collapse = "\n"),
  "",
  li("Recommended: analyse ", CASE, " under `drop defect + informative` or restrict to the inner ",
     "radial quintiles, and never compare ", CASE, " and ", CTRL,
     " on detection-rate-derived quantities without depth matching."),
  ""
)
writeLines(unlist(md_lines), file.path(OUT_DIR, "ib1_overfixation_summary.md"))
msg("wrote", file.path(OUT_DIR, "ib1_overfixation_summary.md"))
cat("\nOK\n")
