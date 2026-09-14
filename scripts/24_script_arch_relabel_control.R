#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 24b -- THE CONTROL step 24 left outstanding: labels or power?
##
## WHAT STEP 24 FOUND, AND WHY IT COULD NOT CONCLUDE. Re-running step 16's co-occurrence
## permutation under the step-19 labels dropped RB1's fraction-of-pairs-at-chance from
## **0.755** (step-16 argmax labels) to **0.503**, i.e. step 16's headline "RB1 is spatially
## homogeneous, essentially every pair at chance" looked much weaker. But step 19 names a
## different NUMBER of cells than the step-16 argmax does, and the permutation's power scales
## with how many cells carry a label. So two explanations were confounded:
##      (i)  the better labels reveal real structure in RB1, or
##      (ii) more labelled cells simply buy statistical power.
## This script separates them.
##
## THE CONTROL. For every region and every shared cell type, take
##      n_type = min(n labelled by step 16, n labelled by step 19)
## and down-sample BOTH label sets to it. Then n_a and n_b are identical between the two runs
## for every pair by construction, so any surviving difference in the enrichment or in the
## fraction at chance is the labels and not the sample size. Note the matching runs in BOTH
## directions: step 19 does not uniformly name more cells -- it names ~1.24x as many B cells
## but only ~0.45x as many T cells, because step 16's T count is the sum of its four T
## subsets. A one-directional "down-sample step 19" control would have been wrong.
##
## THE GEOMETRY IS HELD FIXED, WHICH STEP 24 DID NOT DO. Step 24's `run_one` subset to the
## labelled cells BEFORE building the neighbour graph, so the two label sets saw different
## neighbour graphs and different local-density deciles -- a second, smaller confound sitting
## underneath the first. Here the graph and the deciles are built ONCE per region from ALL
## cells present, and a cell that a given label source does not name simply carries NA and
## contributes nothing. Down-sampling therefore removes LABELS, never cells: it mimics "the
## caller named fewer cells", which is the power question, rather than "the tissue had fewer
## cells", which is a different one.
##
## Two conditions are run so the effect of each change is visible rather than asserted:
##   A  common geometry, FULL labels      -- isolates the geometry fix against step 24
##   B  common geometry, MATCHED n        -- the control; B is the answer
## Condition B is repeated over several seeds, because down-sampling adds variance of its own.
##
## Usage: sbatch --mem=64G --time=4:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/24_script_arch_relabel_control.R [tier] [smoke]
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(RANN); library(Matrix); library(data.table)})
setDTthreads(1); options(width = 200)

AG    <- commandArgs(TRUE)
TIER  <- if (length(AG) >= 1) AG[1] else "matched"
SMOKE <- length(AG) >= 2 && AG[2] == "smoke"
REG_F <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT16  <- "output/16_celltypes/cell_celltypes.csv.gz"
CT19  <- "output/19_celltypes/all_cells_celltypes.csv.gz"
OUT   <- if (SMOKE) file.path("output/24_arch_relabel", TIER, "smoke_control") else
                    file.path("output/24_arch_relabel", TIER)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

## constants verbatim from 16_script_cooccurrence.R / 24_script_arch_relabel.R
R_UM <- 20; K_MAX <- 40L; N_PERM <- 200L; MIN_A <- 100L; N_DEC <- 10L
N_SEED <- 3L; SEED0 <- 240L
NICHE_LABEL <- "Hodgkin chemokine niche"
if (SMOKE) { N_PERM <- 20L; N_SEED <- 2L }

MAP19 <- c("B cell"="B", "T cell"="T", "Myeloid"="myeloid", "Stromal/FDC"="fibro",
           "Endothelial"="endo", "NK cell"="NK", "Plasma cell"="PC", "Plasmablast"="PC")
MAP16 <- c("B"="B", "T_unresolved"="T", "Tfh"="T", "Treg"="T", "CD8_T"="T", "CD4_conv"="T",
           "Monocyte"="myeloid", "Macrophage"="myeloid", "cDC"="myeloid",
           "myeloid_unresolved"="myeloid", "fibro"="fibro", "endo"="endo",
           "NK"="NK", "PC"="PC")
SHARED <- c("B","T","myeloid","fibro","endo","NK","PC")

cat("========================================================================\n")
if (SMOKE) cat("*** SMOKE MODE -> ", OUT, " ***\n", sep = "")
cat("step 24b -- labels or power? tier =", TIER, "\n")
cat("========================================================================\n")

md <- fread(REG_F, select = c("cell_id","samp","fov_id","x_slide_mm","y_slide_mm",
                              "bin_region","is_defect","tier_matched"))
md <- md[tier_matched == TRUE & is_defect == FALSE & bin_region != "unbinned"]
md <- merge(md, fread(CT16, select = c("cell_id","celltype")), by = "cell_id", all.x = TRUE)
md <- merge(md, fread(CT19, select = c("cell_id","celltype_cluster")), by = "cell_id", all.x = TRUE)
md[, lab16 := MAP16[celltype]]
md[, lab19 := fifelse(celltype_cluster == NICHE_LABEL, NA_character_, MAP19[celltype_cluster])]
md[!lab16 %in% SHARED, lab16 := NA_character_]
md[!lab19 %in% SHARED, lab19 := NA_character_]

cat("\n=== named cells per type, and the matched target n = min of the two\n")
cnt <- md[, .(n16 = sum(!is.na(lab16)), n19 = sum(!is.na(lab19))), by = .(samp, bin_region)]
per <- rbindlist(lapply(SHARED, function(t)
  md[, .(type = t, n16 = sum(lab16 == t, na.rm = TRUE), n19 = sum(lab19 == t, na.rm = TRUE)),
     by = .(samp, bin_region)]))
per[, n_match := pmin(n16, n19)][, ratio_19_over_16 := round(n19/pmax(n16,1), 2)]
print(per[bin_region %in% c("background","FL_follicular","interfollicular")][order(samp, bin_region, type)])
sc("control_matched_n.csv", per)

## ---- the statistic, geometry built ONCE per region from ALL cells ----------------------
prep <- function(z) {
  xy <- as.matrix(z[, .(x_slide_mm, y_slide_mm)]) * 1000
  k  <- min(K_MAX, nrow(z) - 1L)
  nn <- RANN::nn2(xy, xy, k = k + 1L)
  idx <- nn$nn.idx[, -1, drop = FALSE]; dst <- nn$nn.dists[, -1, drop = FALSE]
  dcol <- min(20L, ncol(dst)); dens <- dst[, dcol]
  strat <- cut(dens, quantile(dens, seq(0, 1, length.out = N_DEC + 1L)),
               include.lowest = TRUE, labels = FALSE)
  list(idx = idx, keep = dst <= R_UM, strat = strat, n = nrow(z))
}
run_stat <- function(G, lab0) {
  types <- intersect(SHARED, unique(lab0[!is.na(lab0)]))
  if (length(types) < 2L) return(NULL)
  ti <- setNames(seq_along(types), types)
  stat <- function(lab) {
    li <- ti[lab]
    m <- matrix(0, length(types), length(types), dimnames = list(types, types))
    for (j in seq_len(ncol(G$idx))) {
      ok <- G$keep[, j]; if (!any(ok)) next
      a <- li[ok]; b <- li[G$idx[ok, j]]
      good <- !is.na(a) & !is.na(b); if (!any(good)) next
      m <- m + as.matrix(sparseMatrix(i = a[good], j = b[good], x = 1,
             dims = c(length(types), length(types))))
    }
    sweep(m, 1, pmax(as.numeric(table(factor(lab, levels = types))), 1), "/")
  }
  obs <- stat(lab0)
  perms <- lapply(seq_len(N_PERM), function(p) {
    lab <- lab0
    for (s in unique(G$strat[!is.na(G$strat)])) {
      ii <- which(G$strat == s); lab[ii] <- sample(lab0[ii])
    }
    stat(lab)
  })
  nu <- Reduce(`+`, perms) / N_PERM
  ge <- Reduce(`+`, lapply(perms, function(m) m >= obs))
  na <- as.numeric(table(factor(lab0, levels = types)))
  data.table(a = rep(types, times = length(types)), b = rep(types, each = length(types)),
             n_a = rep(na, times = length(types)), n_b = rep(na, each = length(types)),
             obs = as.vector(obs), null = as.vector(nu),
             p_emp = (1 + as.vector(ge)) / (1 + N_PERM))[, enrich := obs/pmax(null,1e-9)][]
}
downsample <- function(lab, target) {   # remove LABELS, never cells
  out <- lab
  for (t in names(target)) {
    ii <- which(lab == t); need <- target[[t]]
    if (length(ii) > need) out[sample(ii, length(ii) - need)] <- NA_character_
  }
  out
}

REGIONS <- md[, .N, by = .(samp, bin_region)][N >= 3*MIN_A][order(samp, -N)]
if (SMOKE) REGIONS <- REGIONS[1:min(3, .N)]
cat("\nregions:", nrow(REGIONS), "\n")

res <- list()
for (i in seq_len(nrow(REGIONS))) {
  SM <- REGIONS$samp[i]; BR <- REGIONS$bin_region[i]
  z <- md[samp == SM & bin_region == BR]
  G <- prep(z)
  tgt <- per[samp == SM & bin_region == BR]
  target <- setNames(tgt$n_match, tgt$type)
  cat(sprintf("\n--- %s / %s: %s cells | named 16/19 = %s/%s\n", SM, BR,
      format(nrow(z), big.mark=","), format(sum(!is.na(z$lab16)), big.mark=","),
      format(sum(!is.na(z$lab19)), big.mark=",")))
  ## condition A -- full labels, common geometry
  for (L in c("lab16","lab19")) {
    o <- run_stat(G, z[[L]])
    if (!is.null(o)) res[[length(res)+1]] <-
      cbind(samp=SM, bin_region=BR, condition="A_full", labels=L, seed=NA_integer_, o)
  }
  ## condition B -- matched n, several seeds
  for (s in seq_len(N_SEED)) {
    set.seed(SEED0 + s)
    for (L in c("lab16","lab19")) {
      o <- run_stat(G, downsample(z[[L]], target))
      if (!is.null(o)) res[[length(res)+1]] <-
        cbind(samp=SM, bin_region=BR, condition="B_matched", labels=L, seed=s, o)
    }
  }
}

R <- rbindlist(res, fill = TRUE)
R <- R[n_a >= MIN_A & n_b >= MIN_A]
R[, q := p.adjust(p_emp, "BH"), by = .(samp, bin_region, condition, labels, seed)]
sc("control_cooccurrence.csv.gz", R)

## ---- the answer, and it is judged THRESHOLD-FREE first --------------------------------
## `frac_at_chance` is what step 24 reported, so it is reproduced for comparability, but it is
## a poor primary statistic: q comes from a permutation p floored at 1/(N_PERM+1), and BH only
## lets pairs through when MANY of them pile up at that floor together. So the count of
## "significant" pairs is partly a property of the floor. The primary comparison here is the
## EFFECT SIZE itself -- how far the enrichments sit from 1 -- which has no threshold in it,
## plus a paired per-pair contrast between the two label sources at matched n.
R[, l2 := log2(pmax(enrich, 1e-9))]
eff <- R[, .(pairs = .N, med_abs_l2 = round(median(abs(l2)), 4),
             p90_abs_l2 = round(quantile(abs(l2), .9), 4)),
         by = .(samp, condition, labels, seed)]
effs <- eff[, .(seeds = .N, med_abs_l2 = round(mean(med_abs_l2), 4),
                lo = round(min(med_abs_l2), 4), hi = round(max(med_abs_l2), 4)),
            by = .(samp, condition, labels)]
cat("\n\n=== PRIMARY, threshold-free: how far do the enrichments sit from chance?\n")
print(effs[order(samp, condition, labels)])
sc("control_effect_size.csv", effs)

## paired, per (region x pair x seed), matched n only
W <- dcast(R[condition == "B_matched"], samp + bin_region + a + b + seed ~ labels,
           value.var = "l2")
W <- W[!is.na(lab16) & !is.na(lab19)]
pair_cmp <- W[, {
  d <- abs(lab19) - abs(lab16)
  st <- suppressWarnings(binom.test(sum(d > 0), sum(d != 0))$p.value)
  .(n = .N, med_diff_abs_l2 = round(median(d), 4),
    frac_lab19_stronger = round(mean(d > 0), 3), p_sign = signif(st, 3),
    rho = round(suppressWarnings(cor(lab16, lab19, method = "spearman")), 3))
}, by = samp]
cat("\n=== PAIRED, matched n: is the SAME pair further from chance under step-19 labels?\n")
print(pair_cmp)
sc("control_paired.csv", pair_cmp)

## ---- the thresholded version, reproduced for comparability with step 24 ----------------
fc <- R[, .(pairs = .N, med_enrich = round(median(enrich), 3),
            frac_at_chance = round(mean(q >= 0.05), 3)),
        by = .(samp, condition, labels, seed)]
summ <- fc[, .(seeds = .N, frac_at_chance = round(mean(frac_at_chance), 3),
               lo = round(min(frac_at_chance), 3), hi = round(max(frac_at_chance), 3),
               pairs = round(mean(pairs))), by = .(samp, condition, labels)]
cat("\n\n=== fraction of pairs at chance, by condition and label source\n")
print(summ[order(samp, condition, labels)])
sc("control_chance_fraction.csv", summ)

cat("\n=== THE COMPARISON THAT ANSWERS IT (condition B, matched n)\n")
B <- summ[condition == "B_matched"]
for (SM in unique(B$samp)) {
  a <- B[samp == SM & labels == "lab16", frac_at_chance]
  b <- B[samp == SM & labels == "lab19", frac_at_chance]
  if (!length(a) || !length(b)) next
  cat(sprintf("  %s: lab16 %.3f vs lab19 %.3f  (difference %+.3f)\n", SM, a, b, b - a))
}
cat("\n  Step 24 reported RB1 0.755 (lab16) vs 0.503 (lab19) on UNMATCHED n and each\n")
cat("  label source's own geometry. Judge it on the threshold-free table and the paired\n")
cat("  test above, not on this one: if the matched difference is near zero there, step 24's\n")
cat("  gap was power; if it persists, the labels carry it.\n")

writeLines(format(Sys.time()), file.path(OUT, "arch_control.done"))
cat("\nOK\n")
