## ============================================================================
## Step 36 -- Is a step-35 composite ONE carved cell or SEVERAL bright cells?
## One-vs-two molecule clouds, per touching member pair.
##
## WHY. Step 35 showed chemokine-bright polygons are contiguous far beyond chance, and the
## molecule overlay on its panels showed two kinds of composite: one CCL17/CCL22 cloud
## continuous across the shared polygon edge (a cell the segmentation carved up), or one
## cloud per member with a gap between (a group of separately bright cells). This step
## classifies every pair without the eye.
##
## STATISTIC (primary candidate): the PROJECTED dip. The two members' molecules are
## projected onto the axis from a's centroid to b's; a 1-D Gaussian KDE (sigma SIG_UM) of
## the projected positions; the polygon boundary sits at b_pos on that axis (from the two
## polygons' support functions); dip = density at b_pos / min(peak on a's side, peak on
## b's side). One cloud cut by the boundary is unimodal, so the boundary is at or below a
## flank and dip ~ 1; two clouds give a saddle at the boundary and dip < 1. Projection
## uses every molecule, which is what fixed two earlier versions: a line through the
## centroids sampled at sigma 2 um merged clouds 10 um apart, and at sigma 1 um it missed
## any cloud not sitting on the line.
## Secondary candidate: the labelled dBIC, BIC(one Gaussian over all) - BIC(one Gaussian
## per member), on the same points. (An unsupervised mixtools EM was tried and failed on
## a third of the two-cell reference -- spurious optima with a near-empty component.)
##
## REFERENCES, count-matched to each real pair (bins on counts):
##   one-cell  an ISOLATED bright polygon (single-member component at the same threshold)
##             matched to the pair's total count, cut by a random line at a random offset
##             within +/- 0.3 r of its centroid: its own molecules, labelled by side.
##   two-cell  two isolated polygons matched to a and b separately; b's molecules are
##             translated so that b's POLYGON is tangent to a's along a random direction
##             (support functions), the boundary at the tangency. Each cloud keeps its real
##             offset from its own centroid.
## The cut per rule is the equal-error point of the two references; the PRIMARY rule is
## whichever separates the references better (decided on the references, not on the
## real pairs), and both classes are written out.
##
## IB1 only: the molecule table (step 29) exists for IB1 only, and RB1 showed contiguity
## only at >= 10.
##
## Usage: Rscript scripts/36_script_molecule_dip.R [smoke]
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(parallel)
})
RNGkind("L'Ecuyer-CMRG"); set.seed(36)
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
SMOKE <- "smoke" %in% commandArgs(trailingOnly = TRUE)
IN  <- "output/35_giant_cells"
OUT <- if (SMOKE) "output/36_molecule_dip_smoke" else "output/36_molecule_dip"
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)

UM_PER_PX   <- 0.120280945
FOV_PX      <- 4256L
SIG_UM      <- 1.0
GRID_UM     <- 0.25
MIN_CONTACT <- 3L
MIN_MOL     <- 5L                        # fewer molecules on a side than this: no statistic
THR         <- c(10L, 20L, 50L)
CNT_BRK     <- c(10, 20, 35, 50, 100, 200, Inf)
MAX_PAIRS   <- if (SMOKE) 300L else Inf
CLASS_COLS  <- c(real = "grey20", `one-cell ref` = "#0b3d7a", `two-cell ref` = "#eb6834")

theme_viz <- function(base = 11) theme_minimal(base) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", plot.title.position = "plot")

## ---- inputs ---------------------------------------------------------------------------
fl   <- fread(file.path(IN, "fov_list.csv"))[samp == "IB1"]
asg  <- fread(file.path(IN, "polygon_assignments.csv.gz"))[samp == "IB1"]
cells <- fread(file.path(IN, "cells_IB1.csv.gz"),
               select = c("fov", "cell_ID", "Area", "CenterX_local_px", "CenterY_local_px", "chemo"))
cells <- merge(cells, fl[, .(fov, x0, y0)], by = "fov")
cells[, `:=`(cx = (x0 + (FOV_PX - 1L) - CenterX_local_px) * UM_PER_PX,      # canvas um
             cy = (y0 + CenterY_local_px) * UM_PER_PX,
             r_um = sqrt(Area / pi) * UM_PER_PX)]                          # vendor Area is px^2
setkey(cells, fov, cell_ID)
tx <- fread("output/29_bright_cells/ib1_chemokine_tx.csv.gz")[cell_ID > 0]
tx <- tx[, .(fov = as.integer(fov), cell_ID = as.integer(cell_ID), x = as.numeric(X_px) * UM_PER_PX, y = as.numeric(Y_px) * UM_PER_PX)]
setkey(tx, fov, cell_ID)
cat(sprintf("assigned molecules %d; cells %d; positive polygons %d\n", nrow(tx), nrow(cells), nrow(asg)))
mols <- function(f, id) as.matrix(tx[.(f, id), .(x, y), nomatch = NULL])

edges <- rbindlist(lapply(unique(asg$fov), function(f) {
  p <- file.path(IN, "adj", sprintf("IB1_%d_edges.csv", f))
  if (!file.exists(p)) return(NULL)
  fread(p)[contact_px >= MIN_CONTACT]
}))
cat(sprintf("touching pairs (any cells) %d\n", nrow(edges)))
cbin <- function(x) cut(x, CNT_BRK, right = FALSE, include.lowest = TRUE)

## ---- real pairs, both partitions, each threshold ----------------------------------------
pairs <- rbindlist(lapply(THR, function(t) {
  A <- asg[thr == t]
  rbindlist(lapply(c("comp", "anch"), function(part) {
    A[, cid := get(part)]
    multi <- A[, .N, by = cid][N >= 2, cid]
    M <- A[cid %in% multi]
    E <- merge(edges, M[, .(fov, a = cell_ID, cid_a = cid, chemo_a = chemo, nuc_a = has_nuc)], by = c("fov", "a"))
    E <- merge(E, M[, .(fov, b = cell_ID, cid_b = cid, chemo_b = chemo, nuc_b = has_nuc)], by = c("fov", "b"))
    E <- E[cid_a == cid_b]
    E[, `:=`(thr = t, partition = part, cid = cid_a)]
    E[, has_seed := cid %in% A[seed == TRUE, cid]]
    E[, .(thr, partition, cid, has_seed, fov, a, b, contact_px, chemo_a, chemo_b, nuc_a, nuc_b)]
  }), use.names = TRUE)
}))
cat("member pairs:\n"); print(pairs[, .N, by = .(thr, partition)])
uniq <- unique(pairs[, .(thr, fov, a, b, chemo_a, chemo_b)])
if (is.finite(MAX_PAIRS)) uniq <- uniq[sample.int(nrow(uniq), min(MAX_PAIRS, nrow(uniq)))]
uniq[, `:=`(bin_total = cbin(chemo_a + chemo_b), bin_a = cbin(chemo_a), bin_b = cbin(chemo_b), bin_min = cbin(pmin(chemo_a, chemo_b)))]
cat(sprintf("unique pairs to score: %d\n", nrow(uniq)))

## isolated = single-member COMPONENT at the pair's own threshold (nothing positive touches it)
iso <- rbindlist(lapply(THR, function(t) {
  A <- asg[thr == t]; A[, n_in := .N, by = comp]
  A[n_in == 1, .(thr = t, fov, cell_ID, chemo, bin = cbin(chemo))]
}))
cat("isolated positive polygons by count bin:\n"); print(iso[, .N, by = .(thr, bin)][order(thr, bin)])

## ---- polygon vertices (relative to the vendor centroid, um) for every cell we touch -------
need <- unique(rbind(uniq[, .(fov, cell_ID = a)], uniq[, .(fov, cell_ID = b)], iso[, .(fov, cell_ID)]))
cat(sprintf("polygons needed: %d cells in %d FOVs\n", nrow(need), uniqueN(need$fov)))
VERT <- new.env()
for (f in unique(need$fov)) {
  P <- fread(file.path(IN, "polygons", sprintf("IB1_%d.csv", f)), header = FALSE, col.names = c("fov", "cell_ID", "xl", "yl"))
  P <- P[cell_ID %in% need[fov == f, cell_ID]]
  o <- fl[fov == f]
  P[, `:=`(x = (o$x0 + (FOV_PX - 1L) - xl) * UM_PER_PX, y = (o$y0 + yl) * UM_PER_PX)]
  P <- merge(P, cells[.(f), .(cell_ID, cx, cy)], by = "cell_ID")
  for (id in unique(P$cell_ID)) assign(sprintf("%d_%d", f, id), as.matrix(P[cell_ID == id, .(x - cx, y - cy)]), envir = VERT)
}
support <- function(f, id, u) {           # polygon extent from the centroid in direction u
  V <- get0(sprintf("%d_%d", f, id), envir = VERT, inherits = FALSE)
  if (is.null(V)) return(NA_real_)
  max(V %*% u)
}

## ---- the statistics --------------------------------------------------------------------
kde1 <- function(t, grid) {
  if (length(t) == 0) return(rep(0, length(grid)))
  colSums(exp(-outer(t, grid, "-")^2 / (2 * SIG_UM^2))) / (sqrt(2 * pi) * SIG_UM)
}
ll_one <- function(X) {
  n <- nrow(X); if (n < 3) return(NA_real_)
  S <- cov(X) * (n - 1) / n; d <- det(S); if (!is.finite(d) || d <= 0) return(NA_real_)
  -n / 2 * (2 * log(2 * pi) + log(d) + 2)
}
proj_stats <- function(Xa, Xb, ca, u, b_pos) {
  ## Xa, Xb: molecules (um) of the two sides; ca: axis origin; u: unit axis; b_pos: boundary on the axis
  na <- nrow(Xa); nb <- nrow(Xb)
  if (na < MIN_MOL || nb < MIN_MOL) return(data.table(dip = NA_real_, dbic = NA_real_, n_a = na, n_b = nb, peak_a = NA_real_, peak_b = NA_real_, f_boundary = NA_real_))
  X <- rbind(Xa, Xb) + matrix(runif(2 * (na + nb), -0.06, 0.06), na + nb)   # molecules sit on integer px
  t <- as.vector((X - matrix(ca, na + nb, 2, byrow = TRUE)) %*% u)
  grid <- seq(min(t) - 3, max(t) + 3, by = GRID_UM)
  f <- kde1(t, grid)
  fa <- max(f[grid <= b_pos]); fb <- max(f[grid >= b_pos]); f0 <- kde1(t, b_pos)
  dip <- if (min(fa, fb) > 0) f0 / min(fa, fb) else NA_real_
  n <- na + nb
  l1 <- ll_one(X); l2 <- ll_one(X[1:na, , drop = FALSE]) + ll_one(X[(na + 1):n, , drop = FALSE]) + na * log(na / n) + nb * log(nb / n)
  dbic <- if (is.na(l1) || is.na(l2)) NA_real_ else (-2 * l1 + 5 * log(n)) - (-2 * l2 + 11 * log(n))
  data.table(dip = dip, dbic = dbic, n_a = na, n_b = nb, peak_a = fa, peak_b = fb, f_boundary = f0)
}
pair_stats <- function(f, a, b) {
  ca <- cells[.(f, a), nomatch = NULL]; cb <- cells[.(f, b), nomatch = NULL]
  if (nrow(ca) != 1 || nrow(cb) != 1) return(NULL)
  d <- c(cb$cx - ca$cx, cb$cy - ca$cy); sep <- sqrt(sum(d^2)); u <- d / sep
  ha <- support(f, a, u); hb <- support(f, b, -u)                    # a's reach toward b, b's reach toward a
  b_pos <- if (is.na(ha) || is.na(hb)) sep / 2 else (ha + (sep - hb)) / 2
  cbind(proj_stats(mols(f, a), mols(f, b), c(ca$cx, ca$cy), u, b_pos), sep_um = sep, b_pos_um = b_pos)
}
one_cell_stats <- function(f, o) {
  co <- cells[.(f, o), nomatch = NULL]
  th <- runif(1, 0, 2 * pi); u <- c(cos(th), sin(th))
  b_pos <- runif(1, -0.3, 0.3) * co$r_um
  X <- mols(f, o); if (nrow(X) == 0) return(NULL)
  t <- as.vector((X - matrix(c(co$cx, co$cy), nrow(X), 2, byrow = TRUE)) %*% u)
  cbind(proj_stats(X[t < b_pos, , drop = FALSE], X[t >= b_pos, , drop = FALSE], c(co$cx, co$cy), u, b_pos), sep_um = NA_real_, b_pos_um = b_pos)
}
two_cell_stats <- function(fa, a, fb, b) {
  ca <- cells[.(fa, a), nomatch = NULL]; cb <- cells[.(fb, b), nomatch = NULL]
  th <- runif(1, 0, 2 * pi); u <- c(cos(th), sin(th))
  ha <- support(fa, a, u); hb <- support(fb, b, -u)
  if (is.na(ha) || is.na(hb)) return(NULL)
  sep <- ha + hb                                                       # polygons tangent along u
  nb_c <- c(ca$cx, ca$cy) + sep * u
  Xb <- mols(fb, b); Xb <- Xb - matrix(c(cb$cx, cb$cy), nrow(Xb), 2, byrow = TRUE) + matrix(nb_c, nrow(Xb), 2, byrow = TRUE)
  cbind(proj_stats(mols(fa, a), Xb, c(ca$cx, ca$cy), u, ha), sep_um = sep, b_pos_um = ha)
}

## ---- score the real pairs --------------------------------------------------------------
t0 <- Sys.time()
dips <- rbindlist(mclapply(seq_len(nrow(uniq)), mc.cores = NC, mc.set.seed = TRUE, FUN = function(i) {
  r <- uniq[i]; d <- pair_stats(r$fov, r$a, r$b)
  if (is.null(d)) return(NULL)
  cbind(r, d)
}))
cat(sprintf("scored %d pairs in %.1f min\n", nrow(dips), as.numeric(Sys.time() - t0, units = "mins")))

## ---- references, count-matched per real pair --------------------------------------------
draw_iso <- function(t, bn) {
  lv <- levels(iso$bin); k <- match(bn, lv)
  for (dk in 0:length(lv)) for (kk in unique(c(k - dk, k + dk))) if (kk >= 1 && kk <= length(lv)) {
    pool <- iso[thr == t & bin == lv[kk]]
    if (nrow(pool)) return(pool[sample.int(nrow(pool), 1)])
  }
  stop("no isolated polygon at threshold ", t)
}
t0 <- Sys.time()
refs <- rbindlist(mclapply(seq_len(nrow(dips)), mc.cores = NC, mc.set.seed = TRUE, FUN = function(i) {
  r <- dips[i]
  o <- draw_iso(r$thr, r$bin_total); d1 <- one_cell_stats(o$fov, o$cell_ID)
  a <- draw_iso(r$thr, r$bin_a); b <- draw_iso(r$thr, r$bin_b); d2 <- two_cell_stats(a$fov, a$cell_ID, b$fov, b$cell_ID)
  rbind(if (!is.null(d1)) data.table(thr = r$thr, pair_i = i, kind = "one-cell ref", chemo_a = o$chemo, chemo_b = NA_integer_, d1),
        if (!is.null(d2)) data.table(thr = r$thr, pair_i = i, kind = "two-cell ref", chemo_a = a$chemo, chemo_b = b$chemo, d2))
}))
cat(sprintf("references in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
refs[, bin_min := dips$bin_min[pair_i]]
fwrite(refs, file.path(OUT, "reference_stats.csv"))

## ---- the cuts: equal-error point between the two references, per rule ------------------
## scores are oriented so that HIGHER = two clouds: -dip, and dbic
equal_error <- function(r1, r2) {
  fin <- c(r1, r2)[is.finite(c(r1, r2))]
  grid <- sort(unique(c(quantile(fin, seq(0.005, 0.995, 0.005)), seq(min(fin), max(fin), length.out = 200))))
  err <- sapply(grid, function(c) c(mean(r1 >= c), mean(r2 < c)))
  k <- which.min(abs(err[1, ] - err[2, ]))
  c(cut = unname(grid[k]), one_wrong = unname(err[1, k]), two_wrong = unname(err[2, k]))
}
refs[, `:=`(s_dip = -dip, s_dbic = dbic)]; dips[, `:=`(s_dip = -dip, s_dbic = dbic)]
cuts <- rbindlist(lapply(c("dip", "dbic"), function(rule) {
  sc <- paste0("s_", rule)
  r1 <- refs[kind == "one-cell ref" & !is.na(get(sc)), get(sc)]; r2 <- refs[kind == "two-cell ref" & !is.na(get(sc)), get(sc)]
  e <- equal_error(r1, r2)
  data.table(rule = rule, cut = if (rule == "dip") -e[["cut"]] else e[["cut"]],
             one_cell_misclassified = e[["one_wrong"]], two_cell_misclassified = e[["two_wrong"]],
             mean_error = (e[["one_wrong"]] + e[["two_wrong"]]) / 2, n_one_ref = length(r1), n_two_ref = length(r2),
             median_one_ref = median(refs[kind == "one-cell ref", get(rule)], na.rm = TRUE),
             median_two_ref = median(refs[kind == "two-cell ref", get(rule)], na.rm = TRUE))
}))
PRIMARY <- cuts[which.min(mean_error), rule]
cuts[, primary := rule == PRIMARY]
cat("cuts (primary = the rule that separates the references better):\n"); print(cuts)
fwrite(cuts, file.path(OUT, "cut.csv"))
cut_dip <- cuts[rule == "dip", cut]; cut_dbic <- cuts[rule == "dbic", cut]
s_primary <- paste0("s_", PRIMARY); cut_s <- if (PRIMARY == "dip") -cut_dip else cut_dbic

cut_by_bin <- refs[!is.na(get(s_primary)), {
  a1 <- get(s_primary)[kind == "one-cell ref"]; a2 <- get(s_primary)[kind == "two-cell ref"]
  e <- equal_error(a1, a2)
  .(rule = PRIMARY, cut = if (PRIMARY == "dip") -e[["cut"]] else e[["cut"]], n = .N,
    one_wrong = e[["one_wrong"]], two_wrong = e[["two_wrong"]],
    one_wrong_at_pooled_cut = mean(a1 >= cut_s), two_wrong_at_pooled_cut = mean(a2 < cut_s))
}, by = .(thr, bin_min)][order(thr, bin_min)]
cat("cut by count bin:\n"); print(cut_by_bin)
fwrite(cut_by_bin, file.path(OUT, "cut_by_count_bin.csv"))

## ---- classify pairs and composites -----------------------------------------------------
dips[, class := fifelse(is.na(get(s_primary)), NA_character_, fifelse(get(s_primary) >= cut_s, "two_clouds", "one_cloud"))]
dips[, class_dip  := fifelse(is.na(dip), NA_character_, fifelse(dip < cut_dip, "two_clouds", "one_cloud"))]
dips[, class_dbic := fifelse(is.na(dbic), NA_character_, fifelse(dbic >= cut_dbic, "two_clouds", "one_cloud"))]
P <- merge(pairs, dips[, .(thr, fov, a, b, dip, dbic, n_a, n_b, peak_a, peak_b, f_boundary, sep_um, b_pos_um,
                            bin_min, bin_total, class, class_dip, class_dbic)], by = c("thr", "fov", "a", "b"))
fwrite(P, file.path(OUT, "pair_class.csv"))

comp <- P[!is.na(class), .(n_pairs = .N, n_one = sum(class == "one_cloud"), frac_one = mean(class == "one_cloud"),
                           min_dip = min(dip), has_seed = has_seed[1]),
          by = .(thr, partition, cid, fov)]
comp[, class := fifelse(n_one == n_pairs, "all_one_cloud", fifelse(n_one == 0, "all_separate", "mixed"))]
fwrite(comp, file.path(OUT, "composite_class.csv"))

summ <- rbind(
  P[!is.na(class), .(unit = "pair", n = .N, frac_one_cloud = mean(class == "one_cloud"), median_dip = median(dip, na.rm = TRUE),
                     frac_one_cloud_dip = mean(class_dip == "one_cloud", na.rm = TRUE),
                     frac_one_cloud_dbic = mean(class_dbic == "one_cloud", na.rm = TRUE)),
    by = .(thr, partition, seed_composites = has_seed)],
  comp[, .(unit = "composite", n = .N, frac_one_cloud = mean(class == "all_one_cloud"), median_dip = median(min_dip),
           frac_mixed = mean(class == "mixed"), frac_all_separate = mean(class == "all_separate")),
       by = .(thr, partition, seed_composites = has_seed)], fill = TRUE)[order(unit, thr, partition, -seed_composites)]
cat("summary:\n"); print(summ)
fwrite(summ, file.path(OUT, "summary.csv"))

by_bin <- P[!is.na(class) & partition == "comp", .(n = .N, frac_one_cloud = mean(class == "one_cloud"), median_dip = median(dip)),
            by = .(thr, bin_min)][order(thr, bin_min)]
cat("pairs by count bin (components):\n"); print(by_bin)
fwrite(by_bin, file.path(OUT, "pairs_by_count_bin.csv"))

## nucleus-bearing vs nucleus-free pairs: is the dense-nucleus member a separate cell more often?
by_nuc <- P[!is.na(class) & partition == "anch", .(n = .N, frac_one_cloud = mean(class == "one_cloud")),
            by = .(thr, pair_type = fifelse(nuc_a & nuc_b, "both nucleated", fifelse(nuc_a | nuc_b, "nucleus + fragment", "two fragments")))][order(thr, pair_type)]
cat("pairs by nucleus status (anchored):\n"); print(by_nuc)
fwrite(by_nuc, file.path(OUT, "pairs_by_nucleus.csv"))

## the 24 composites on the step-35 panel sheet (thr 20, anchored, largest by area): the eye check
sc <- fread(file.path(IN, "seed_composites.csv.gz"))[samp == "IB1" & thr == 20 & partition == "anchored" & members >= 2]
top <- sc[order(-area_um2)][1:min(24, .N), cid]
chk <- P[thr == 20 & partition == "anch" & cid %in% top][order(match(cid, top))]
chk_out <- chk[, .(cid, fov, a, b, chemo_a, chemo_b, nuc_a, nuc_b, dip = round(dip, 2), dbic = round(dbic, 1), class, class_dip, class_dbic)]
cat("panel-sheet composites (eye check):\n"); print(chk_out, nrows = 200)
fwrite(chk_out, file.path(OUT, "panel_sheet_pairs.csv"))

## ---- figures ---------------------------------------------------------------------------
hist_d <- rbind(dips[, .(thr, kind = "real", dip, dbic)], refs[, .(thr, kind, dip, dbic)])
hist_d[, kind := factor(kind, names(CLASS_COLS))]
lab_thr <- labeller(thr = function(x) paste0("positive >= ", x, " counts"))
p1 <- ggplot(hist_d[!is.na(dip)], aes(dip, colour = kind, fill = kind)) +
  geom_density(alpha = 0.15, bw = 0.03) +
  geom_vline(xintercept = cut_dip, linetype = 2) +
  facet_wrap(~thr, labeller = lab_thr) +
  scale_colour_manual(values = CLASS_COLS, name = NULL) + scale_fill_manual(values = CLASS_COLS, name = NULL) +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  labs(title = sprintf("%s: projected molecule-density dip at the polygon boundary (sigma %g um)", if (PRIMARY == "dip") "PRIMARY" else "Secondary", SIG_UM),
       subtitle = sprintf("dip = density at the boundary / smaller side's peak. Dashed: equal-error cut %.2f (one-cell ref called two-clouds %.0f%%, two-cell ref called one-cloud %.0f%%)",
                          cut_dip, 100 * cuts[rule == "dip", one_cell_misclassified], 100 * cuts[rule == "dip", two_cell_misclassified]),
       x = "dip (1 = no saddle at the boundary)", y = "density") +
  theme_viz()
p1b <- ggplot(hist_d[!is.na(dbic)], aes(pmin(pmax(dbic, -40), 120), colour = kind, fill = kind)) +
  geom_density(alpha = 0.15, bw = 3) +
  geom_vline(xintercept = cut_dbic, linetype = 2) +
  facet_wrap(~thr, labeller = lab_thr) +
  scale_colour_manual(values = CLASS_COLS, name = NULL) + scale_fill_manual(values = CLASS_COLS, name = NULL) +
  labs(title = sprintf("%s: labelled dBIC, one Gaussian over both sides vs one per side", if (PRIMARY == "dbic") "PRIMARY" else "Secondary"),
       subtitle = sprintf("Dashed: equal-error cut %.1f (one-cell ref called two-clouds %.0f%%, two-cell ref called one-cloud %.0f%%). Clipped to [-40, 120].",
                          cut_dbic, 100 * cuts[rule == "dbic", one_cell_misclassified], 100 * cuts[rule == "dbic", two_cell_misclassified]),
       x = "dBIC (positive favours two clouds)", y = "density") +
  theme_viz()
ggsave(file.path(OUT, "figures/fig36_statistic_distributions.png"), p1 / p1b, width = 11, height = 8, dpi = 150, bg = "white")

bb <- rbind(dips[!is.na(dip), .(thr, kind = "real", bin_min, dip)], refs[!is.na(dip), .(thr, kind, bin_min, dip)])
bb[, kind := factor(kind, names(CLASS_COLS))]
p2 <- ggplot(bb, aes(bin_min, dip, colour = kind)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = cut_dip, linetype = 2) +
  facet_wrap(~thr, labeller = lab_thr) +
  scale_colour_manual(values = CLASS_COLS, name = NULL) +
  labs(title = "Dip by the dimmer member's count", subtitle = "references are count-matched to each real pair",
       x = "CCL17+CCL22 counts of the dimmer member (bin)", y = "dip") +
  theme_viz() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT, "figures/fig36_dip_by_count.png"), p2, width = 11, height = 4.2, dpi = 150, bg = "white")

cs <- comp[, .N, by = .(thr, partition, has_seed, class)][, share := N / sum(N), by = .(thr, partition, has_seed)]
cs[, class := factor(class, c("all_one_cloud", "mixed", "all_separate"))]
cs[, seed := fifelse(has_seed, "composites with a seed (>= 50)", "other composites")]
p3 <- ggplot(cs, aes(factor(thr), share, fill = class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 3, colour = "white") +
  facet_grid(seed ~ partition, labeller = labeller(partition = c(comp = "connected components", anch = "nucleus-anchored"))) +
  scale_fill_manual(values = c(all_one_cloud = "#0b3d7a", mixed = "grey55", all_separate = "#eb6834"), name = NULL) +
  labs(title = sprintf("Composite class from its member pairs (%s rule)", PRIMARY),
       subtitle = sprintf("one_cloud = the members' molecules are continuous across the boundary; reference error ~%.0f%% each way", 100 * cuts[primary == TRUE, mean_error]),
       x = "positivity threshold (counts)", y = "share of multi-member composites") +
  theme_viz()
ggsave(file.path(OUT, "figures/fig36_composite_class.png"), p3, width = 8, height = 6, dpi = 150, bg = "white")

cat("done ->", OUT, "\n")
