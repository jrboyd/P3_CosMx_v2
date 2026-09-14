# =============================================================================
# 02 — spatial FOV quality control (formalized from dev_FOV_investigate.R)
#
# Emits per-cell QC flags for the two spatial artifacts measured in this project.
# Runs BEFORE step 10, because that is where the artifact does its damage
# (clustering -> PCA -> HierAnnot). Reads the flat metadata directly, so it needs
# neither an 8 GB readRDS nor the BPCells store.
#
# ARTIFACT 1 — within-FOV position bias (BOTH samples). Counts per cell fall to
#   -52 % at the high-x_FOV_px edge; asymmetric high_x >> low_y > high_y > low_x;
#   extends ~2000 px. Not segmentation (area flat, holds in unsplit cells), not
#   decoding (negprobe/falsecode fractions flat), not acquisition order (all Order
#   correlations vanish after removing a smooth spatial trend).
#   -> handled by FOUR automatically fitted per-edge cuts (section 4) plus a
#      whole-FOV over-correlation rule for FOVs anomalous beyond the systematic
#      effect (section 6).
#
# ARTIFACT 2 — staining/section failure, IB1 ONLY. InSituType cluster `i` is 59 %
#   of IB1 (vs 8 % RB1), has no enriched genes (max log2FC +0.76 vs +1.65..+6.75
#   elsewhere) and a panel-average profile. Per-FOV frac_i correlates -0.58 with
#   Mean.Membrane and -0.56 with Mean.CD45 but POSITIVELY with DAPI: nuclei
#   stained, membrane/CD45 antibodies did not.
#   -> handled by a 3-channel logistic rule fit on IB1 (section 5).
#
# WHY THE STAINING RULE IS IB1-ONLY AND WHY THAT IS ASSERTED, NOT ASSUMED.
# RB1 has no staining problem and its channel scale differs enough (sample-wide
# Mean.Membrane 662 vs 1688; Membrane/DAPI 0.20 vs 0.69) that pooling it would
# corrupt the fit. So the scope is deliberate. But in the prototype the scope arose
# from `is.na -> "ok"`, which ALSO silently converts a failed join into a pass.
# Here: out-of-scope samples get the literal "not_evaluated", in-scope samples must
# receive a prediction or the script stops, and RB1 gets the full diagnostic panel
# as the positive control that justifies the exemption.
#
# COORDINATES. Verified identical: the flat file's CenterX/Y_local_px equal the
# Seurat meta.data x_FOV_px / y_FOV_px exactly (20,592/20,592 cells checked). So
# "high x" means the same thing here as in the measured profiles.
#
# Usage:  Rscript 09b_spatial_qc.R          (processes both samples in one run)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)      # loaded up front: the prototype used `+` before this
})

set.seed(1)                                    # covers every RNG path below
options(width = 200)

# ------------------------------------------------------------------- CONFIG ----
OUT_DIR   <- "output/02_filter_QC"
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES   <- c("IB1", "RB1")
NETFLAT   <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles"
# The repo copy of RB1_metadata_file.csv.gz is a TRUNCATED gzip (fails `gzip -t`).
# Always read from /netfiles.

FOV_PX    <- 4256L          # FOV side, px. Pitch 0.512 mm, no overlap.
PX_UM     <- 0.120280945

EDGE_TOL      <- 0.10       # accept a bin once within 10 % of the interior plateau
EDGE_BIN      <- 25L        # px per bin for the per-edge profile fit
EDGE_PLATEAU  <- 1000L      # bins beyond this px define the interior plateau
EDGE_CUT_MODE <- "shared_max"   # "shared_max" (identical geometry both samples,
#  conservative) or "per_sample"

STAINING_SCOPE <- "IB1"                       # explicit, not emergent
CHANNELS       <- c("Mean.DAPI", "Mean.Membrane", "Mean.CD45")  # deliberate 3
I_FRAC_CUT     <- 0.30                        # frac_i defining a "junk" FOV, no effect on filter, only labelling
CV_K           <- 5L

CORR_CUT  <- 0.20   # |Spearman rho| of nCount vs local px, cell level
DENS_CUT  <- 0.60   # |Spearman rho| of per-bin cell count vs bin index
CORR_NMIN <- 200L   # min post-cut cells in a FOV to judge it

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
save_fig <- function(p, name, w = 9, h = 7, dpi = 130) {
  ggsave(file.path(FIG_DIR, name), p, width = w, height = h, dpi = dpi, limitsize = FALSE)
  invisible(file.path(FIG_DIR, name))
}
theme_set(theme_bw(base_size = 9))

# ------------------------------------------------------------------ 1. READ ----
READ_COLS <- c("cell", "fov", "CenterX_local_px", "CenterY_local_px",
               "CenterX_global_px", "CenterY_global_px",
               "nCount_RNA", "nFeature_RNA", "nCount_negprobes",
               "Area.um2", "qcCellsPassed",
               "Mean.DAPI", "Max.DAPI", "Mean.Membrane", "Max.Membrane",
               "Mean.CD45", "Max.CD45", "Mean.PanCK", "Max.PanCK",
               "Mean.G", "Max.G")

read_sample <- function(S) {
  gz  <- file.path(NETFLAT, S, paste0(S, "_metadata_file.csv.gz"))
  stopifnot(file.exists(gz))
  hdr <- strsplit(readLines(gzfile(gz), 1L), ",", fixed = TRUE)[[1]]
  ist <- grep("InSituType.*clusters$", hdr, value = TRUE)
  stopifnot(length(ist) == 1L)
  keep <- c(intersect(READ_COLS, hdr), ist)
  d <- fread(cmd = paste("zcat", shQuote(gz)), select = keep, showProgress = FALSE)
  setnames(d, ist, "inSitu")
  d[, sample_id := S]
  # canonical local coords; verified equal to Seurat's x_FOV_px / y_FOV_px
  setnames(d, c("CenterX_local_px", "CenterY_local_px"), c("x_FOV_px", "y_FOV_px"))
  d[, fov := as.integer(fov)]
  d[]
}

msg("reading flat metadata for", paste(SAMPLES, collapse = ", "))
md <- rbindlist(lapply(SAMPLES, read_sample), use.names = TRUE, fill = TRUE)
msg("cells:", nrow(md))
print(md[, .N, by = sample_id])
stopifnot(!anyNA(md$x_FOV_px), !anyNA(md$y_FOV_px), !anyNA(md$fov))
# fov is NOT unique across samples -> always key on (sample_id, fov)
stopifnot(md[, uniqueN(fov), by = sample_id][, all(V1 > 1)])

# distance from each of the four edges
md[, `:=`(d_low_x  = x_FOV_px,
          d_high_x = (FOV_PX - 1L) - x_FOV_px,
          d_low_y  = y_FOV_px,
          d_high_y = (FOV_PX - 1L) - y_FOV_px)]
EDGES <- c("low_x", "high_x", "low_y", "high_y")

# --------------------------------------------- 2. per-FOV channel summary ------
CHAN_COLS <- grep("^(Mean|Max)\\.", names(md), value = TRUE)
fov_chan <- md[, c(lapply(.SD, mean, na.rm = TRUE), .(n_cells = .N)),
               by = .(sample_id, fov), .SDcols = CHAN_COLS]
chan_long <- melt(fov_chan, id.vars = c("sample_id", "fov", "n_cells"),
                  variable.name = "stat_chan", value.name = "signal")
chan_long[, c("stat", "channel") := tstrsplit(as.character(stat_chan), ".", fixed = TRUE)]
# z and relative signal WITHIN sample (the prototype pooled samples, which with a
# 2.6x channel-scale difference encodes "which sample" rather than "which FOV")
chan_long[, `:=`(
  z_signal   = (signal - mean(signal, na.rm = TRUE)) / sd(signal, na.rm = TRUE),
  rel_signal = signal / max(signal, na.rm = TRUE)
), by = .(sample_id, channel, stat)]

# ------------------------------------------------- 3. per-FOV cluster-i burden --
# The InSituType letters ARE comparable across samples -- verified two ways:
#   (a) 39_ H1: cross-sample correlation of per-cluster marker profiles gives 9/10
#       same-letter best matches, mean diagonal r = 0.866. Under independent
#       per-slide fits letter assignment would be arbitrary.
#   (b) RB1's cluster `d` is only 0.9 % of RB1 yet is the argmax for CCL17, CCL22
#       AND BATF3 in RB1 (4.6x / 5.8x / 2.2x its own sample mean) -- the same
#       defining trio as IB1's abundant `d` (6.8x / 6.6x / 4.1x). A shared profile
#       set, not coincidence.
# So cluster `i` is ONE shared "low-information" profile, and IB1 having 59 % of its
# cells there against RB1's 8 % is a statement about IB1's data quality measured on
# a yardstick common to both samples. It is still used here only as a QC indicator,
# never as a cell-type label.
fov_i <- md[, .(n_cells = .N, n_i = sum(inSitu == "i", na.rm = TRUE)),
            by = .(sample_id, fov)][, frac_i := n_i / n_cells][]
fov_i[, z_frac_i := (frac_i - mean(frac_i)) / sd(frac_i), by = sample_id]
msg("frac_i quantiles by sample:")
print(fov_i[, as.list(round(quantile(frac_i, c(0, .1, .25, .5, .75, .9, 1)), 3)), by = sample_id])

# ================================= 4. AUTOMATIC PER-EDGE CUT FIT (4 cuts) ======
# The prototype's min_px / max_px were read by eye off plots built from an UNSEEDED
# slice_sample(n = 30), so the thresholds driving the filter were not reproducible.
# Here each of the four edges is fit automatically from the pooled depth profile.
#
# Fitting x and y separately is not enough: the asymmetry is WITHIN an axis
# (high_x -35.7 % at 175 px vs low_x -6.4 %), so all four edges get their own cut.
fit_edge_cut <- function(dist, ncount, tol = EDGE_TOL) {
  b <- data.table(d = (dist %/% EDGE_BIN) * EDGE_BIN + EDGE_BIN %/% 2L, y = ncount)
  prof <- b[, .(y = mean(y), n = .N), by = d][order(d)][n >= 50]
  plateau <- prof[d > EDGE_PLATEAU, median(y)]
  if (!is.finite(plateau) || plateau <= 0) return(list(cut = NA_integer_, plateau = NA_real_,
                                                       method = "no_plateau", prof = prof))
  target <- (1 - tol) * plateau
  # saturating fit: y = plateau * (1 - a * exp(-d / tau))
  fit <- try(nls(y ~ plateau * (1 - a * exp(-d / tau)), data = prof,
                 start = list(a = 0.5, tau = 200),
                 control = nls.control(warnOnly = TRUE, maxiter = 200)), silent = TRUE)
  cut <- NA_integer_; method <- "crossing"
  if (!inherits(fit, "try-error")) {
    co <- coef(fit)
    if (all(is.finite(co)) && co[["a"]] > 0 && co[["tau"]] > 0) {
      # solve plateau*(1 - a*exp(-d/tau)) = target  ->  d = -tau*log((1-target/plateau)/a)
      arg <- (1 - target / plateau) / co[["a"]]
      if (arg > 0) { cut <- as.integer(ceiling(-co[["tau"]] * log(arg))); method <- "nls" }
    }
  }
  if (is.na(cut) || cut < 0 || cut > FOV_PX / 2) {
    # fallback: first bin from which every subsequent bin is above target
    ok <- prof$y >= target
    run <- rev(cumprod(rev(as.integer(ok))))       # 1 while all-later-are-ok
    cut <- if (any(run == 1L)) as.integer(prof$d[which(run == 1L)[1]]) else NA_integer_
    method <- "crossing"
  }
  list(cut = cut, plateau = plateau, method = method, prof = prof)
}

msg("== fitting per-edge cuts (tol =", EDGE_TOL, ") ==")
edge_fits <- rbindlist(lapply(SAMPLES, function(S) rbindlist(lapply(EDGES, function(E) {
  d <- md[sample_id == S]
  f <- fit_edge_cut(d[[paste0("d_", E)]], d$nCount_RNA)
  data.table(sample_id = S, edge = E, cut_px = f$cut, plateau = f$plateau, method = f$method)
}))))
print(edge_fits)
edge_prof <- rbindlist(lapply(SAMPLES, function(S) rbindlist(lapply(EDGES, function(E) {
  d <- md[sample_id == S]
  p <- fit_edge_cut(d[[paste0("d_", E)]], d$nCount_RNA)$prof
  p[, `:=`(sample_id = S, edge = E)][]
}))))

EDGE_CUTS <- if (EDGE_CUT_MODE == "shared_max") {
  z <- edge_fits[, .(cut_px = max(cut_px, na.rm = TRUE)), by = edge]
  msg("EDGE_CUT_MODE = shared_max -> identical geometry in both samples")
  setNames(z$cut_px, z$edge)
} else {
  msg("EDGE_CUT_MODE = per_sample")
  NULL
}
if (!is.null(EDGE_CUTS)) print(EDGE_CUTS)

cut_for <- function(S, E) if (is.null(EDGE_CUTS)) edge_fits[sample_id == S & edge == E]$cut_px else EDGE_CUTS[[E]]

# per-cell edge verdict (single label, explicit precedence -- the prototype's
# sequential ifelse silently kept only the last matching label)
md[, edge_fail := NA_character_]
for (S in SAMPLES) for (E in EDGES) {
  cut <- cut_for(S, E); if (!is.finite(cut)) next
  md[sample_id == S & is.na(edge_fail) & get(paste0("d_", E)) < cut,
     edge_fail := paste0(E, "_lt_", cut, "px")]
}
msg("cells failing the edge cuts:")
print(md[, .(n = .N, failed = sum(!is.na(edge_fail)),
             frac_failed = round(mean(!is.na(edge_fail)), 4)), by = sample_id])

# ======================== 5. STAINING RULE — IB1 ONLY, EXPLICIT SCOPE ==========
msg("== staining rule; scope =", paste(STAINING_SCOPE, collapse = ", "), "==")
chan_wide <- dcast(chan_long[stat == "Mean" & channel %in% sub("^Mean\\.", "", CHANNELS)],
                   sample_id + fov ~ paste0("z_", channel), value.var = "z_signal")
fov_tab <- merge(merge(fov_i, chan_wide, by = c("sample_id", "fov")),
                 fov_chan[, c("sample_id", "fov", CHAN_COLS), with = FALSE],
                 by = c("sample_id", "fov"))
fov_tab[, junk := frac_i > I_FRAC_CUT]

# Predictors are the THREE channel z-scores, named explicitly. A `grep("^z_")`
# here also catches z_frac_i -- the z-score of the outcome itself -- which is
# target leakage: `junk` is a deterministic function of frac_i, so the fit
# separates perfectly (SE ~3e4), probabilities saturate to {0,1}, the Youden cut
# collapses to 0 and every FOV is flagged. Select by name, and assert it.
zcols <- paste0("z_", sub("^Mean\\.", "", CHANNELS))
stopifnot(all(zcols %in% names(fov_tab)))
stopifnot(!any(grepl("frac_i|junk", zcols)))          # no outcome-derived predictor

train <- fov_tab[sample_id %in% STAINING_SCOPE]
stopifnot(nrow(train) > 50, train[, uniqueN(junk)] == 2L)
form <- as.formula(paste("junk ~", paste(zcols, collapse = " + ")))
gfit <- glm(form, data = train, family = binomial())
msg("logistic coefficients (z-scored channels, ", paste(STAINING_SCOPE, collapse = "/"), "):")
print(round(summary(gfit)$coefficients, 4))
# separation guard: runaway coefficients mean a predictor encodes the outcome
if (max(abs(coef(gfit)[-1])) > 25)
  stop("logistic fit looks separated (|coef| > 25); check for target leakage in: ",
       paste(zcols, collapse = ", "))

# Youden-optimal cut on the fitted probability, over a proper interior grid
train[, p_junk := predict(gfit, newdata = train, type = "response")]
grid <- sort(unique(c(seq(0.01, 0.99, by = 0.005), round(train$p_junk, 3))))
grid <- grid[grid > 0 & grid < 1]
you <- vapply(grid, function(t) {
  tp <- train[, sum(p_junk >= t & junk)]; fn <- train[, sum(p_junk < t & junk)]
  fp <- train[, sum(p_junk >= t & !junk)]; tn <- train[, sum(p_junk < t & !junk)]
  tp / max(tp + fn, 1) + tn / max(tn + fp, 1) - 1
}, 0)
P_CUT <- grid[which.max(you)]
# AUC, so we can see whether the three channels discriminate at all
auc <- with(train, {
  r <- rank(p_junk); n1 <- sum(junk); n0 <- sum(!junk)
  (sum(r[junk]) - n1 * (n1 + 1) / 2) / (n1 * n0)
})
msg(sprintf("Youden cut = %.4f  J = %.3f  AUC = %.3f  (%d junk / %d clean FOVs)",
            P_CUT, max(you), auc, train[, sum(junk)], train[, sum(!junk)]))

# k-fold CV of the same procedure (the prototype reported an in-sample matrix)
folds <- sample(rep_len(seq_len(CV_K), nrow(train)))
cv <- rbindlist(lapply(seq_len(CV_K), function(k) {
  tr <- train[folds != k]; te <- train[folds == k]
  if (tr[, uniqueN(junk)] < 2L) return(NULL)
  g <- glm(form, data = tr, family = binomial())
  te[, p := predict(g, newdata = te, type = "response")]
  data.table(fold = k, n = nrow(te),
             err = te[, mean((p >= P_CUT) != junk)],
             sens = te[junk == TRUE, mean(p >= P_CUT)],
             spec = te[junk == FALSE, mean(p < P_CUT)])
}))
msg(sprintf("%d-fold CV: error %.3f  sens %.3f  spec %.3f",
            CV_K, mean(cv$err, na.rm = TRUE), mean(cv$sens, na.rm = TRUE),
            mean(cv$spec, na.rm = TRUE)))
print(cv)

# apply: in-scope samples get a real verdict, out-of-scope get "not_evaluated".
# Level names are assigned BY NAME (the prototype used levels(pred)[1] <- ...,
# which breaks silently if the labels change).
fov_tab[, p_junk := NA_real_]
fov_tab[sample_id %in% STAINING_SCOPE,
        p_junk := predict(gfit, newdata = .SD, type = "response"), .SDcols = zcols]
fov_tab[, qc_spatial_staining := fifelse(
  !(sample_id %in% STAINING_SCOPE), "not_evaluated",
  fifelse(p_junk >= P_CUT, "dim_membrane", "ok"))]
# in-scope FOVs MUST have a verdict; a missing prediction is an error, not a pass
stopifnot(!any(is.na(fov_tab$qc_spatial_staining)))
stopifnot(fov_tab[sample_id %in% STAINING_SCOPE, !any(is.na(p_junk))])
stopifnot(all(fov_tab$qc_spatial_staining %in% c("ok", "dim_membrane", "not_evaluated")))
msg("staining verdict by sample (FOVs):")
print(dcast(fov_tab[, .N, by = .(sample_id, qc_spatial_staining)],
            sample_id ~ qc_spatial_staining, value.var = "N", fill = 0))

# ============ 6. WHOLE-FOV OVER-CORRELATION RULE (post-cut cells only) ========
# Complete the wiring the prototype intended: `toss` becomes its own
# qc_spatial_fovpx level. Computed on POST-CUT cells so it flags FOVs anomalous
# BEYOND the systematic edge effect section 4 already removes.
#
# Decision uses effect size only. The prototype's p-value cut (1e-25) is
# n-dependent -- FOV cell counts span ~1.1k-4.3k, so it preferentially discards
# cell-dense FOVs for the same underlying gradient. The analytic p is reported but
# not used.
msg("== over-correlation rule on post-cut cells ==")
post <- md[is.na(edge_fail)]
corr_fov <- rbindlist(lapply(c("x_FOV_px", "y_FOV_px"), function(ax) {
  post[, {
    if (.N < CORR_NMIN) .(axis = ax, n = .N, rho = NA_real_, p = NA_real_, dens_rho = NA_real_)
    else {
      ct <- suppressWarnings(cor.test(get(ax), nCount_RNA, method = "spearman"))
      # density gradient: per-bin cell COUNT vs bin index. An acquisition gradient
      # drops counts per cell at constant density; an off-tissue edge drops density.
      bin <- (get(ax) %/% 532L)
      db  <- data.table(bin = bin)[, .N, by = bin][order(bin)]
      dr  <- if (nrow(db) >= 4L) suppressWarnings(cor(db$bin, db$N, method = "spearman")) else NA_real_
      .(axis = ax, n = .N, rho = unname(ct$estimate), p = ct$p.value, dens_rho = dr)
    }
  }, by = .(sample_id, fov)]
}))
corr_wide <- dcast(corr_fov, sample_id + fov ~ axis, value.var = c("rho", "p", "dens_rho", "n"))
corr_wide[, max_abs_rho := pmax(abs(rho_x_FOV_px), abs(rho_y_FOV_px), na.rm = TRUE)]
corr_wide[, max_abs_dens := pmax(abs(dens_rho_x_FOV_px), abs(dens_rho_y_FOV_px), na.rm = TRUE)]
corr_wide[, fov_flag := fifelse(
  is.na(max_abs_rho), "insufficient_cells",
  fifelse(max_abs_rho >= CORR_CUT & max_abs_dens >= DENS_CUT, "FOV_partial_tissue",
          fifelse(max_abs_rho >= CORR_CUT, "FOV_discarded_overcorrelation", "ok")))]
msg("FOV-level over-correlation verdict:")
print(dcast(corr_wide[, .N, by = .(sample_id, fov_flag)],
            sample_id ~ fov_flag, value.var = "N", fill = 0))

# ------------------------------------------------- 7. ASSEMBLE THE VERDICT -----
fov_verdict <- merge(fov_tab[, .(sample_id, fov, frac_i, z_frac_i, p_junk,
                                 qc_spatial_staining)],
                     corr_wide[, .(sample_id, fov, max_abs_rho, max_abs_dens, fov_flag)],
                     by = c("sample_id", "fov"), all = TRUE)
md <- merge(md, fov_verdict, by = c("sample_id", "fov"), all.x = TRUE)
stopifnot(!any(is.na(md$qc_spatial_staining)))    # no silent default-to-pass

# precedence: whole-FOV rejection > per-cell edge > ok
md[, qc_spatial_fovpx := fifelse(
  fov_flag %in% c("FOV_discarded_overcorrelation", "FOV_partial_tissue"), fov_flag,
  fifelse(!is.na(edge_fail), edge_fail, "ok"))]
md[, qc_spatial_i_fraction := fifelse(is.na(frac_i), NA_character_,
                                      fifelse(frac_i > I_FRAC_CUT, "i_too_high", "ok"))]
# qc_spatial_passing = the SPATIAL verdict alone (edge cuts + whole-FOV rules +
# staining). qc_spatial_i_fraction stays reported-not-enforced, to avoid selecting
# on the outcome the staining rule is calibrated against.
md[, qc_spatial_passing := qc_spatial_fovpx == "ok" &
     qc_spatial_staining %in% c("ok", "not_evaluated")]
# qc_passing_all additionally requires NanoString's own per-cell verdict. This is
# the column downstream steps should consume; the two components stay separate so
# it is always visible which rule dropped a cell. Cost of the extra AND is small
# because the spatial filter already removes most qcCellsPassed == FALSE cells:
# 1,016 further cells in IB1 (0.25 %) and 3,345 in RB1 (0.24 %).
md[, qcCellsPassed_lgl := qcCellsPassed %in% c(TRUE, "TRUE")]
md[, qc_passing_all := qc_spatial_passing & qcCellsPassed_lgl]

msg("=========== FINAL VERDICT ===========")
print(md[, .(cells = .N,
             spatial_pass = sum(qc_spatial_passing),
             passing_all  = sum(qc_passing_all),
             frac_all     = round(mean(qc_passing_all), 4)), by = sample_id])
msg("qc_spatial_fovpx breakdown:")
print(dcast(md[, .N, by = .(sample_id, qc_spatial_fovpx)],
            qc_spatial_fovpx ~ sample_id, value.var = "N", fill = 0))

# ------------------------------------------------------------ 8. OUTPUTS -------
cells_out <- md[, .(cell, cell_id_prefixed = paste0(sample_id, "_", cell),
                    sample_id, fov,
                    qc_spatial_fovpx, qc_spatial_staining, qc_spatial_i_fraction,
                    qcCellsPassed, qc_spatial_passing, qc_passing_all)]
fwrite(cells_out, file.path(OUT_DIR, "qc_spatial_cells.csv.gz"))
fwrite(merge(fov_verdict, fov_chan, by = c("sample_id", "fov"), all.x = TRUE),
       file.path(OUT_DIR, "qc_spatial_fov.csv"))
fwrite(edge_fits,  file.path(OUT_DIR, "qc_spatial_edge_fits.csv"))
fwrite(edge_prof,  file.path(OUT_DIR, "qc_spatial_edge_profiles.csv"))
fwrite(chan_long,  file.path(OUT_DIR, "qc_spatial_fov_channels.csv"))

attr_rows <- rbindlist(list(
  md[, .(rule = "total",                  n = .N),                                   by = sample_id],
  md[, .(rule = "fail_edge_cut",          n = sum(!is.na(edge_fail))),               by = sample_id],
  md[, .(rule = "fail_overcorrelation",   n = sum(qc_spatial_fovpx == "FOV_discarded_overcorrelation")), by = sample_id],
  md[, .(rule = "fail_partial_tissue",    n = sum(qc_spatial_fovpx == "FOV_partial_tissue")),            by = sample_id],
  md[, .(rule = "fail_staining",          n = sum(qc_spatial_staining == "dim_membrane")),               by = sample_id],
  md[, .(rule = "fail_qcCellsPassed",     n = sum(!qcCellsPassed_lgl)),               by = sample_id],
  md[, .(rule = "passing_spatial_only",   n = sum(qc_spatial_passing)),              by = sample_id],
  md[, .(rule = "passing_all",            n = sum(qc_passing_all)),                  by = sample_id]))
fwrite(dcast(attr_rows, rule ~ sample_id, value.var = "n", fill = 0),
       file.path(OUT_DIR, "qc_spatial_attrition.csv"))
msg("attrition:"); print(dcast(attr_rows, rule ~ sample_id, value.var = "n", fill = 0))

params <- list(
  script = "09b_spatial_qc.R", timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  seed = 1L, samples = SAMPLES, netflat = NETFLAT,
  fov_px = FOV_PX, px_um = PX_UM,
  edge = list(tol = EDGE_TOL, bin = EDGE_BIN, plateau_from = EDGE_PLATEAU,
              mode = EDGE_CUT_MODE, fits = edge_fits,
              applied = if (is.null(EDGE_CUTS)) "per_sample" else as.list(EDGE_CUTS)),
  staining = list(scope = STAINING_SCOPE, channels = CHANNELS,
                  i_frac_cut = I_FRAC_CUT, formula = deparse(form),
                  coefficients = coef(gfit), p_cut = P_CUT,
                  youden_J = max(you), cv_k = CV_K,
                  cv_error = mean(cv$err, na.rm = TRUE),
                  cv_sens = mean(cv$sens, na.rm = TRUE),
                  cv_spec = mean(cv$spec, na.rm = TRUE)),
  verdict = list(spatial_only = "qc_spatial_passing",
                 recommended  = "qc_passing_all (spatial AND qcCellsPassed)"),
  overcorrelation = list(corr_cut = CORR_CUT, dens_cut = DENS_CUT, n_min = CORR_NMIN,
                         note = "decision on |rho| only; analytic p reported but n-dependent"))
saveRDS(params, file.path(OUT_DIR, "qc_spatial_params.rds"))
writeLines(jsonlite::toJSON(params, auto_unbox = TRUE, digits = 6, pretty = TRUE),
           file.path(OUT_DIR, "qc_spatial_params.json"))

# ------------------------------------------------------------- 9. FIGURES ------
msg("figures")
fov_pos <- md[, .(xmin = min(CenterX_global_px), xmax = max(CenterX_global_px),
                  ymin = min(CenterY_global_px), ymax = max(CenterY_global_px)),
              by = .(sample_id, fov)]
fv <- merge(fov_verdict, fov_pos, by = c("sample_id", "fov"))
rect <- function(dt, fill, title, div = FALSE, lim = NULL) {
  p <- ggplot(dt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                      fill = .data[[fill]])) +
    geom_rect() + coord_fixed() + facet_wrap(~ sample_id) +
    labs(title = title, x = NULL, y = NULL)
  if (div) p + scale_fill_gradientn(colours = c("blue", "grey90", "red"), limits = lim)
  else     p + scale_fill_viridis_c()
}

save_fig(rect(fv, "frac_i", "cluster-i fraction per FOV"), "fov_frac_i.png", 11, 6)
save_fig(rect(fv, "max_abs_rho", "max |rho| nCount ~ local px (post-cut)"),
         "fov_overcorrelation.png", 11, 6)
save_fig(ggplot(fv, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                        fill = qc_spatial_staining)) +
           geom_rect() + coord_fixed() + facet_wrap(~ sample_id) +
           labs(title = "staining verdict (IB1 in scope; RB1 not_evaluated)"),
         "fov_staining_verdict.png", 11, 6)

chan_map <- chan_long[stat == "Mean"]
chan_map <- merge(chan_map, fov_pos, by = c("sample_id", "fov"))
chan_map[, z_cap := pmax(pmin(z_signal, 4), -4)]
save_fig(ggplot(chan_map, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = z_cap)) +
           geom_rect() + coord_fixed() + facet_grid(channel ~ sample_id) +
           scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-4, 4)) +
           labs(title = "per-FOV mean channel signal, z WITHIN sample"),
         "fov_channels_z.png", 10, 12)

ep <- merge(edge_prof, edge_fits, by = c("sample_id", "edge"))
save_fig(ggplot(ep, aes(d, y)) + geom_line() + geom_point(size = .5) +
           geom_hline(aes(yintercept = plateau), lty = 3) +
           geom_hline(aes(yintercept = (1 - EDGE_TOL) * plateau), lty = 2, colour = "red") +
           geom_vline(aes(xintercept = cut_px), colour = "blue") +
           facet_grid(sample_id ~ edge, scales = "free_y") + xlim(0, 1200) +
           labs(title = sprintf("per-edge depth profile; blue = fitted cut at tol = %.2f", EDGE_TOL),
                x = "px from edge", y = "mean nCount_RNA"),
         "edge_profiles_fitted.png", 12, 6)

st <- fov_tab[, .(sample_id, frac_i, p_junk, z_Membrane, z_CD45, z_DAPI,
                  qc_spatial_staining)]
save_fig((ggplot(st, aes(z_Membrane, frac_i, colour = qc_spatial_staining)) +
            geom_point(size = .8) + geom_hline(yintercept = I_FRAC_CUT, lty = 2) +
            facet_wrap(~ sample_id) + labs(title = "frac_i vs Membrane z")) /
           (ggplot(st[!is.na(p_junk)], aes(p_junk, frac_i, colour = qc_spatial_staining)) +
              geom_point(size = .8) + geom_vline(xintercept = P_CUT, lty = 2) +
              geom_hline(yintercept = I_FRAC_CUT, lty = 2) +
              labs(title = "fitted junk probability vs frac_i (IB1)")),
         "staining_rule.png", 10, 8)

# RB1 control panel: where RB1 falls relative to IB1's fitted boundary
ctrl <- copy(fov_tab)
ctrl[, p_all := predict(gfit, newdata = ctrl, type = "response")]
save_fig(ggplot(ctrl, aes(p_all, frac_i, colour = sample_id)) +
           geom_point(size = .9, alpha = .8) +
           geom_vline(xintercept = P_CUT, lty = 2) +
           geom_hline(yintercept = I_FRAC_CUT, lty = 2) +
           labs(title = "control: IB1's fitted rule scored on BOTH samples",
                subtitle = "RB1 is exempt by design; this shows whether the exemption is justified",
                x = "fitted junk probability (IB1 model)", y = "frac_i"),
         "control_rb1_vs_ib1_boundary.png", 8, 6)

aw <- dcast(attr_rows, rule ~ sample_id, value.var = "n", fill = 0)
awl <- melt(aw, id.vars = "rule", variable.name = "sample_id", value.name = "n")
save_fig(ggplot(awl, aes(reorder(rule, n), n, fill = sample_id)) +
           geom_col(position = "dodge") + coord_flip() +
           labs(title = "attrition by rule", x = NULL),
         "attrition.png", 8, 5)

msg("done ->", OUT_DIR)

library(magrittr)
library(tidyverse)
fov_tab$qc_spatial_staining %>% table
md$qc_spatial_fovpx %>% table
md$qc_passing_all
md$qcCellsPassed
md$qc_spatial_i_fraction %>% table
md$qc_spatial_staining

md_fovpx = md %>% group_by(fov, qc_spatial_fovpx, sample_id) %>% summarise(N = length(qc_spatial_fovpx))
md_fovpx = md_fovpx %>% group_by(fov, sample_id) %>% mutate(total = sum(N))
md_fovpx = md_fovpx %>% mutate(fraction = N / total)

md_fov_rect =md %>% group_by(sample_id, fov) %>% summarise(xmin = min(CenterX_global_px), xmax = max(CenterX_global_px), ymin = min(CenterY_global_px), ymax = max(CenterY_global_px))
md_fovpx = merge(md_fovpx, md_fov_rect)

md_fovpx$stat_group = "fraction"
md_fovpx$stat_group[md_fovpx$qc_spatial_fovpx == "ok"] = "ok"
md_fovpx$stat_group[grepl("FOV", md_fovpx$qc_spatial_fovpx)] = "FOV"

fov_discard = md_fovpx %>% filter(stat_group == "FOV")
fov_discard = select(fov_discard, fov, sample_id)

p_part = ggplot(md_fovpx %>% filter(stat_group == "fraction"), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fraction)) +
  geom_rect() +
  facet_grid(qc_spatial_fovpx~sample_id) +
  coord_fixed() +
  labs(title = "FOV partially filtered")
p_part
save_fig(p_part, "fov_qc_position_filters.png", w = 5, h = 10)

p_ok = ggplot(md_fovpx %>% filter(stat_group == "ok"), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fraction)) +
  geom_rect() +
  facet_grid(qc_spatial_fovpx~sample_id) +
  coord_fixed() +
  labs(title = "Fraction OK by FOV")
save_fig(p_ok, "fov_qc_ok.png", w = 6, h = 3.4)

p_discard = ggplot(md_fovpx %>% filter(stat_group == "FOV"), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = qc_spatial_fovpx)) +
  annotate("rect", xmin = md_fov_rect$xmin, xmax = md_fov_rect$xmax, ymin = md_fov_rect$ymin, ymax = md_fov_rect$ymax, fill = "gray") +
  geom_rect() +
  facet_grid(.~sample_id) +
  coord_fixed() +
  labs(title = "FOV completely discarded")
save_fig(p_discard, "fov_total_discard.png", w = 5, h = 10)

