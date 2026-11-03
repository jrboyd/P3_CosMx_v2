## ============================================================================
## Step 28 -- Bright chemokine cells: where they are, whether they are the HRS-like
## ones, and whether IL13 / BATF3 / CD30 co-occur in the same cells.
##
## WHY. The step-16 hotspots are 100 um bins at IB1's 99th percentile of POOLED
## CCL17+CCL22 cp10k with >= 40 cells. A ~58-cell bin carries ~19,000 counts, so one
## cell with 100 chemokine counts adds ~52 cp10k -- under half the cut -- while two or
## three moderate cells add as much. The rule is a DENSITY detector: it finds clusters
## of moderate cells and dilutes isolated bright ones (a 1 mm preview around
## IB1_55_162 had 10.4% of the square's chemokine inside its 3 hotspot boxes). Step 27
## found IB1's chemokine sources carry an HRS-compatible excess over RB1's (IL13 4.5x,
## BATF3 3.5x, CCR7 2.8x) but as a population average over 23% of hotspot cells, far
## more than HRS cells can be. This step works per cell, sample-wide.
##
## DESIGN (fixed before running)
##   cells    matched tier, defect-free, both samples, all regions
##            (output/16_regions/matched/cell_regions.csv.gz)
##   tiers    CCL17+CCL22 counts: <10 (reference), 10-19, 20-49, 50-99, 100+
##   depth    nCount_RNA MINUS the chemokine counts (a bright cell's own chemokine
##            would otherwise inflate its depth and deflate every excess)
##   part 1   per tier: share of cells and molecules, area, how many in hotspots;
##            per hotspot: share of its chemokine from cells >= 50
##   part 2   per tier, per gene: obs / exp vs the <10 cells of the same sample in the
##            same depth decile; FOV cluster bootstrap; contrast 100+ vs 10-19
##   part 3   among sources (>= 10): co-positivity (>= 1 count) of IL13, BATF3,
##            TNFRSF8 pairs and triple, O/E against independence within depth-quintile
##            x tier strata. NULL = random panel gene pairs/triples with each member's
##            detection rate matched within 1.25x, same O/E -- whatever co-occurrence
##            they share is generic (size, activation, leakage). Controls CCR7-LAMP3,
##            CD79A-MS4A1, CD3D-CD3E. FOV bootstrap for the observed O/E.
##
## CAVEATS carried: leakage can put an HRS cell's transcripts into its neighbour, so
## co-occurrence localises a programme to a cell OR its immediate contact; the region
## remains the unit of inference and n = 1 patient.
##
##   sbatch --mem=48G --time=2:00:00 scripts/sbatch/run_rscript.sh \
##       scripts/28_script_bright_sources.R [smoke]
## ============================================================================
suppressPackageStartupMessages(library(data.table))

args   <- commandArgs(trailingOnly = TRUE)
SMOKE  <- "smoke" %in% args
FLAT   <- "input/CosMx/export1/flatFiles"
REG    <- "output/16_regions/matched/cell_regions.csv.gz"
OUT    <- if (SMOKE) "output/28_bright_sources_smoke" else "output/28_bright_sources"
N_BOOT <- if (SMOKE) 20L else 1000L
N_NULL <- if (SMOKE) 50L else 1000L
SRC_MIN <- 10L
TIER_BREAKS <- c(-Inf, 9.5, 19.5, 49.5, 99.5, Inf)
TIERS <- c("<10", "10-19", "20-49", "50-99", "100+")
MATCH_FOLD <- 1.25
set.seed(28)

SETS <- list(
  chemokine = c("CCL17", "CCL22"),
  DC        = c("LAMP3", "CCR7", "CCL19", "IDO1", "CD1C", "CLEC9A", "XCR1", "CLEC10A",
                "FCER1A", "CD1A", "CD207", "MARCKSL1"),
  HRS       = c("TNFRSF8", "FUT4", "IL13", "LTA", "CSF2"),
  ambiguous = c("CD274", "PDCD1LG2", "CD40", "BATF3", "IRF4", "CD80", "CD86", "CD83", "JAK2"),
  macro     = c("CD68", "CD163", "C1QA", "CD14", "LYZ"),
  Bprog     = c("CD19", "CD79A", "MS4A1"),
  T         = c("CD3E", "CD3D"))
GENES  <- unique(unlist(SETS))
TEST   <- setdiff(GENES, SETS$chemokine)
SET_OF <- setNames(rep(names(SETS), lengths(SETS)), unlist(SETS))
HRS3   <- c("IL13", "BATF3", "TNFRSF8")
PAIRS  <- list(c("IL13", "BATF3"), c("IL13", "TNFRSF8"), c("BATF3", "TNFRSF8"),
               c("CCR7", "LAMP3"), c("CD79A", "MS4A1"), c("CD3D", "CD3E"))
## c_2_ is IB1, c_1_ is RB1 (CLAUDE.md). Used only to BUILD ids from each sample's own
## file; the join then checks `samp` agrees.
SLIDE <- c(IB1 = 2L, RB1 = 1L)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
r3 <- function(x) round(x, 3)

## ---- 0. cells and panel counts ---------------------------------------------
reg <- fread(cmd = paste("zcat", REG),
             select = c("cell_id", "samp", "fov_id", "region", "bx", "by_", "is_defect",
                        "tier_matched", "nCount_RNA", "Area.um2"))
reg <- reg[tier_matched == TRUE & is_defect == FALSE]

exp_file <- function(s) file.path(FLAT, s, paste0(s, "_exprMat_file.csv.gz"))
read_panel <- function(s) {
  hdr <- names(fread(cmd = sprintf("zcat %s | head -1", exp_file(s))))
  stopifnot(all(GENES %in% hdr))
  cols <- match(c("fov", "cell_ID", GENES), hdr)
  ex <- fread(cmd = sprintf("zcat %s | cut -d, -f%s%s", exp_file(s),
                            paste(cols, collapse = ","), if (SMOKE) " | head -300001" else ""))
  ex[, `:=`(cell_id = sprintf("c_%d_%d_%d", SLIDE[[s]], fov, cell_ID), samp_file = s)]
  ex
}
d <- merge(reg, rbindlist(lapply(names(SLIDE), read_panel)), by = "cell_id")
stopifnot(all(d$samp == d$samp_file))
if (!SMOKE) stopifnot(nrow(d) == nrow(reg))
cat(sprintf("cells: %s\n", paste(d[, .N, by = samp][, sprintf("%s %d", samp, N)], collapse = ", ")))
d[, chemo := CCL17 + CCL22]
d[, depth0 := nCount_RNA - chemo]
d <- d[depth0 > 0]
d[, tier := factor(cut(chemo, TIER_BREAKS, labels = TIERS), levels = TIERS)]
d[, depth_q := cut(depth0, unique(quantile(depth0, 0:10 / 10)), include.lowest = TRUE,
                   labels = FALSE), by = samp]
d[, hot := region == "cHL_hotspot"]

## ---- 1. where the bright cells are -------------------------------------------
p1 <- d[, .(cells = .N, molecules = sum(chemo), in_hotspot = sum(hot),
            area_med = r3(median(Area.um2)), depth0_med = as.numeric(median(depth0))),
        keyby = .(samp, tier)]
p1[, `:=`(pct_cells = r3(100 * cells / sum(cells)), pct_molecules = r3(100 * molecules / sum(molecules)),
          pct_in_hotspot = r3(100 * in_hotspot / cells)), by = samp]
fwrite(p1, file.path(OUT, "tier_summary.csv"))
hs <- d[samp == "IB1" & hot == TRUE,
        .(cells = .N, chemo = sum(chemo), n50 = sum(chemo >= 50), n100 = sum(chemo >= 100),
          n10_49 = sum(chemo >= 10 & chemo < 50),
          share_from_50plus = r3(sum(chemo[chemo >= 50]) / max(sum(chemo), 1))), by = .(bx, by_)]
fwrite(hs, file.path(OUT, "hotspot_bright_cells.csv"))
cat("\n== part 1: tiers\n"); print(p1)
cat(sprintf("\nIB1 hotspots (n = %d): %d have NO cell >= 50, %d none >= 100; median share of a hotspot's chemokine from cells >= 50: %.2f\n",
            nrow(hs), hs[n50 == 0, .N], hs[n100 == 0, .N], median(hs$share_from_50plus)))
cat(sprintf("IB1 cells >= 50: %d, of which %d (%.1f%%) inside a hotspot; cells >= 100: %d, %d inside\n",
            d[samp == "IB1" & chemo >= 50, .N], d[samp == "IB1" & chemo >= 50 & hot, .N],
            100 * d[samp == "IB1" & chemo >= 50, mean(hot)],
            d[samp == "IB1" & chemo >= 100, .N], d[samp == "IB1" & chemo >= 100 & hot, .N]))

## ---- 2. per-tier excess, FOV bootstrap -----------------------------------------
G <- c("dep", TEST)
agg <- d[, c(list(dep = sum(as.numeric(depth0))), lapply(.SD, function(z) sum(as.numeric(z)))),
         by = .(samp, fov_id, tier, depth_q), .SDcols = TEST]

tier_lr <- function(A, w) {          # A: one sample's agg; w: data.table(fov_id, w)
  X <- A[w, on = "fov_id", nomatch = 0]
  X[, (G) := lapply(.SD, `*`, w), .SDcols = G]
  ref  <- X[tier == "<10", lapply(.SD, sum), by = depth_q, .SDcols = G]
  rate <- as.matrix(ref[, ..TEST]) / ref$dep
  rownames(rate) <- ref$depth_q
  out <- lapply(c(TIERS[-1], "10+"), function(t) {
    tg <- if (t == "10+") X[tier != "<10"] else X[tier == t]
    tg <- tg[, lapply(.SD, sum), by = depth_q, .SDcols = G][as.character(depth_q) %in% rownames(rate)]
    ## An EMPTY tier (RB1 has no 100+ cell) must be NA, not log(0.5/0.5) = 0 -- that
    ## reads as "no excess" and fakes every contrast against it.
    if (!nrow(tg)) return(setNames(rep(NA_real_, length(TEST)), TEST))
    o <- colSums(as.matrix(tg[, ..TEST]))
    e <- colSums(tg$dep * rate[as.character(tg$depth_q), , drop = FALSE])
    log((o + 0.5) / (e + 0.5))
  })
  M <- do.call(rbind, out); rownames(M) <- c(TIERS[-1], "10+"); M
}
pboot <- function(x) { x <- x[is.finite(x)]; n <- length(x)
  min(1, 2 * min((1 + sum(x <= 0)) / (1 + n), (1 + sum(x >= 0)) / (1 + n))) }

ex_rows <- list(); top_rows <- list()
for (s in names(SLIDE)) {
  A <- agg[samp == s]
  fovs <- unique(A$fov_id)
  M0 <- tier_lr(A, data.table(fov_id = fovs, w = 1))
  obs_tab <- d[samp == s, c(list(n = .N), lapply(.SD, sum)), by = tier, .SDcols = TEST]
  BB <- array(NA_real_, c(N_BOOT, nrow(M0), ncol(M0)), dimnames = list(NULL, rownames(M0), colnames(M0)))
  for (b in seq_len(N_BOOT)) {
    w <- data.table(fov_id = sample(fovs, replace = TRUE))[, .(w = .N), by = fov_id]
    BB[b, , ] <- tier_lr(A, w)
  }
  for (t in rownames(M0)) {
    se <- apply(BB[, t, ], 2, sd, na.rm = TRUE)
    n_t <- if (t == "10+") d[samp == s & tier != "<10", .N] else d[samp == s & tier == t, .N]
    ex_rows[[paste(s, t)]] <- data.table(samp = s, tier = t, n_cells = n_t, gene = colnames(M0),
      excess = exp(M0[t, ]), lo = exp(apply(BB[, t, ], 2, quantile, .025, na.rm = TRUE)),
      hi = exp(apply(BB[, t, ], 2, quantile, .975, na.rm = TRUE)),
      p_boot = apply(BB[, t, ], 2, pboot), p_wald = 2 * pnorm(-abs(M0[t, ] / se)))
  }
  dl <- BB[, "100+", ] - BB[, "10-19", ]
  d0 <- M0["100+", ] - M0["10-19", ]
  top_rows[[s]] <- data.table(samp = s, gene = colnames(M0), ratio_100_vs_10 = exp(d0),
    lo = exp(apply(dl, 2, quantile, .025, na.rm = TRUE)), hi = exp(apply(dl, 2, quantile, .975, na.rm = TRUE)),
    p_boot = apply(dl, 2, pboot), p_wald = 2 * pnorm(-abs(d0 / apply(dl, 2, sd, na.rm = TRUE))))
}
ex <- rbindlist(ex_rows); ex[, set := SET_OF[gene]]
ex[, q := p.adjust(p_wald, "BH"), by = .(samp, tier)]
top <- rbindlist(top_rows); top[, set := SET_OF[gene]]
top[, q := p.adjust(p_wald, "BH"), by = samp]
fwrite(ex[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 4) else z)], file.path(OUT, "tier_excess.csv"))
fwrite(top[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 4) else z)], file.path(OUT, "tier_top_vs_low.csv"))
cat("\n== part 2: excess by tier (obs/exp vs <10 cells, depth-matched)\n")
print(dcast(ex[gene %in% c(HRS3, "IRF4", "CCR7", "LAMP3", "XCR1", "CD1C", "IDO1", "PDCD1LG2",
                           "CD274", "FUT4", "CD80", "LYZ", "CD79A")],
            samp + gene ~ tier, value.var = "excess", fun.aggregate = function(x) round(x, 2))[order(samp, gene)])
cat("\n100+ vs 10-19, q < 0.05:\n")
print(top[q < 0.05, .(samp, set, gene, ratio = round(ratio_100_vs_10, 2), lo = round(lo, 2),
                      hi = round(hi, 2), q = signif(q, 2))][order(samp, -ratio)])

## ---- 3. co-occurrence among sources -------------------------------------------
src <- d[chemo >= SRC_MIN]
src[, sq := cut(depth0, unique(quantile(depth0, 0:5 / 5)), include.lowest = TRUE, labels = FALSE), by = samp]
src[, stratum := paste(sq, tier)]
## genome-wide counts for the source cells only (needed for the matched null)
read_src_all <- function(s) {
  keys <- file.path(OUT, sprintf("src_keys_%s.csv", s))
  fwrite(src[samp == s, .(fov, cell_ID)], keys, col.names = FALSE)
  cmd <- sprintf("zcat %s | awk -F, 'NR==FNR{k[$1\",\"$2]; next} FNR==1 || (($1\",\"$2) in k)' %s -",
                 exp_file(s), keys)
  x <- fread(cmd = cmd)
  x[, cell_id := sprintf("c_%d_%d_%d", SLIDE[[s]], fov, cell_ID)]
  file.remove(keys)
  x
}
co_rows <- list(); null_rows <- list(); kof3 <- list()
for (s in names(SLIDE)) {
  S <- src[samp == s]
  if (nrow(S) < 50) { cat("too few sources in", s, "\n"); next }
  X <- read_src_all(s)
  X <- X[match(S$cell_id, X$cell_id)]
  stopifnot(!anyNA(X$cell_id))
  allg <- setdiff(names(X), c("fov", "cell_ID", "cell_id"))
  allg <- allg[!grepl("^(Negative|NegPrb|SystemControl|FalseCode)", allg) &
               !allg %in% SETS$chemokine]
  P <- as.matrix(X[, ..allg]) > 0
  det <- colMeans(P)
  st <- S$stratum
  ust <- unique(st); sidx <- match(st, ust)
  ## O/E for a gene combination over a set of rows (a bootstrap draw may repeat rows).
  ## Per-stratum positive rates are re-estimated from THOSE rows, so a resample or a
  ## tier subset is scored against its own independence expectation.
  oe <- function(gs, rows = seq_len(nrow(P))) {
    Pr <- P[rows, gs, drop = FALSE]
    o  <- sum(rowSums(Pr) == length(gs))
    sr <- tabulate(sidx[rows], length(ust))
    ## rowsum() drops strata absent from `rows`, so place its rows into a full matrix
    tmp <- rowsum(Pr * 1, sidx[rows])
    pr <- matrix(0, length(ust), length(gs))
    pr[as.integer(rownames(tmp)), ] <- tmp
    pr <- pr / pmax(sr, 1)
    e  <- sum(sr * apply(pr, 1, prod))
    c(o = o, e = e, oe = (o + 0.5) / (e + 0.5))
  }
  pool_for <- function(g) { pool <- names(det)[det >= det[g] / MATCH_FOLD & det <= det[g] * MATCH_FOLD]
    setdiff(pool, c(HRS3, g)) }
  fovs <- unique(S$fov_id); fidx <- split(seq_len(nrow(S)), S$fov_id)
  combos <- c(PAIRS, list(HRS3))
  for (gs in combos) {
    if (!all(gs %in% allg)) next
    obs <- oe(gs)
    pools <- lapply(gs, pool_for)
    nul <- replicate(N_NULL, { pick <- vapply(pools, function(p) sample(p, 1), "")
      if (anyDuplicated(pick)) NA_real_ else oe(pick)[["oe"]] })
    nul <- nul[is.finite(nul)]
    bt <- replicate(N_BOOT, { rows <- unlist(fidx[sample(fovs, replace = TRUE)], use.names = FALSE)
      oe(gs, rows)[["oe"]] })
    co_rows[[paste(s, paste(gs, collapse = "+"))]] <- data.table(
      samp = s, combo = paste(gs, collapse = "+"), n_sources = nrow(S),
      det = paste(sprintf("%.3f", det[gs]), collapse = "/"),
      obs = obs[["o"]], exp = r3(obs[["e"]]), oe = r3(obs[["oe"]]),
      oe_lo = r3(quantile(bt, .025)), oe_hi = r3(quantile(bt, .975)),
      null_median = r3(median(nul)), null_p95 = r3(quantile(nul, .95)),
      p_null = signif((1 + sum(nul >= obs[["oe"]])) / (1 + length(nul)), 3),
      n_pool = paste(lengths(pools), collapse = "/"))
    null_rows[[paste(s, paste(gs, collapse = "+"))]] <- data.table(samp = s,
      combo = paste(gs, collapse = "+"), null_oe = nul)
  }
  ## how many of the three each source carries, by tier
  S[, k := as.integer(P[, "IL13"]) + P[, "BATF3"] + P[, "TNFRSF8"]]
  kof3[[s]] <- S[, .(n = .N, k0 = r3(mean(k == 0)), k1 = r3(mean(k == 1)),
                     k2 = r3(mean(k == 2)), k3 = r3(mean(k == 3)),
                     area_med = r3(median(Area.um2))), keyby = .(samp, tier)]
  ## and co-occurrence within the bright tiers alone (>= 50) vs moderate (10-49)
  for (lab in c("10-49", "50+")) {
    rows <- which(if (lab == "50+") S$chemo >= 50 else S$chemo < 50)
    for (gs in c(PAIRS[1:3], list(HRS3))) {
      v <- oe(gs, rows)
      co_rows[[paste(s, lab, paste(gs, collapse = "+"))]] <- data.table(samp = s,
        combo = paste0(paste(gs, collapse = "+"), " [", lab, "]"), n_sources = length(rows),
        det = paste(sprintf("%.3f", colMeans(P[rows, gs, drop = FALSE])), collapse = "/"),
        obs = v[["o"]], exp = r3(v[["e"]]), oe = r3(v[["oe"]]))
    }
  }
  rm(X, P); gc()
}
co <- rbindlist(co_rows, fill = TRUE)
fwrite(co, file.path(OUT, "cooccurrence.csv"))
fwrite(rbindlist(null_rows), file.path(OUT, "cooccurrence_null.csv.gz"))
k3 <- rbindlist(kof3)
fwrite(k3, file.path(OUT, "hrs3_k_of_3_by_tier.csv"))
cat("\n== part 3: co-occurrence among sources (O/E vs independence within depth x tier)\n")
print(co)
cat("\nIL13/BATF3/TNFRSF8: share of sources positive for k of 3, by tier\n")
print(k3)
cat("\ndone:", OUT, "\n")
