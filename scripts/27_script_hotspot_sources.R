## ============================================================================
## Step 27 -- Who makes the chemokine in the Hodgkin hotspots: HRS or dendritic cells?
##
## The hotspots are defined by CCL17+CCL22, which HRS cells secrete -- but so do
## mature dendritic cells, and the hotspots are DC territory (R4). Per-cell HRS calling
## is closed (step 10), so this asks a narrower question: in the cells that carry the
## chemokine ("sources"), which programme is in excess over their own neighbours?
##
## A one-window preview (IB1_55_162, 27 sources) found a shared HRS/activated-DC
## programme up (BATF3 5.9x, PDCD1LG2, CCR7, CD83, CD274), with the genes that would
## split the two -- DC-specific LAMP3/CD1C/XCR1, HRS-specific TNFRSF8/FUT4/IRF4 --
## flat or down, all on a handful of counts. This pools every hotspot.
##
## DESIGN (fixed before running; gene sets and thresholds are the preview's, unchanged)
##   cells    defect-free (`is_defect == FALSE`), matched tier, from
##            output/16_regions/matched/cell_regions.csv.gz
##   units    IB1: the 80 step-16 hotspots (100 um bins, `region == "cHL_hotspot"`)
##            RB1: its own 105 top-1% chemokine bins (`chemo_top`), the construction-
##            identical control. RB1 has no Hodgkin component, so its sources CANNOT be
##            HRS cells: whatever they look like is what a chemokine-making non-HRS
##            cell looks like on this panel.
##   source   CCL17 + CCL22 >= SRC_MIN counts
##   excess   per gene, observed source counts / expected, where expected uses the
##            NON-source cells of the same sample's unit set, depth-matched within
##            nCount deciles (rate per count in the decile x the source cell's depth).
##   inference unit-level bootstrap (hotspots / bins resampled), percentile CI, p with
##            the (1+count)/(1+n) correction plus a Wald p on the log ratio, BH within
##            sample. The readout is the IB1-minus-RB1 log excess per gene: within-
##            sample ratios, so the between-slide offset (fact 5) enters only through
##            how each slide's sources differ from their OWN neighbours.
##   sensitivity  seam-free hotspots only; SRC_MIN 5 and 20.
##
## CAVEATS carried, not solved: transcripts leak into neighbouring cells (fact 6, step
## 18), so a source's excess can be partly its neighbour's; the HRS-specific genes sit
## at the detection floor everywhere (TNFRSF8 max 6 counts), so their silence is weak.
##
##   sbatch --mem=32G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##       scripts/27_script_hotspot_sources.R [smoke]
## ============================================================================
suppressPackageStartupMessages(library(data.table))

args  <- commandArgs(trailingOnly = TRUE)
SMOKE <- "smoke" %in% args
FLAT  <- "input/CosMx/export1/flatFiles"
REG   <- "output/16_regions/matched/cell_regions.csv.gz"
OUT   <- if (SMOKE) "output/27_hotspot_sources_smoke" else "output/27_hotspot_sources"
N_BOOT  <- if (SMOKE) 50L else 2000L
SRC_MIN <- 10L
SRC_SENS <- c(5L, 20L)
N_DEPTH  <- 10L
set.seed(27)

SETS <- list(
  chemokine = c("CCL17", "CCL22"),
  DC        = c("LAMP3", "CCR7", "CCL19", "IDO1", "CD1C", "CLEC9A", "XCR1", "CLEC10A",
                "FCER1A", "CD1A", "CD207", "MARCKSL1"),
  HRS       = c("TNFRSF8", "FUT4", "IL13", "LTA", "CSF2"),
  ambiguous = c("CD274", "PDCD1LG2", "CD40", "BATF3", "IRF4", "CD80", "CD86", "CD83", "JAK2"),
  macro     = c("CD68", "CD163", "C1QA", "CD14", "LYZ"),
  Bprog     = c("CD19", "CD79A", "MS4A1"),
  T         = c("CD3E"))
GENES <- unique(unlist(SETS))
TEST  <- setdiff(GENES, SETS$chemokine)          # chemokine defines the sources: circular
SET_OF <- setNames(rep(names(SETS), lengths(SETS)), unlist(SETS))
## The cell-id slide token is c_2_ for IB1 and c_1_ for RB1 (CLAUDE.md). Used only to
## BUILD ids from each sample's own file; the join below then checks `samp` agrees.
SLIDE <- c(IB1 = 2L, RB1 = 1L)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- 1. units and cells ---------------------------------------------------
reg <- fread(cmd = paste("zcat", REG),
             select = c("cell_id", "samp", "region", "bx", "by_", "is_defect",
                        "tier_matched", "chemo_top", "seam_flag", "nCount_RNA", "Area.um2"))
reg <- reg[tier_matched == TRUE & is_defect == FALSE &
           ((samp == "IB1" & region == "cHL_hotspot") | (samp == "RB1" & chemo_top == TRUE))]
reg[, unit := sprintf("%s_%d_%d", samp, bx, by_)]
cat(sprintf("units: IB1 %d hotspots (%d seam-free), RB1 %d top-1%% bins\n",
            reg[samp == "IB1", uniqueN(unit)],
            reg[samp == "IB1" & seam_flag == FALSE, uniqueN(unit)],
            reg[samp == "RB1", uniqueN(unit)]))
stopifnot(reg[samp == "IB1", uniqueN(unit)] == 80L)

## ---- 2. counts for the panel genes, per sample -------------------------------
read_counts <- function(s) {
  f   <- file.path(FLAT, s, paste0(s, "_exprMat_file.csv.gz"))
  hdr <- names(fread(cmd = sprintf("zcat %s | head -1", f)))
  stopifnot(all(GENES %in% hdr))
  cols <- match(c("fov", "cell_ID", GENES), hdr)
  ## The trailing head -n bounds the smoke run; it is a no-op cap otherwise.
  cmd <- sprintf("zcat %s | cut -d, -f%s%s", f, paste(cols, collapse = ","),
                 if (SMOKE) " | head -400001" else "")
  ex <- fread(cmd = cmd)
  ex[, cell_id := sprintf("c_%d_%d_%d", SLIDE[[s]], fov, cell_ID)]
  ex[, c("fov", "cell_ID") := NULL]
  ex[, samp_file := s]
  ex
}
ex <- rbindlist(lapply(names(SLIDE), read_counts))
d  <- merge(reg, ex, by = "cell_id")
stopifnot(all(d$samp == d$samp_file))            # id construction did not swap samples
cat(sprintf("joined %d of %d unit cells to counts%s\n", nrow(d), nrow(reg),
            if (SMOKE) " (SMOKE: partial matrix)" else ""))
if (!SMOKE) stopifnot(nrow(d) == nrow(reg))
d[, chemo := CCL17 + CCL22]
d[, depth_q := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:N_DEPTH / N_DEPTH)),
                   include.lowest = TRUE, labels = FALSE), by = samp]

## ---- 3. the excess statistic ----------------------------------------------
## For a cell set D (one sample, some units), source threshold k:
##   obs_g = sum over sources of count_g
##   exp_g = sum over sources of nCount * (rate_g among non-sources in that depth decile)
## Returns per-gene obs, exp, and log(obs/exp) with +0.5 continuity on both.
excess <- function(D, k) {
  D[, src := chemo >= k]
  ref <- D[src == FALSE, c(list(n = sum(nCount_RNA)), lapply(.SD, sum)),
           by = depth_q, .SDcols = TEST]
  srcd <- D[src == TRUE, .(n = sum(nCount_RNA)), by = depth_q]
  m <- merge(srcd, ref, by = "depth_q", suffixes = c("_src", "_ref"))
  exp_g <- vapply(TEST, function(g) sum(m$n_src * m[[g]] / m$n_ref), numeric(1))
  obs_g <- vapply(TEST, function(g) sum(D[src == TRUE][[g]]), numeric(1))
  data.table(gene = TEST, obs = obs_g, exp = exp_g,
             lr = log((obs_g + 0.5) / (exp_g + 0.5)))
}

boot_excess <- function(D, k) {
  units <- unique(D$unit)
  idx <- split(seq_len(nrow(D)), D$unit)
  B <- matrix(NA_real_, N_BOOT, length(TEST), dimnames = list(NULL, TEST))
  for (b in seq_len(N_BOOT)) {
    pick <- sample(units, replace = TRUE)
    Db <- D[unlist(idx[pick], use.names = FALSE)]
    B[b, ] <- excess(Db, k)$lr
  }
  B
}

pboot <- function(x) {                            # two-sided, corrected (Conventions)
  n <- sum(is.finite(x)); x <- x[is.finite(x)]
  min(1, 2 * min((1 + sum(x <= 0)) / (1 + n), (1 + sum(x >= 0)) / (1 + n)))
}

run <- function(D, k, label) {
  out <- list(); boots <- list()
  for (s in c("IB1", "RB1")) {
    Ds <- D[samp == s]
    e  <- excess(copy(Ds), k)
    B  <- boot_excess(Ds, k)
    boots[[s]] <- B
    e[, `:=`(samp = s, n_units = uniqueN(Ds$unit), n_src = sum(Ds$chemo >= k),
             n_src_units = uniqueN(Ds[chemo >= k]$unit),
             lo = apply(B, 2, quantile, 0.025, na.rm = TRUE),
             hi = apply(B, 2, quantile, 0.975, na.rm = TRUE),
             se = apply(B, 2, sd, na.rm = TRUE),
             p_boot = apply(B, 2, pboot))]
    e[, p_wald := 2 * pnorm(-abs(lr / se))]
    e[, q := p.adjust(p_wald, "BH")]
    out[[s]] <- e
  }
  res <- rbindlist(out)
  ## IB1 minus RB1: independent resamples, so difference the bootstrap draws directly.
  dB <- boots$IB1 - boots$RB1
  dif <- merge(res[samp == "IB1", .(gene, lr_IB1 = lr)], res[samp == "RB1", .(gene, lr_RB1 = lr)],
               by = "gene")
  dif[, `:=`(d_lr = lr_IB1 - lr_RB1,
             d_lo = apply(dB, 2, quantile, 0.025, na.rm = TRUE)[gene],
             d_hi = apply(dB, 2, quantile, 0.975, na.rm = TRUE)[gene],
             d_se = apply(dB, 2, sd, na.rm = TRUE)[gene],
             d_p_boot = apply(dB, 2, pboot)[gene])]
  dif[, d_p_wald := 2 * pnorm(-abs(d_lr / d_se))]
  dif[, d_q := p.adjust(d_p_wald, "BH")]
  res[, analysis := label]; dif[, analysis := label]
  list(res = res, dif = dif)
}

## ---- 4. primary + sensitivity ------------------------------------------------
an <- list(primary = list(D = d, k = SRC_MIN),
           seam_free = list(D = d[samp == "RB1" | seam_flag == FALSE], k = SRC_MIN))
for (k in SRC_SENS) an[[sprintf("src_min_%d", k)]] <- list(D = d, k = k)
R <- lapply(names(an), function(a) { cat("running", a, "\n"); run(an[[a]]$D, an[[a]]$k, a) })
res <- rbindlist(lapply(R, `[[`, "res")); dif <- rbindlist(lapply(R, `[[`, "dif"))
res[, set := SET_OF[gene]]; dif[, set := SET_OF[gene]]
num <- function(x) round(x, 4)
fwrite(res[, .(analysis, samp, set, gene, obs, exp = num(exp), excess = num(exp(lr)),
               lo = num(exp(lo)), hi = num(exp(hi)), p_boot = signif(p_boot, 3),
               p_wald = signif(p_wald, 3), q = signif(q, 3),
               n_units, n_src_units, n_src)],
       file.path(OUT, "source_excess.csv"))
fwrite(dif[, .(analysis, set, gene, ratio_IB1_over_RB1 = num(exp(d_lr)),
               lo = num(exp(d_lo)), hi = num(exp(d_hi)), p_boot = signif(d_p_boot, 3),
               p_wald = signif(d_p_wald, 3), q = signif(d_q, 3))],
       file.path(OUT, "source_excess_IB1_vs_RB1.csv"))

## ---- 5. descriptives -----------------------------------------------------------
d[, src := chemo >= SRC_MIN]
desc <- d[, .(cells = .N, sources = sum(src),
              chemo_share_in_sources = round(sum(chemo[src]) / sum(chemo), 3),
              area_src = round(median(Area.um2[src]), 1), area_ref = round(median(Area.um2[!src]), 1),
              depth_src = as.numeric(median(nCount_RNA[src])),
              depth_ref = as.numeric(median(nCount_RNA[!src])),
              area_src_p90 = round(quantile(Area.um2[src], .9), 1),
              area_src_max = round(max(Area.um2[src]), 1)), by = samp]
fwrite(desc, file.path(OUT, "source_descriptives.csv"))
fwrite(d[src == TRUE, c("cell_id", "samp", "unit", "seam_flag", "nCount_RNA", "Area.um2",
                        "chemo", GENES), with = FALSE],
       file.path(OUT, "source_cells.csv.gz"))

cat("\n== descriptives\n"); print(desc)
cat("\n== primary: source excess (obs/exp), by sample\n")
w <- dcast(res[analysis == "primary"], set + gene ~ samp,
           value.var = c("obs", "lr", "q"))
w[, `:=`(ex_IB1 = round(exp(lr_IB1), 2), ex_RB1 = round(exp(lr_RB1), 2),
         q_IB1 = signif(q_IB1, 2), q_RB1 = signif(q_RB1, 2))]
print(w[order(set, -ex_IB1), .(set, gene, obs_IB1, ex_IB1, q_IB1, obs_RB1, ex_RB1, q_RB1)])
cat("\n== IB1 / RB1 ratio of excess, all analyses\n")
print(dcast(dif, set + gene ~ analysis, value.var = "d_lr",
            fun.aggregate = function(x) round(exp(x), 2))[order(set, gene)])
cat("\nprimary, q < 0.05 for IB1 vs RB1:\n")
print(dif[analysis == "primary" & d_q < 0.05,
          .(set, gene, ratio = round(exp(d_lr), 2), lo = round(exp(d_lo), 2),
            hi = round(exp(d_hi), 2), q = signif(d_q, 2))][order(-ratio)])
cat("\ndone:", OUT, "\n")
