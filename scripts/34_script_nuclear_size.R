## ============================================================================
## Step 34 -- Nuclear size and shape of chemokine-bright cells, from the vendor's DAPI
## segmentation (NucArea, NucAspectRatio), and whether the HRS-leaning genes track it.
##
## WHY. Across steps 28-33 the one morphology feature that has followed the candidates
## is the nucleus: bright cells are larger (area 58 -> 108 um^2 across tiers) and their
## DNA profile dips less at 5-10 um than size-matched controls' (steps 31-33), i.e. a
## bigger nucleus relative to the segmented cell. HRS cells are defined partly by that
## (large, often multilobed nuclei with huge nucleoli). CD45 profiles separated nothing
## in the RNA (step 33); this asks the same question of nuclear size.
##
## GATE. The vendor NucArea must agree with a nuclear area measured directly from the
##   DNA image INSIDE each cell's vendor polygon (scripts/34_extract_polygons.R +
##   34_nuclear_in_polygon.py, the 1,000 step-31 cells): Spearman >= 0.5 over cells whose
##   polygon holds a nucleus. Fail -> every downstream table is marked NOT INTERPRETABLE.
##   (A first, centroid-anchored image measure -- 34_nuclear_image_check.py -- was
##   abandoned: candidate centroids often sit on DNA-dim pixels between nuclei, so a region
##   grown from the centroid is a dim blob or a flood; 725 of 1,000 were merge-flagged.)
##   REVISION NOTE. The first pre-declared gate had two parts that the smoke run showed
##   were badly designed, and it was replaced BEFORE the full run: (a) "lymphocyte nuclei
##   are smaller than myeloid/fibro nuclei at AUC >= 0.65" -- in lymphoid tissue those
##   lineages' nuclei are similar (all 32-36 um^2 here) and the labels are 70% accurate,
##   so the premise, not the metric, failed (AUC 0.43-0.48); it is kept as a descriptive
##   table. (b) "NucArea vs the step-31 DNA profile at 5-10 um" -- a typical nucleus
##   (r ~ 3.3 um) does not reach that band, which reads cytoplasm vs neighbours' nuclei
##   instead (cell area correlated -0.47 with it). The image measurement replaces it.
##
## DESIGN
##   cells   matched tier, defect-free, both samples; tiers by CCL17+CCL22 as step 28
##   size    NucArea and Area in um^2 (px x 0.120281^2), nuc_frac = NucArea / Area
##   desc    by tier x sample; and IB1 >= 50 vs < 10 within cell-area deciles (is the
##           nucleus bigger than the cell size alone predicts?)
##   genes   IB1 sources (>= 10) in bands 10-19 / 20-49 / 50-99 / 100+; quintiles of
##           NucArea, nuc_frac and NucAspectRatio within band; per gene obs / exp vs < 10
##           cells of the same depth0 decile (step 28's reference); FOV cluster bootstrap;
##           top vs bottom quintile; BH over genes x bands. Gauge CD3E / CD3D as before.
##
##   sbatch --mem=48G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##       scripts/34_script_nuclear_size.R [smoke]
## ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
args  <- commandArgs(trailingOnly = TRUE)
SMOKE <- "smoke" %in% args
FLAT  <- "input/CosMx/export1/flatFiles"
REG   <- "output/16_regions/matched/cell_regions.csv.gz"
S31   <- "output/31_radial_cd45"
OUT   <- if (SMOKE) "output/34_nuclear_size_smoke" else "output/34_nuclear_size"
FIG   <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
N_BOOT <- if (SMOKE) 20L else 1000L
UM2_PER_PX2 <- 0.120280945^2
GATE_AUC <- 0.65; GATE_RHO <- 0.2
TIER_BREAKS <- c(-Inf, 9.5, 19.5, 49.5, 99.5, Inf)
TIERS <- c("<10", "10-19", "20-49", "50-99", "100+")
SLIDE <- c(IB1 = 2L, RB1 = 1L)       # c_2_ = IB1 (CLAUDE.md); checked against `samp` below
set.seed(34)

SETS <- list(
  DC        = c("LAMP3", "CCR7", "CCL19", "IDO1", "CD1C", "CLEC9A", "XCR1", "CLEC10A",
                "FCER1A", "CD1A", "CD207", "MARCKSL1"),
  HRS       = c("TNFRSF8", "FUT4", "IL13", "LTA", "CSF2"),
  ambiguous = c("CD274", "PDCD1LG2", "CD40", "BATF3", "IRF4", "CD80", "CD86", "CD83", "JAK2"),
  macro     = c("CD68", "CD163", "C1QA", "CD14", "LYZ"),
  Bprog     = c("CD19", "CD79A", "MS4A1"),
  T         = c("CD3E", "CD3D"))
TEST <- unlist(SETS, use.names = FALSE); GENES <- c("CCL17", "CCL22", TEST)
SET_OF <- setNames(rep(names(SETS), lengths(SETS)), TEST)
METRICS <- c(nuc_area = "Nuclear area (um^2)", nuc_frac = "Nuclear fraction of cell area",
             nuc_ar = "Nuclear aspect ratio", dapi_rel = "Mean.DAPI / FOV median (nuclear paleness, low = pale)")

## ---- house style (canonical copy: 07_script_lightchain_neighborhoods.R:126-147) ----
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
theme_viz <- function(base=11, xgrid=FALSE) {
  th <- theme_minimal(base_size=base) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.major.y=element_line(colour=GRID,linewidth=.25),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.line.x=element_line(colour=AXIS,linewidth=.3), axis.ticks=element_blank(),
          axis.text=element_text(colour=MUTED,size=base-2),
          axis.title=element_text(colour=INK2,size=base-1),
          plot.title=element_text(colour=INK,face="bold",size=base+2),
          plot.subtitle=element_text(colour=INK2,size=base-1),
          plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
          legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left", legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF);cat("wrote",file.path(FIG,n),"\n")}
SAMP_COLS <- c(IB1 = "#2a78d6", RB1 = "#eb6834")
SET_COL <- c(HRS = "#0b3d7a", ambiguous = "#7a5aa6", DC = "#2a9d8f", other = "#9a968d")
OR_BREAKS <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)
pboot <- function(x) { x <- x[is.finite(x)]; n <- length(x)
  min(1, 2 * min((1 + sum(x <= 0)) / (1 + n), (1 + sum(x >= 0)) / (1 + n))) }

## ---- 1. cells: nuclear metrics + panel counts, both samples -----------------------
reg <- fread(cmd = paste("zcat", REG),
             select = c("cell_id", "samp", "fov_id", "is_defect", "tier_matched", "nCount_RNA"))
reg <- reg[tier_matched == TRUE & is_defect == FALSE]
read_sample <- function(s) {
  m <- fread(file.path(FLAT, s, paste0(s, "_metadata_file.csv.gz")),
             select = c("fov", "cell_ID", "Area", "AspectRatio", "NucArea", "NucAspectRatio", "Mean.DAPI"))
  m[, cell_id := sprintf("c_%d_%d_%d", SLIDE[[s]], fov, cell_ID)]
  f <- file.path(FLAT, s, paste0(s, "_exprMat_file.csv.gz"))
  hdr <- names(fread(cmd = sprintf("zcat %s | head -1", f)))
  stopifnot(all(GENES %in% hdr))
  ex <- fread(cmd = sprintf("zcat %s | cut -d, -f%s%s", f, paste(match(c("fov", "cell_ID", GENES), hdr), collapse = ","),
                            if (SMOKE) " | head -300001" else ""))
  ex[, cell_id := sprintf("c_%d_%d_%d", SLIDE[[s]], fov, cell_ID)]
  x <- merge(m[, .(cell_id, Area, AspectRatio, NucArea, NucAspectRatio, dapi = `Mean.DAPI`)], ex[, c("cell_id", GENES), with = FALSE], by = "cell_id")
  x[, samp_file := s]
  x
}
d <- merge(reg, rbindlist(lapply(names(SLIDE), read_sample)), by = "cell_id")
stopifnot(all(d$samp == d$samp_file))
if (!SMOKE) stopifnot(nrow(d) == nrow(reg))
d[, `:=`(cell_area = Area * UM2_PER_PX2, nuc_area = NucArea * UM2_PER_PX2,
         nuc_frac = NucArea / Area, nuc_ar = NucAspectRatio, chemo = CCL17 + CCL22)]
d[, depth0 := nCount_RNA - chemo]
d <- d[depth0 > 0 & is.finite(nuc_area) & nuc_area > 0]
## Mean.DAPI relative to the FOV median: DAPI staining varies by field (RB1 <10 median 3235 vs
## IB1 1772), so paleness is only meaningful against the cell's own field
d[, dapi_rel := dapi / median(dapi), by = .(samp, fov_id)]
d[, tier := factor(cut(chemo, TIER_BREAKS, labels = TIERS), levels = TIERS)]
cat(sprintf("cells: %s\n", paste(d[, .N, by = samp][, sprintf("%s %d", samp, N)], collapse = ", ")))
cat("NucArea missing/zero dropped:", reg[, .N] - nrow(d), "\n")

## ---- 2. GATE -------------------------------------------------------------------------
l16 <- fread(cmd = "zcat output/16_celltypes/cell_celltypes.csv.gz", select = c("cell_id", "samp", "lineage16"))
d <- merge(d, l16[, .(cell_id, lineage16)], by = "cell_id", all.x = TRUE)
auc_of <- function(pos, neg) { r <- rank(c(pos, neg)); np <- as.numeric(length(pos)); nn <- as.numeric(length(neg))
  (sum(r[seq_len(np)]) - np * (np + 1) / 2) / (np * nn) }
gate_a <- d[lineage16 %in% c("T", "B", "myeloid", "endo", "fibro"),
            .(auc_big_vs_lymph = auc_of(nuc_area[lineage16 %in% c("myeloid", "endo", "fibro")],
                                        nuc_area[lineage16 %in% c("T", "B")]),
              n_lymph = sum(lineage16 %in% c("T", "B")), n_big = sum(lineage16 %in% c("myeloid", "endo", "fibro"))), by = samp]
lin_tab <- d[!is.na(lineage16), .(n = .N, nuc_area_med = round(median(nuc_area), 1), cell_area_med = round(median(cell_area), 1),
                                  nuc_frac_med = round(median(nuc_frac), 3)), by = .(samp, lineage16)][order(samp, -nuc_area_med)]
fwrite(lin_tab, file.path(OUT, "gate_nuclear_by_lineage.csv"))
## the gate: vendor NucArea vs the polygon-bounded image nuclear area (same cells)
img <- fread(file.path("output/34_nuclear_size", "nuclei_in_polygon.csv"))[ok == TRUE]
setnames(img, c("nuc_area", "nuc_frac", "nuc_ar"), c("img_nuc_area", "img_nuc_frac", "img_nuc_ar"))
img <- merge(img, d[, .(cell_id, nuc_area, cell_area, nuc_frac, nuc_ar, dapi, fov_id)], by = "cell_id")
## is the vendor's per-cell Mean.DAPI the DNA level of the polygon in the picture?
cat(sprintf("\nvendor Mean.DAPI vs image dna_rel (polygon median DNA, background-referenced): Spearman %.3f (n = %d)\n",
            img[, cor(dapi, dna_rel, method = "spearman")], nrow(img)))
gate_b <- img[has_nucleus == TRUE, .(n = .N, rho_area = cor(nuc_area, img_nuc_area, method = "spearman"),
                  rho_ar = cor(nuc_ar, img_nuc_ar, method = "spearman"),
                  rho_cell_vs_img = cor(cell_area, img_nuc_area, method = "spearman"),
                  rho_polyarea_vs_cellarea = cor(cell_area, poly_area, method = "spearman"))]
GATE_PASS <- gate_b$rho_area >= 0.5
## what the polygon holds, per group: nucleus present? DNA at the centroid? nucleus pressed
## against the outside? lobation and nucleoli where a nucleus exists. FOV bootstrap.
fovs_i <- unique(img$fov_id); idx_i <- split(seq_len(nrow(img)), img$fov_id)
summ <- function(rows) { x <- img[rows]
  rbindlist(lapply(c("candidate", "control"), function(g) { y <- x[group == g]; z <- y[has_nucleus == TRUE]
    data.table(group = g, n = nrow(y), frac_no_nucleus = mean(!y$has_nucleus), dna_centroid_med = median(y$dna_centroid),
               frac_centroid_dim = mean(y$dna_centroid < 1), bright_outside_med = median(y$bright_outside),
               poly_area_med = median(y$poly_area), nuc_area_med = median(z$img_nuc_area), nuc_frac_med = median(y$img_nuc_frac),
               n_maxima_mean = mean(z$n_maxima), n_holes_mean = mean(z$n_holes), nuc_ar_med = median(z$img_nuc_ar),
               multi_comp_frac = mean(z$n_nuc_components > 1),
               dna_rel_med = median(y$dna_rel), frac_dna_rel_lt_0.5 = mean(y$dna_rel < 0.5),
               pale_area_med = median(y$pale_area), pale_frac_med = median(y$pale_frac),
               frac_no_pale = mean(y$pale_area < 8), pale_ar_med = median(y$pale_ar, na.rm = TRUE),
               pale_cv_med = median(y$pale_cv, na.rm = TRUE), pale_holes_mean = mean(y$pale_n_holes, na.rm = TRUE),
               vendor_dapi_med = median(y$dapi)) })) }
S0 <- summ(seq_len(nrow(img)))
SB <- rbindlist(lapply(seq_len(N_BOOT), function(b) summ(unlist(idx_i[as.character(sample(fovs_i, replace = TRUE))], use.names = FALSE))[, b := b]))
stats <- setdiff(names(S0), c("group", "n"))
polysum <- rbindlist(lapply(stats, function(st) {
  dd <- SB[, .(d = get(st)[group == "candidate"] - get(st)[group == "control"]), by = b]$d
  data.table(measure = st, candidate = round(S0[group == "candidate"][[st]], 3), control = round(S0[group == "control"][[st]], 3),
             diff = round(S0[group == "candidate"][[st]] - S0[group == "control"][[st]], 3),
             lo = round(quantile(dd, .025, na.rm = TRUE), 3), hi = round(quantile(dd, .975, na.rm = TRUE), 3)) }))
fwrite(polysum, file.path(OUT, "polygon_nuclei_candidate_vs_control.csv"))
cat("\n== what the vendor polygon holds: candidate vs area-matched control (FOV bootstrap on the difference)\n"); print(polysum)
fwrite(cbind(gate_b[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)],
             data.table(gate_rho = 0.5, gate_pass = GATE_PASS)), file.path(OUT, "gate_verdict.csv"))
fwrite(gate_a[, .(samp, auc_big_vs_lymph = round(auc_big_vs_lymph, 3), n_lymph, n_big)],
       file.path(OUT, "descriptive_lineage_auc.csv"))
cat("\n== GATE: vendor NucArea vs image-measured nuclear area, step-31 cells (unflagged)\n"); print(gate_b)
cat(sprintf("GATE %s\n", if (GATE_PASS) "PASS" else "FAIL -- downstream NOT INTERPRETABLE"))
cat("\n(descriptive) NucArea by step-16 lineage, and myeloid/endo/fibro-over-T/B AUC\n"); print(gate_a); print(lin_tab)

TITLE_PREFIX <- if (GATE_PASS) "" else "GATE FAILED, NOT INTERPRETABLE: "

## ---- 3. descriptives: nuclear size by tier and sample; bright vs dim within cell-area decile
q <- function(x, p) as.numeric(quantile(x, p))
desc <- d[, .(n = .N, nuc_area_med = round(median(nuc_area), 1), nuc_area_q25 = round(q(nuc_area, .25), 1),
              nuc_area_q75 = round(q(nuc_area, .75), 1), cell_area_med = round(median(cell_area), 1),
              nuc_frac_med = round(median(nuc_frac), 3), nuc_ar_med = round(median(nuc_ar), 3),
              dapi_med = round(median(dapi), 1)),
          keyby = .(samp, tier)]
fwrite(desc, file.path(OUT, "nuclear_by_tier.csv"))
cat("\n== nuclear metrics by tier\n"); print(desc)
ib <- d[samp == "IB1"]
ib[, adec := cut(cell_area, unique(quantile(cell_area, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
am <- ib[, .(n_dim = sum(chemo < 10), n_bright = sum(chemo >= 50),
             nuc_area_dim = median(nuc_area[chemo < 10]), nuc_area_bright = median(nuc_area[chemo >= 50]),
             nuc_frac_dim = median(nuc_frac[chemo < 10]), nuc_frac_bright = median(nuc_frac[chemo >= 50]),
             dapi_dim = median(dapi[chemo < 10]), dapi_bright = median(dapi[chemo >= 50]),
             cell_area_med = median(cell_area)), keyby = adec][n_bright >= 10]
am[, `:=`(nuc_area_ratio = round(nuc_area_bright / nuc_area_dim, 3), nuc_frac_diff = round(nuc_frac_bright - nuc_frac_dim, 4),
          dapi_ratio = round(dapi_bright / dapi_dim, 3))]
fwrite(am, file.path(OUT, "bright_vs_dim_within_cell_area_decile.csv"))
cat("\n== IB1 >= 50 vs < 10, within cell-area decile\n"); print(am)

## ---- 4. gene excess by nuclear-metric quintile within brightness band (IB1) ----------
ib[, depth_q := cut(depth0, unique(quantile(depth0, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
ref <- ib[chemo < 10, c(list(dep = sum(as.numeric(depth0))), lapply(.SD, function(z) sum(as.numeric(z)))),
          by = depth_q, .SDcols = TEST]
for (g in TEST) ref[, (paste0("r_", g)) := get(g) / dep]
src <- merge(ib[chemo >= 10], ref[, c("depth_q", paste0("r_", TEST)), with = FALSE], by = "depth_q")
for (g in TEST) src[, (paste0("e_", g)) := depth0 * get(paste0("r_", g))]
src[, band := factor(as.character(tier), levels = TIERS[-1])]
fwrite(src[, c("cell_id", "fov_id", "band", "chemo", "depth0", "cell_area", "nuc_area", "nuc_frac", "nuc_ar", "dapi_rel", TEST), with = FALSE],
       file.path(OUT, "source_cells_nuclear.csv.gz"))
obs_m <- as.matrix(src[, ..TEST]); exp_m <- as.matrix(src[, paste0("e_", TEST), with = FALSE])
fovs <- unique(src$fov_id); idx <- split(seq_len(nrow(src)), src$fov_id)

run_metric <- function(metric) {
  src[, rq := cut(get(metric), unique(quantile(get(metric), 0:5 / 5)), include.lowest = TRUE, labels = FALSE), by = band]
  key <- paste(src$band, src$rq)
  lr_of <- function(rows) {                     # (band x quintile) x gene
    k <- key[rows]
    t(sapply(split(rows, k), function(i) log((colSums(obs_m[i, , drop = FALSE]) + 0.5) / (colSums(exp_m[i, , drop = FALSE]) + 0.5))))
  }
  L0 <- lr_of(seq_len(nrow(src)))
  BB <- lapply(seq_len(N_BOOT), function(b) lr_of(unlist(idx[as.character(sample(fovs, replace = TRUE))], use.names = FALSE)))
  keys <- rownames(L0)
  arr <- array(NA_real_, c(length(BB), length(keys), length(TEST)), dimnames = list(NULL, keys, TEST))
  for (b in seq_along(BB)) arr[b, rownames(BB[[b]]), ] <- BB[[b]]
  ex <- rbindlist(lapply(keys, function(k) data.table(metric = metric, band = sub(" .*", "", k), rq = as.integer(sub(".* ", "", k)),
    n = sum(key == k), gene = TEST, set = SET_OF[TEST], excess = exp(L0[k, ]),
    lo = exp(apply(arr[, k, ], 2, quantile, .025, na.rm = TRUE)), hi = exp(apply(arr[, k, ], 2, quantile, .975, na.rm = TRUE)))))
  tvb <- rbindlist(lapply(levels(src$band), function(bd) {
    kt <- paste(bd, 5); kb <- paste(bd, 1)
    if (!all(c(kt, kb) %in% keys)) return(NULL)
    d0 <- L0[kt, ] - L0[kb, ]; db <- arr[, kt, ] - arr[, kb, ]
    g <- rowMeans(db[, c("CD3E", "CD3D")]); g0 <- mean(d0[c("CD3E", "CD3D")])
    data.table(metric = metric, band = bd, gene = TEST, set = SET_OF[TEST], ratio = exp(d0),
               lo = exp(apply(db, 2, quantile, .025, na.rm = TRUE)), hi = exp(apply(db, 2, quantile, .975, na.rm = TRUE)),
               se = apply(db, 2, sd, na.rm = TRUE), p_boot = apply(db, 2, pboot),
               net_vs_gauge = exp(d0 - g0), net_lo = exp(apply(db - g, 2, quantile, .025, na.rm = TRUE)),
               net_hi = exp(apply(db - g, 2, quantile, .975, na.rm = TRUE)))
  }))
  tvb[, p_wald := 2 * pnorm(-abs(log(ratio) / se))]
  tvb[, q := p.adjust(p_wald, "BH")]
  list(ex = ex, tvb = tvb, med = src[, .(metric = metric, med = median(get(metric))), by = .(band, rq)])
}
R <- lapply(names(METRICS), run_metric)
ex  <- rbindlist(lapply(R, `[[`, "ex")); tvb <- rbindlist(lapply(R, `[[`, "tvb")); med <- rbindlist(lapply(R, `[[`, "med"))
ex <- merge(ex, med, by = c("metric", "band", "rq"))
num <- c("excess", "lo", "hi"); ex[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(ex[order(metric, band, gene, rq)], file.path(OUT, "excess_by_quintile.csv"))
numt <- c("ratio", "lo", "hi", "net_vs_gauge", "net_lo", "net_hi"); tvb[, (numt) := lapply(.SD, round, 3), .SDcols = numt]
fwrite(tvb[order(metric, band, set, gene)], file.path(OUT, "top_vs_bottom_quintile.csv"))
cat("\n== top vs bottom quintile of NUCLEAR AREA within band (ratio of excess), HRS / gauge / DC genes\n")
print(tvb[metric == "nuc_area" & gene %in% c("IL13", "TNFRSF8", "LTA", "BATF3", "IRF4", "CD3E", "CD3D", "LAMP3", "XCR1", "CCR7"),
          .(band, set, gene, ratio, lo, hi, q = signif(q, 2), net = net_vs_gauge)][order(band, set, gene)])
cat("\nq < 0.05 across all metrics:\n"); print(tvb[q < 0.05, .(metric, band, set, gene, ratio, lo, hi, q = signif(q, 2), net = net_vs_gauge)][order(metric, band, -ratio)])

## ---- 5. figures ----------------------------------------------------------------------
dd <- d[, if (.N > 20000) .SD[sample(.N, 20000)] else .SD, by = .(samp, tier)]
dl <- melt(dd[, .(samp, tier, nuc_area, nuc_frac, nuc_ar, dapi_rel)], id.vars = c("samp", "tier"), variable.name = "metric")
dl[, metric_lab := factor(METRICS[as.character(metric)], levels = METRICS)]
f1 <- ggplot(dl, aes(tier, value, fill = samp)) +
  geom_violin(colour = NA, scale = "width", width = .8, position = position_dodge(width = .85)) +
  geom_boxplot(aes(group = interaction(tier, samp)), width = .15, outlier.shape = NA, fill = NA, colour = INK2, linewidth = .3,
               position = position_dodge(width = .85)) +
  facet_wrap(~ metric_lab, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = SAMP_COLS, name = NULL) +
  labs(title = paste0(TITLE_PREFIX, "Nuclear size and shape by chemokine brightness"),
       subtitle = "Vendor DAPI segmentation (NucArea, NucAspectRatio); defect-free matched-tier cells; the < 10 tier is a random 20,000 per sample for display",
       x = "CCL17+CCL22 counts per cell", y = NULL,
       caption = sprintf("RB1 has 9 cells >= 50 and none >= 100. Gate: NucArea myeloid/endo/fibro over T/B AUC %s; Spearman(NucArea, step-31 DNA at 5-10 um) = %.2f. Source: output/34_nuclear_size/nuclear_by_tier.csv.",
                         paste(sprintf("%s %.2f", gate_a$samp, gate_a$auc_big_vs_lymph), collapse = " / "), gate_b$rho_nuc)) +
  theme_viz(10) + theme(strip.text = element_text(colour = INK2, face = "bold", hjust = 0))
sf("fig34_nuclear_by_tier.png", f1, 16, 4.8)
cat("\n== top vs bottom quintile of DAPI_REL within band (quintile 5 = BRIGHTEST DAPI; a pale-nucleus signal reads as ratio < 1 for HRS genes)\n")
print(tvb[metric == "dapi_rel" & gene %in% c("IL13", "TNFRSF8", "LTA", "BATF3", "IRF4", "CD3E", "CD3D", "LAMP3", "XCR1", "CCR7", "LYZ"),
          .(band, set, gene, ratio, lo, hi, q = signif(q, 2), net = net_vs_gauge)][order(band, set, gene)])

e2 <- ex[metric == "nuc_area"]
e2[, band := factor(band, levels = TIERS[-1])]
e2[, panel := factor(fcase(set == "HRS", "HRS-leaning", set == "ambiguous", "HRS or activated APC",
                           set == "DC", "Dendritic-cell lineage", default = "Macrophage / B / T (gauge)"),
                     levels = c("HRS-leaning", "HRS or activated APC", "Dendritic-cell lineage", "Macrophage / B / T (gauge)"))]
e2[, col := fifelse(set %in% names(SET_COL), set, "other")]
SHOW <- c("IL13", "TNFRSF8", "LTA", "FUT4", "BATF3", "IRF4", "CD80", "CD274", "LAMP3", "XCR1", "CCR7", "CD1C", "LYZ", "CD79A", "CD3E", "CD3D")
e2[, shown := gene %in% SHOW]
f2 <- ggplot(e2, aes(rq, excess, group = gene)) +
  geom_hline(yintercept = 1, colour = AXIS, linewidth = .4) +
  geom_line(data = e2[shown == FALSE], colour = "#d8d6cf", linewidth = .3) +
  geom_line(data = e2[shown == TRUE], aes(colour = col), linewidth = .7) +
  geom_point(data = e2[shown == TRUE], aes(colour = col), size = 1.1) +
  ggrepel::geom_text_repel(data = e2[shown == TRUE & rq == 5], aes(label = gene, colour = col), size = 2.4, nudge_x = .35,
                           direction = "y", hjust = 0, segment.size = .15, min.segment.length = 0, seed = 1) +
  facet_grid(panel ~ band, scales = "free_y") +
  scale_y_continuous(trans = "log2", breaks = OR_BREAKS, labels = function(x) paste0(x, "x")) +
  scale_x_continuous(breaks = 1:5, limits = c(1, 6.4)) +
  scale_colour_manual(values = SET_COL, guide = "none") +
  labs(title = paste0(TITLE_PREFIX, "Gene excess across nuclear-area quintiles, within each chemokine brightness band"),
       subtitle = "IB1 sources. Excess = observed / expected counts, expected from depth-matched cells with < 10 chemokine counts. Quintile 5 = largest nuclei",
       x = "Nuclear-area quintile within band", y = "Excess (log scale)",
       caption = "A rise in the HRS-leaning row counts only where it exceeds the CD3E/CD3D gauge. Source: output/34_nuclear_size/excess_by_quintile.csv.") +
  theme_viz(9) + theme(strip.text = element_text(colour = INK2, face = "bold"))
sf("fig34_excess_by_nuclear_area.png", f2, 13, 9)

t3 <- tvb[metric == "nuc_area"]; t3[, band := factor(band, levels = TIERS[-1])]
t3[, gene := factor(gene, levels = rev(TEST))]
GCOL2 <- c(HRS = "#0b3d7a", ambiguous = "#7a5aa6", DC = "#2a9d8f", macro = "#9a968d", Bprog = "#9a968d", T = "#d03b3b")
f3 <- ggplot(t3, aes(ratio, gene, colour = set)) +
  geom_vline(xintercept = 1, colour = AXIS, linewidth = .4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = .5, orientation = "y") +
  geom_point(size = 1.8) +
  facet_wrap(~ band, nrow = 1) +
  scale_x_continuous(trans = "log2", breaks = OR_BREAKS, labels = function(x) paste0(x, "x")) +
  scale_colour_manual(values = GCOL2, name = NULL) +
  labs(title = paste0(TITLE_PREFIX, "Largest vs smallest nuclear-area quintile, within chemokine band"),
       subtitle = "Ratio of excess (top / bottom quintile), 95% FOV cluster bootstrap; T genes (red) are the leakage gauge",
       x = "Ratio (log scale)", y = NULL) + theme_viz(9, xgrid = TRUE)
sf("fig34_top_vs_bottom.png", f3, 13, 6)
cat("\ndone:", OUT, "\n")
