## ============================================================================
## Step 30 -- "CD45-dim cell among CD45-bright neighbours" as an HRS score, and whether
## HRS-leaning genes rise with it WITHIN chemokine brightness bands (IB1).
##
## WHY. Classic HRS cells are CD45-negative and typically rosetted by CD45+ small
## lymphocytes; that pattern is visible in the CD45 morphology channel and does not
## depend on RNA. Step 28 found the HRS-leaning genes climb with chemokine brightness; if
## the rosette pattern marks HRS cells, the same genes should ALSO climb with it inside a
## brightness band, i.e. on an axis independent of how the cells were selected.
##
## SCORE (fixed before running)
##   own   = the cell's vendor Mean.CD45 (mean stain intensity inside its segmentation)
##   nbr   = 75th percentile of Mean.CD45 over defect-free cells within NBR_UM centre-to-
##           centre (>= MIN_NBR of them; else NA). The 75th percentile catches a ring of
##           bright lymphocytes without needing every neighbour to be bright.
##   raw   = log2((nbr + 1) / (own + 1))
##   score = raw - median(raw) in the cell's Area decile (all defect-free IB1 cells).
##           CD45 is membranous, so a LARGE cell has a lower MEAN intensity from geometry
##           alone, and bright chemokine cells are large: unadjusted, the score would
##           rediscover cell size. As a local ratio it also cancels field-level staining
##           variation (IB1 has region-dependent staining, CLAUDE.md fact 2/9).
##
## TEST, within each chemokine band (10-19, 20-49, 50-99, 100+): quintiles of score
##   within band; per gene, pooled excess = obs / exp, exp = depth0 x the gene's rate in
##   < 10 cells of the same depth0 decile (step 28's reference; depth0 excludes the
##   chemokine counts). Primary contrast = top vs bottom quintile; FOV cluster bootstrap,
##   percentile p with the (1+count)/(1+n) correction plus Wald p, BH over genes x bands.
##
## CONTROLS
##   CD3E, CD3D  LEAKAGE GAUGES. A rosette IS T cells, and their transcripts leak into the
##               cell they surround (fact 6); IL13, LTA and TNFRSF8 are also made by
##               activated T cells. An HRS-leaning rise is read as evidence only where it
##               clearly exceeds the CD3E/CD3D rise in the same band.
##   PTPRC       the CD45 gene. Should FALL with the score (protein-dim cell) unless the
##               neighbours' PTPRC leaks in. Also: sample-wide, Mean.CD45 must track
##               PTPRC at all, or the stain is not measuring CD45.
##   LAMP3, XCR1, CCR7   DC-side genes, for direction.
##
##   sbatch --mem=48G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##       scripts/30_script_cd45_rosette.R [smoke]
## ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(RANN) })

args  <- commandArgs(trailingOnly = TRUE)
SMOKE <- "smoke" %in% args
FLAT  <- "input/CosMx/export1/flatFiles/IB1"
REG   <- "output/16_regions/matched/cell_regions.csv.gz"
OUT   <- if (SMOKE) "output/30_cd45_rosette_smoke" else "output/30_cd45_rosette"
FIG   <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
N_BOOT <- if (SMOKE) 20L else 1000L
if (SMOKE) options(warn = 1)        # print warnings as they happen, so a smoke run shows them
NBR_UM <- 30; MAX_NBR <- 12L; MIN_NBR <- 3L; NBR_Q <- 0.75
set.seed(30)

HRS  <- c("IL13", "BATF3", "TNFRSF8", "IRF4", "LTA")
DC   <- c("LAMP3", "XCR1", "CCR7")
GAUGE <- c("CD3E", "CD3D", "PTPRC")
TEST <- c(HRS, DC, GAUGE)
GENES <- c("CCL17", "CCL22", TEST)
GROUP <- c(setNames(rep("HRS-leaning", length(HRS)), HRS), setNames(rep("Dendritic cell", length(DC)), DC),
           CD3E = "T leakage gauge", CD3D = "T leakage gauge", PTPRC = "CD45 gene (PTPRC)")
BANDS <- c(9.5, 19.5, 49.5, 99.5, Inf); BAND_LAB <- c("10-19", "20-49", "50-99", "100+")

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
GCOL <- c("HRS-leaning" = "#0b3d7a", "Dendritic cell" = "#2a9d8f",
          "T leakage gauge" = "#9a968d", "CD45 gene (PTPRC)" = "#d03b3b")
LOG_BREAKS <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)   # explicit (CLAUDE.md, step 25)

## ---- 1. cells: IB1, matched tier, defect-free; CD45 intensity; panel counts ---------
reg <- fread(cmd = paste("zcat", REG),
             select = c("cell_id", "samp", "fov_id", "x_slide_mm", "y_slide_mm", "is_defect",
                        "tier_matched", "nCount_RNA", "Area.um2"))
reg <- reg[samp == "IB1" & tier_matched == TRUE & is_defect == FALSE]
## c_2_ is IB1 (CLAUDE.md); ids are BUILT from IB1's own files and then joined on cell_id
## to rows already filtered on samp == "IB1", so a wrong prefix would fail the join.
meta <- fread(file.path(FLAT, "IB1_metadata_file.csv.gz"),
              select = c("fov", "cell_ID", "Mean.CD45", "Mean.Membrane"))
meta[, cell_id := sprintf("c_2_%d_%d", fov, cell_ID)]
hdr <- names(fread(cmd = sprintf("zcat %s/IB1_exprMat_file.csv.gz | head -1", FLAT)))
stopifnot(all(GENES %in% hdr))
ex <- fread(cmd = sprintf("zcat %s/IB1_exprMat_file.csv.gz | cut -d, -f%s%s", FLAT,
                          paste(match(c("fov", "cell_ID", GENES), hdr), collapse = ","),
                          if (SMOKE) " | head -300001" else ""))
ex[, cell_id := sprintf("c_2_%d_%d", fov, cell_ID)]
d <- merge(reg, meta[, .(cell_id, cd45 = `Mean.CD45`, membrane = `Mean.Membrane`)], by = "cell_id")
d <- merge(d, ex[, c("cell_id", GENES), with = FALSE], by = "cell_id")
if (!SMOKE) stopifnot(nrow(d) == nrow(reg))
cat(sprintf("IB1 defect-free matched cells with intensity + counts: %d\n", nrow(d)))
d[, chemo := CCL17 + CCL22]
d[, depth0 := nCount_RNA - chemo]
d <- d[depth0 > 0]

## ---- 2. the rosette score -------------------------------------------------------------
xy <- as.matrix(d[, .(x_slide_mm * 1000, y_slide_mm * 1000)])
nn <- nn2(xy, xy, k = MAX_NBR + 1L, searchtype = "radius", radius = NBR_UM)
idx <- nn$nn.idx[, -1, drop = FALSE]                 # drop self (first column)
idx[idx == 0L] <- NA_integer_
cdv <- d$cd45
nbr_vals <- matrix(cdv[idx], nrow = nrow(idx))
n_nbr <- rowSums(!is.na(nbr_vals))
nbr_q <- apply(nbr_vals, 1, function(v) if (sum(!is.na(v)) >= MIN_NBR) quantile(v, NBR_Q, na.rm = TRUE) else NA_real_)
d[, `:=`(n_nbr = n_nbr, cd45_nbr = nbr_q, raw = log2((nbr_q + 1) / (cd45 + 1)))]
d[, area_q := cut(Area.um2, unique(quantile(Area.um2, 0:10 / 10, na.rm = TRUE)),
                  include.lowest = TRUE, labels = FALSE)]
d[, rosette := raw - median(raw, na.rm = TRUE), by = area_q]
cat(sprintf("cells with >= %d neighbours within %g um: %d of %d\n", MIN_NBR, NBR_UM,
            d[!is.na(rosette), .N], nrow(d)))

## ---- 2b. GATE: does per-cell Mean.CD45 measure the cell's own CD45 at all? ----------
## (i) against the CD45 gene, and (ii) the stronger test: CD45-NEGATIVE lineages
## (step-16 argmax fibro, endo) against CD45-POSITIVE ones (T, B). A per-cell mean
## intensity can be dominated by a bright neighbour's membrane spilling over the
## boundary, which would make the rosette score meaningless. Fixed rule: lineage AUC
## < GATE_AUC means the score is reported as NOT INTERPRETABLE, whatever the genes do.
GATE_AUC <- 0.70
d[, ptprc_cp10k := 1e4 * PTPRC / depth0]
rho_stain <- d[, cor(cd45, ptprc_cp10k, method = "spearman")]
cat(sprintf("Mean.CD45 vs PTPRC cp10k: Spearman %.3f\n", rho_stain))
l16 <- fread(cmd = "zcat output/16_celltypes/cell_celltypes.csv.gz",
             select = c("cell_id", "samp", "lineage16"))[samp == "IB1"]
d <- merge(d, l16[, .(cell_id, lineage16)], by = "cell_id", all.x = TRUE)
auc_of <- function(pos, neg) {             # P(a CD45+ lineage cell is brighter), Mann-Whitney
  r <- rank(c(pos, neg)); np <- as.numeric(length(pos)); nn <- as.numeric(length(neg))  # no int overflow
  (sum(r[seq_len(np)]) - np * (np + 1) / 2) / (np * nn)
}
lin <- d[lineage16 %in% c("T", "B", "fibro", "endo")]
gate_auc <- auc_of(lin[lineage16 %in% c("T", "B")]$cd45, lin[lineage16 %in% c("fibro", "endo")]$cd45)
gate_auc_rosette <- auc_of(lin[lineage16 %in% c("fibro", "endo") & !is.na(rosette)]$rosette,
                           lin[lineage16 %in% c("T", "B") & !is.na(rosette)]$rosette)
lin_tab <- d[!is.na(lineage16), .(n = .N, cd45_med = median(cd45), rosette_med = median(rosette, na.rm = TRUE)),
             by = lineage16][order(-cd45_med)]
fwrite(lin_tab, file.path(OUT, "gate_cd45_by_lineage.csv"))
GATE_PASS <- gate_auc >= GATE_AUC
## Diagnostics, so a FAIL can be read: is it this slide's staining (then RB1 passes) or
## the per-cell measurement in general (then RB1 fails too)? And does the vendor's other
## per-cell summary, Max.CD45, do better than the mean? Neither changes the gate.
diag_auc <- function(s) {
  m <- fread(file.path("input/CosMx/export1/flatFiles", s, paste0(s, "_metadata_file.csv.gz")),
             select = c("fov", "cell_ID", "Mean.CD45", "Max.CD45"))
  m[, cell_id := sprintf("c_%d_%d_%d", c(IB1 = 2L, RB1 = 1L)[[s]], fov, cell_ID)]
  r <- fread(cmd = paste("zcat", REG), select = c("cell_id", "samp", "is_defect", "tier_matched"))[
    samp == s & tier_matched == TRUE & is_defect == FALSE]
  l <- fread(cmd = "zcat output/16_celltypes/cell_celltypes.csv.gz",
             select = c("cell_id", "samp", "lineage16"))[samp == s & lineage16 %in% c("T", "B", "fibro", "endo")]
  x <- merge(merge(r, l[, .(cell_id, lineage16)], by = "cell_id"), m, by = "cell_id")
  pos <- x$lineage16 %in% c("T", "B")
  data.table(samp = s, n_TB = sum(pos), n_fibro_endo = sum(!pos),
             auc_mean = auc_of(x$`Mean.CD45`[pos], x$`Mean.CD45`[!pos]),
             auc_max = auc_of(x$`Max.CD45`[pos], x$`Max.CD45`[!pos]))
}
diag <- rbindlist(lapply(c("IB1", "RB1"), diag_auc))
fwrite(diag, file.path(OUT, "gate_diagnostics.csv"))
cat("\n== gate diagnostics: AUC of CD45 intensity, T/B over fibro/endo, by sample and summary\n")
print(diag)
fwrite(data.table(rho_ptprc = rho_stain, auc_cd45_TB_vs_fibro_endo = gate_auc,
                  auc_rosette_fibro_endo_vs_TB = gate_auc_rosette, gate_auc = GATE_AUC,
                  gate_pass = GATE_PASS), file.path(OUT, "gate_verdict.csv"))
TITLE_PREFIX <- if (GATE_PASS) "" else "GATE FAILED, NOT INTERPRETABLE: "
cat("\n== GATE: Mean.CD45 by step-16 lineage\n"); print(lin_tab)
cat(sprintf("AUC, Mean.CD45 higher in T/B than fibro/endo: %.3f (gate >= %.2f) -> %s\n",
            gate_auc, GATE_AUC, if (GATE_PASS) "PASS" else "FAIL: rosette score NOT interpretable"))
cat(sprintf("AUC, rosette score higher in fibro/endo than T/B: %.3f\n", gate_auc_rosette))

## ---- 3. expected counts from the < 10 reference, depth0 deciles ----------------------
d[, depth_q := cut(depth0, unique(quantile(depth0, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
ref <- d[chemo < 10, c(list(dep = sum(as.numeric(depth0))), lapply(.SD, function(z) sum(as.numeric(z)))),
         by = depth_q, .SDcols = TEST]
for (g in TEST) ref[, (paste0("r_", g)) := get(g) / dep]
d <- merge(d, ref[, c("depth_q", paste0("r_", TEST)), with = FALSE], by = "depth_q")
for (g in TEST) d[, (paste0("e_", g)) := depth0 * get(paste0("r_", g))]
d[, (paste0("r_", TEST)) := NULL]

src <- d[chemo >= 10 & !is.na(rosette)]
src[, band := factor(cut(chemo, BANDS, labels = BAND_LAB), levels = BAND_LAB)]
src[, rq := cut(rosette, unique(quantile(rosette, 0:5 / 5)), include.lowest = TRUE,
                labels = FALSE), by = band]
src[, hrs_score := rowMeans(sapply(c("IL13", "BATF3", "TNFRSF8"), function(g)
  log2((get(g) + 0.5) / (get(paste0("e_", g)) + 0.5))))]
fwrite(src[, c("cell_id", "fov_id", "band", "chemo", "Area.um2", "cd45", "cd45_nbr", "n_nbr",
               "raw", "rosette", "rq", "hrs_score", TEST), with = FALSE][order(band, -rosette)],
       file.path(OUT, "rosette_source_cells.csv.gz"))

## ---- 4. excess by band x quintile, top vs bottom, FOV bootstrap ----------------------
agg <- src[, lapply(.SD, sum), by = .(fov_id, band, rq), .SDcols = c(TEST, paste0("e_", TEST))]
excess_of <- function(A) {
  S <- A[, lapply(.SD, function(z) sum(z * w)), by = .(band, rq),
         .SDcols = c(TEST, paste0("e_", TEST))]
  rbindlist(lapply(TEST, function(g)
    S[, .(band, rq, gene = g, obs = get(g), exp = get(paste0("e_", g)),
          lr = log((get(g) + 0.5) / (get(paste0("e_", g)) + 0.5)))]))
}
pboot <- function(x) { x <- x[is.finite(x)]; n <- length(x)
  min(1, 2 * min((1 + sum(x <= 0)) / (1 + n), (1 + sum(x >= 0)) / (1 + n))) }
E0 <- excess_of(agg[, w := 1])
fovs <- unique(agg$fov_id)
BB <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
  w <- data.table(fov_id = sample(fovs, replace = TRUE))[, .(w = .N), by = fov_id]
  excess_of(merge(agg[, !"w"], w, by = "fov_id"))[, .(band, rq, gene, lr, b = b)]
}))
ci <- BB[, .(lo = quantile(lr, .025, na.rm = TRUE), hi = quantile(lr, .975, na.rm = TRUE)),
         by = .(band, rq, gene)]
ex <- merge(E0, ci, by = c("band", "rq", "gene"))
ex[, `:=`(excess = exp(lr), lo = exp(lo), hi = exp(hi), group = GROUP[gene])]
nq <- src[, .(n_cells = .N, rosette_med = median(rosette)), by = .(band, rq)]
ex <- merge(ex, nq, by = c("band", "rq"))
fwrite(ex[order(band, gene, rq), .(band, rq, n_cells, rosette_med = round(rosette_med, 3), group, gene,
                                   obs, exp = round(exp, 2), excess = round(excess, 3),
                                   lo = round(lo, 3), hi = round(hi, 3))],
       file.path(OUT, "excess_by_rosette_quintile.csv"))

topq <- function(T) T[, .(top = lr[rq == max(rq)], bot = lr[rq == min(rq)]), by = .(band, gene)][, d := top - bot]
D0 <- topq(E0)
DB <- BB[, .(top = lr[rq == max(rq)], bot = lr[rq == min(rq)]), by = .(b, band, gene)][, d := top - bot]
tvb <- DB[, .(lo = quantile(d, .025, na.rm = TRUE), hi = quantile(d, .975, na.rm = TRUE),
              se = sd(d, na.rm = TRUE), p_boot = pboot(d)), by = .(band, gene)]
tvb <- merge(D0[, .(band, gene, d)], tvb, by = c("band", "gene"))
tvb[, `:=`(ratio = exp(d), lo = exp(lo), hi = exp(hi), p_wald = 2 * pnorm(-abs(d / se)),
           group = GROUP[gene])]
tvb[, q := p.adjust(p_wald, "BH")]
## the leakage-gauge comparison: each HRS gene's top/bottom ratio over the band's mean
## CD3E/CD3D ratio, bootstrapped jointly
gauge_b <- DB[gene %in% c("CD3E", "CD3D"), .(g = mean(d)), by = .(b, band)]
net <- merge(DB[gene %in% c(HRS, DC)], gauge_b, by = c("b", "band"))[, net := d - g]
net0 <- merge(D0[gene %in% c(HRS, DC)], D0[gene %in% c("CD3E", "CD3D"), .(g = mean(d)), by = band],
              by = "band")[, net := d - g]
netci <- net[, .(net_lo = quantile(net, .025, na.rm = TRUE), net_hi = quantile(net, .975, na.rm = TRUE),
                 net_p = pboot(net)), by = .(band, gene)]
tvb <- merge(tvb, merge(net0[, .(band, gene, net)], netci, by = c("band", "gene")),
             by = c("band", "gene"), all.x = TRUE)
tvb[, `:=`(net_ratio = exp(net), net_lo = exp(net_lo), net_hi = exp(net_hi))]
num <- c("ratio", "lo", "hi", "net_ratio", "net_lo", "net_hi")
fwrite(tvb[order(band, group, gene), c("band", "group", "gene", num, "p_boot", "p_wald", "q", "net_p"), with = FALSE][
  , (num) := lapply(.SD, round, 3), .SDcols = num],
  file.path(OUT, "top_vs_bottom_quintile.csv"))

## within-band rank correlation of the per-cell HRS-leaning score with the rosette score
rho <- src[, .(rho = cor(rosette, hrs_score, method = "spearman"), n = .N), by = band]
rhob <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
  pick <- sample(fovs, replace = TRUE)
  ## a FOV drawn twice must contribute its cells twice: duplication is the point
  s <- src[data.table(fov_id = pick), on = "fov_id", nomatch = 0, allow.cartesian = TRUE]
  s[, .(rho = cor(rosette, hrs_score, method = "spearman")), by = band][, b := b]
}))
rho <- merge(rho, rhob[, .(lo = quantile(rho, .025, na.rm = TRUE), hi = quantile(rho, .975, na.rm = TRUE),
                           p_boot = pboot(rho)), by = band], by = "band")
fwrite(rho, file.path(OUT, "rho_hrs_score_vs_rosette.csv"))

## is the rosette score itself higher in brighter bands?
bandsum <- d[!is.na(rosette)][, band := fifelse(chemo < 10, "<10", as.character(cut(chemo, BANDS, labels = BAND_LAB)))][
  , .(n = .N, rosette_med = median(rosette), cd45_med = median(cd45), area_med = median(Area.um2),
      frac_top_decile = mean(rosette >= quantile(d$rosette, .9, na.rm = TRUE))), by = band]
fwrite(bandsum, file.path(OUT, "rosette_by_band.csv"))

cat("\n== rosette score by chemokine band\n"); print(bandsum)
cat("\n== top vs bottom rosette quintile, within band (ratio of excess); net = over the CD3E/CD3D gauge\n")
print(tvb[order(band, group, gene), .(band, group, gene, ratio = round(ratio, 2), lo = round(lo, 2),
                                       hi = round(hi, 2), q = signif(q, 2), net = round(net_ratio, 2),
                                       net_lo = round(net_lo, 2), net_hi = round(net_hi, 2))])
cat("\n== Spearman(rosette, per-cell HRS-leaning score) within band\n"); print(rho)

## ---- 5. figures --------------------------------------------------------------------------
bb <- d[!is.na(rosette)][, band := factor(fifelse(chemo < 10, "<10", as.character(cut(chemo, BANDS, labels = BAND_LAB))),
                                           levels = c("<10", BAND_LAB))]
bbs <- rbind(bb[band == "<10"][sample(.N, min(.N, 20000L))], bb[band != "<10"])
f1 <- ggplot(bbs, aes(band, rosette)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_violin(fill = "#dfe6ef", colour = NA, scale = "width", width = .8) +
  geom_boxplot(width = .14, outlier.shape = NA, colour = "#0b3d7a", fill = NA, linewidth = .35) +
  coord_cartesian(ylim = quantile(bbs$rosette, c(.01, .99))) +
  labs(title = paste0(TITLE_PREFIX, "CD45 rosette score by chemokine brightness"),
       subtitle = sprintf("Score = log2(neighbour CD45 / own CD45), neighbour = 75th pct of Mean.CD45 within %g um, minus the median for the cell's area decile", NBR_UM),
       x = "CCL17+CCL22 counts per cell", y = "CD45 rosette score (log2)",
       caption = sprintf("IB1 defect-free matched cells; the <10 band is a random 20,000 for display. Higher = dimmer than its neighbours. Mean.CD45 vs PTPRC cp10k, all cells: Spearman %.2f.", rho_stain)) +
  theme_viz(10)
sf("fig30_rosette_by_band.png", f1, 8, 4.6)

ex[, band := factor(band, levels = BAND_LAB)]
ex[, group := factor(group, levels = names(GCOL))]
f2 <- ggplot(ex, aes(rq, excess, colour = group, group = gene)) +
  geom_hline(yintercept = 1, colour = AXIS, linewidth = .4) +
  geom_line(linewidth = .6) + geom_point(size = 1.2) +
  ggrepel::geom_text_repel(data = ex[rq == max(rq)], aes(label = gene), size = 2.5, nudge_x = .4,
                           direction = "y", hjust = 0, segment.size = .15, min.segment.length = 0, seed = 1,
                           show.legend = FALSE) +
  facet_grid(group ~ band, scales = "free_y") +
  scale_y_continuous(trans = "log2", breaks = LOG_BREAKS, labels = function(x) paste0(x, "x")) +
  scale_x_continuous(breaks = 1:5, limits = c(1, 6.2)) +
  scale_colour_manual(values = GCOL, guide = "none") +
  labs(title = paste0(TITLE_PREFIX, "Gene excess across CD45 rosette quintiles, within each chemokine brightness band"),
       subtitle = "Excess = observed / expected counts, expected from depth-matched cells with < 10 chemokine counts. Quintile 5 = most CD45-dim against CD45-bright neighbours",
       x = "CD45 rosette quintile within band", y = "Excess (log scale)",
       caption = "A rise in the HRS-leaning row counts only where it exceeds the T leakage gauge row: rosetting T cells leak transcripts into the cell they surround. Source: output/30_cd45_rosette/excess_by_rosette_quintile.csv.") +
  theme_viz(9) + theme(strip.text = element_text(colour = INK2, face = "bold"))
sf("fig30_excess_by_rosette.png", f2, 13, 9)

tv <- copy(tvb)[, band := factor(band, levels = BAND_LAB)]
tv[, gene := factor(gene, levels = rev(TEST))]
f3 <- ggplot(tv, aes(ratio, gene, colour = GROUP[as.character(gene)])) +
  geom_vline(xintercept = 1, colour = AXIS, linewidth = .4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = .5, orientation = "y") +
  geom_point(size = 2) +
  facet_wrap(~ band, nrow = 1) +
  scale_x_continuous(trans = "log2", breaks = LOG_BREAKS, labels = function(x) paste0(x, "x")) +
  scale_colour_manual(values = GCOL, name = NULL) +
  labs(title = paste0(TITLE_PREFIX, "Top vs bottom CD45 rosette quintile, within chemokine band"),
       subtitle = "Ratio of excess (top / bottom quintile), 95% FOV cluster bootstrap",
       x = "Ratio (log scale)", y = NULL,
       caption = "Compare HRS-leaning genes against the CD3E/CD3D leakage gauges in the same band; net ratios are in output/30_cd45_rosette/top_vs_bottom_quintile.csv.") +
  theme_viz(9, xgrid = TRUE)
sf("fig30_top_vs_bottom.png", f3, 12, 4.8)

cat("\ndone:", OUT, "\n")
