## ============================================================================
## Step 29 -- Bright chemokine cells: the per-cell list, figures for step 28's two
## observations, and the ROI candidates that scripts/29_bright_cell_rois.py renders.
##
## WHAT IT WRITES (output/29_bright_cells/)
##   source_cells.csv.gz     every source cell (CCL17+CCL22 >= 10) in both samples:
##                           canvas px, tier, counts, the HRS-leaning score, k-of-3
##   roi_candidates.csv      ranked candidates per ROI group (IB1); the renderer takes
##                           the first 20 per group that sit on imaged tissue
##   ib1_chemokine_tx.csv.gz every CCL17/CCL22 molecule in IB1, canvas px (contours)
##   figures/                plotting only, from step 28's saved tables + the list
##
## HRS-LEANING SCORE (fixed before looking at any ROI)
##   per gene g in IL13, BATF3, TNFRSF8: exp_g = depth0 x rate_g, rate_g from the
##   same sample's < 10 cells in the same depth0 decile (step 28's reference; depth0
##   excludes the chemokine counts). score = mean_g log2((obs_g + 0.5) / (exp_g + 0.5)).
##   It RANKS candidates; it is not a call. Step 28 showed a k-of-3 gate mostly selects
##   depth, and that among >= 50 cells co-occurrence is barely above independence.
##
## ROI GROUPS (IB1): tiers high = 100+, medium = 50-99, low = 10-49 counts, each as
##   *_random (seeded draw) and *_hrs (descending score). Cells in a group are >= 150 um
##   apart, so no two of its windows overlap; *_random excludes that tier's *_hrs picks.
##
##   sbatch --mem=48G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##       scripts/29_script_bright_cells.R [smoke]
## ============================================================================
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

args  <- commandArgs(trailingOnly = TRUE)
SMOKE <- "smoke" %in% args
FLAT  <- "input/CosMx/export1/flatFiles"
REG   <- "output/16_regions/matched/cell_regions.csv.gz"
S28   <- "output/28_bright_sources"
OUT   <- if (SMOKE) "output/29_bright_cells_smoke" else "output/29_bright_cells"
FIG   <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
set.seed(29)

SRC_MIN <- 10L
HRS3    <- c("IL13", "BATF3", "TNFRSF8")
CTX     <- c("IRF4", "LTA", "CCR7", "LAMP3", "CD3E", "CD79A", "LYZ")
GENES   <- c("CCL17", "CCL22", HRS3, CTX)
TIER_BREAKS <- c(-Inf, 9.5, 19.5, 49.5, 99.5, Inf)
TIERS   <- c("<10", "10-19", "20-49", "50-99", "100+")
ROI_TIERS <- list(high = c(100, Inf), medium = c(50, 99), low = c(10, 49))
N_CAND  <- 60L                     # candidates per group; the renderer keeps 20
MIN_SEP_PX <- 150 / (COSMX_MM_PER_PX * 1000)   # 150 um in canvas px
SLIDE   <- c(IB1 = 2L, RB1 = 1L)   # c_2_ = IB1 (CLAUDE.md); checked against `samp`

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
## Gene-set hues. Hodgkin blue for the HRS-leaning set (the project's Hodgkin hue);
## FL orange is reserved and deliberately not used.
SET_COL <- c(HRS = "#0b3d7a", ambiguous = "#7a5aa6", DC = "#2a9d8f", other = "#9a968d")
OR_BREAKS <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)  # explicit log breaks (CLAUDE.md, step 25)

## ---- 1. cells, panel counts, reference rates ---------------------------------
reg <- fread(cmd = paste("zcat", REG),
             select = c("cell_id", "samp", "fov_id", "region", "is_defect", "tier_matched",
                        "nCount_RNA", "Area.um2", "d_fl_um", "d_nb_um"))
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
d[, chemo := CCL17 + CCL22]
d[, depth0 := nCount_RNA - chemo]
d <- d[depth0 > 0]
d[, tier := factor(cut(chemo, TIER_BREAKS, labels = TIERS), levels = TIERS)]
d[, depth_q := cut(depth0, unique(quantile(depth0, 0:10 / 10)), include.lowest = TRUE,
                   labels = FALSE), by = samp]

## reference rate per count, per sample x depth decile, from the < 10 cells
ref <- d[tier == "<10", c(list(dep = sum(as.numeric(depth0))),
                          lapply(.SD, function(z) sum(as.numeric(z)))),
         by = .(samp, depth_q), .SDcols = HRS3]
for (g in HRS3) ref[, (paste0("rate_", g)) := get(g) / dep]
src <- d[chemo >= SRC_MIN]
src <- merge(src, ref[, c("samp", "depth_q", paste0("rate_", HRS3)), with = FALSE],
             by = c("samp", "depth_q"), all.x = TRUE)
stopifnot(!anyNA(src[[paste0("rate_", HRS3[1])]]))
for (g in HRS3) {
  src[, (paste0("exp_", g)) := depth0 * get(paste0("rate_", g))]
  src[, (paste0("lr_", g)) := log2((get(g) + 0.5) / (get(paste0("exp_", g)) + 0.5))]
}
src[, hrs_score := rowMeans(.SD), .SDcols = paste0("lr_", HRS3)]
src[, k3 := (IL13 > 0) + (BATF3 > 0) + (TNFRSF8 > 0)]
src[, hot := region == "cHL_hotspot"]

## canvas px from FOV-local centroids (the analytic transform; CLAUDE.md "Images")
meta <- rbindlist(lapply(names(SLIDE), function(s) {
  m <- fread(file.path(FLAT, s, paste0(s, "_metadata_file.csv.gz")),
             select = c("fov", "cell_ID", "CenterX_local_px", "CenterY_local_px"))
  m[, cell_id := sprintf("c_%d_%d_%d", SLIDE[[s]], fov, cell_ID)]
  m[cell_id %in% src[samp == s]$cell_id][, samp := s]
}))
for (s in names(SLIDE)) {
  i <- meta$samp == s
  xy <- cosmx_local_to_px(s, meta$fov[i], meta$CenterX_local_px[i], meta$CenterY_local_px[i])
  meta[i, `:=`(X_px = xy$X_px, Y_px = xy$Y_px)]
}
src <- merge(src, meta[, .(cell_id, X_px, Y_px)], by = "cell_id")
keep <- c("cell_id", "samp", "fov", "cell_ID", "fov_id", "X_px", "Y_px", "region", "hot",
          "Area.um2", "nCount_RNA", "depth0", "chemo", "tier", "CCL17", "CCL22",
          HRS3, paste0("exp_", HRS3), "hrs_score", "k3", CTX)
src <- src[, ..keep]
src[, (paste0("exp_", HRS3)) := lapply(.SD, round, 3), .SDcols = paste0("exp_", HRS3)]
src[, hrs_score := round(hrs_score, 4)]
fwrite(src[order(samp, -chemo)], file.path(OUT, "source_cells.csv.gz"))
cat(sprintf("source cells: %s\n", paste(src[, .N, by = samp][, sprintf("%s %d", samp, N)], collapse = ", ")))

## ---- 2. ROI candidates (IB1) ----------------------------------------------------
pick_spaced <- function(C, n, exclude = character()) {
  C <- C[!cell_id %in% exclude]
  got <- integer()
  for (i in seq_len(nrow(C))) {
    if (length(got) && any((C$X_px[got] - C$X_px[i])^2 + (C$Y_px[got] - C$Y_px[i])^2 <
                           MIN_SEP_PX^2)) next
    got <- c(got, i)
    if (length(got) == n) break
  }
  C[got]
}
## ---- 2a. every IB1 chemokine molecule, canvas px ------------------------------
## One awk pass over the flat transcript file; FOV-local px through the same analytic
## transform as the cell centroids. Unassigned molecules (cell_ID 0) are kept -- they
## are part of the field. Read BEFORE candidate selection so the smoke run, which only
## reads the file's first lines, can restrict its candidates to the FOVs those lines
## cover -- otherwise its ROIs carry no molecules and the contour code never runs.
tx_cmd <- sprintf("zcat %s/IB1/IB1_tx_file.csv.gz%s | awk -F, 'NR==1 || $9==\"CCL17\" || $9==\"CCL22\"' | cut -d, -f1,2,4,5,9,10",
                  FLAT, if (SMOKE) " | head -3000000" else "")
tx <- fread(cmd = tx_cmd)
stopifnot(identical(names(tx), c("fov", "cell_ID", "x_local_px", "y_local_px", "target", "CellComp")))
xy <- cosmx_local_to_px("IB1", tx$fov, tx$x_local_px, tx$y_local_px)
tx[, `:=`(X_px = round(xy$X_px, 1), Y_px = round(xy$Y_px, 1))]
fwrite(tx[, .(X_px, Y_px, target, fov, cell_ID)], file.path(OUT, "ib1_chemokine_tx.csv.gz"))
cat(sprintf("IB1 chemokine molecules: %d (%d unassigned)\n", nrow(tx), tx[cell_ID == 0, .N]))
## cross-check against the step-15 per-FOV pull, where it exists (same molecules, two routes)
chk <- "output/15_imagej/field_test/IB1_55_162_tx.csv"
if (file.exists(chk) && !SMOKE) {
  a <- fread(chk)
  b <- tx[X_px >= min(a$X_px) & X_px <= max(a$X_px) & Y_px >= min(a$Y_px) & Y_px <= max(a$Y_px)]
  cat(sprintf("cross-check vs per-FOV pull: %d vs %d molecules in the same box\n", nrow(a), nrow(b)))
}
TX_FOVS <- unique(tx$fov)
## The smoke read stops mid-file, so its LAST FOV is cut part-way: its cells would have
## only some of their molecules (the renderer's own-count check flags exactly that).
if (SMOKE) TX_FOVS <- setdiff(TX_FOVS, max(TX_FOVS))

## ---- 2b. candidates ------------------------------------------------------------
cand <- list()
for (t in names(ROI_TIERS)) {
  r <- ROI_TIERS[[t]]
  C <- src[samp == "IB1" & chemo >= r[1] & chemo <= r[2]]
  if (SMOKE) C <- C[fov %in% TX_FOVS]
  h <- pick_spaced(C[order(-hrs_score, -chemo)], N_CAND)
  ## the random group must not reuse the HRS group's first 20 (what is actually shown)
  rnd <- pick_spaced(C[sample(.N)], N_CAND, exclude = h$cell_id[seq_len(min(20, nrow(h)))])
  cand[[paste0(t, "_hrs")]]    <- h[, `:=`(group = paste0(t, "_hrs"), rank = .I)]
  cand[[paste0(t, "_random")]] <- rnd[, `:=`(group = paste0(t, "_random"), rank = .I)]
  cat(sprintf("%-7s tier: %5d cells; candidates hrs %d, random %d; hrs-score of hrs #20 %.2f vs tier median %.2f\n",
              t, nrow(C), nrow(h), nrow(rnd), h$hrs_score[min(20, nrow(h))], median(C$hrs_score)))
}
cand <- rbindlist(cand)
setcolorder(cand, c("group", "rank"))
fwrite(cand, file.path(OUT, "roi_candidates.csv"))

## ---- 4. figures ------------------------------------------------------------------
ex <- fread(file.path(S28, "tier_excess.csv"))
ex <- ex[samp == "IB1" & tier %in% TIERS[-1]]
ex[, tier := factor(tier, levels = TIERS[-1])]
ex[, panel := factor(fcase(set == "HRS", "HRS-leaning (pre-declared)",
                           set == "ambiguous", "HRS or activated APC (shared)",
                           set == "DC", "Dendritic-cell lineage",
                           default = "Macrophage / B / T"),
                     levels = c("HRS-leaning (pre-declared)", "HRS or activated APC (shared)",
                                "Dendritic-cell lineage", "Macrophage / B / T"))]
ex[, col := fifelse(set %in% names(SET_COL), set, "other")]
## Label the genes that matter for the claim; the rest stay as thin grey context.
SHOW <- c("IL13", "TNFRSF8", "LTA", "FUT4", "BATF3", "IRF4", "CD80", "CD83",
          "LAMP3", "XCR1", "CLEC9A", "CD1C", "IDO1", "CCR7", "LYZ", "CD68", "CD79A", "CD3E")
ex[, shown := gene %in% SHOW]
lab <- ex[tier == "100+" & shown == TRUE]
p1 <- ggplot(ex, aes(tier, excess, group = gene)) +
  geom_hline(yintercept = 1, colour = AXIS, linewidth = .4) +
  geom_line(data = ex[shown == FALSE], colour = "#d8d6cf", linewidth = .35) +
  geom_errorbar(data = ex[shown == TRUE], aes(ymin = lo, ymax = hi, colour = col),
                width = .08, linewidth = .3, alpha = .6) +
  geom_line(data = ex[shown == TRUE], aes(colour = col), linewidth = .8) +
  geom_point(data = ex[shown == TRUE], aes(colour = col), size = 1.4) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene, colour = col), size = 2.9,
                           nudge_x = .35, direction = "y", hjust = 0, segment.size = .2,
                           min.segment.length = 0, seed = 1) +
  scale_y_continuous(trans = "log2", breaks = OR_BREAKS, labels = function(x) paste0(x, "x")) +
  scale_x_discrete(expand = expansion(add = c(.3, 1.1))) +
  scale_colour_manual(values = SET_COL, guide = "none") +
  facet_wrap(~ panel, nrow = 1) +
  labs(title = "HRS-leaning genes climb with chemokine brightness; dendritic-cell lineage genes do not",
       subtitle = "IB1 source cells by CCL17+CCL22 count tier. Excess = observed / expected counts, expected from cells with < 10 chemokine counts in the same depth decile",
       x = "CCL17+CCL22 counts per cell", y = "Excess over depth-matched reference (log scale)",
       caption = paste("Depth excludes the chemokine counts. Bars: 95% FOV cluster bootstrap. Gene sets pre-declared in step 27; CCR7 was declared DC but is also expressed by HRS cells,",
                       "and it climbs.\nRB1 has 9 cells >= 50 and none >= 100, so it has no comparable bright tiers. Source: output/28_bright_sources/tier_excess.csv.")) +
  theme_viz(10) + theme(strip.text = element_text(colour = INK2, face = "bold", hjust = 0))
sf("fig29_brightness_gene_excess.png", p1, 13, 5.2)

co  <- fread(file.path(S28, "cooccurrence.csv"))[samp == "IB1" & !grepl("\\[", combo)]
nul <- fread(file.path(S28, "cooccurrence_null.csv.gz"))[samp == "IB1"]
ORDER <- c("IL13+BATF3", "IL13+TNFRSF8", "BATF3+TNFRSF8", "IL13+BATF3+TNFRSF8",
           "CD79A+MS4A1", "CD3D+CD3E", "CCR7+LAMP3")
KIND <- c(rep("HRS-leaning genes", 4), rep("Lineage controls", 2), "Activated-DC pair")
co[, combo := factor(combo, levels = ORDER)]; nul[, combo := factor(combo, levels = ORDER)]
co[, kind := factor(KIND[match(combo, ORDER)], levels = unique(KIND))]
nul[, kind := factor(KIND[match(combo, ORDER)], levels = unique(KIND))]
KCOL <- c("HRS-leaning genes" = "#0b3d7a", "Lineage controls" = "#52514e", "Activated-DC pair" = "#2a9d8f")
p2 <- ggplot() +
  geom_hline(yintercept = 1, colour = AXIS, linewidth = .4) +
  geom_violin(data = nul, aes(combo, null_oe), fill = "#e9e7e1", colour = NA, scale = "width",
              width = .7) +
  geom_boxplot(data = nul, aes(combo, null_oe), width = .12, outlier.shape = NA,
               colour = MUTED, fill = NA, linewidth = .3) +
  geom_errorbar(data = co, aes(combo, ymin = oe_lo, ymax = oe_hi, colour = kind), width = .15,
                linewidth = .6) +
  geom_point(data = co, aes(combo, oe, colour = kind), size = 2.8) +
  geom_text(data = co, aes(combo, oe_hi, label = sprintf("%.2f", oe), colour = kind),
            vjust = -0.8, size = 3) +
  geom_text(data = co, aes(combo, 0.72, label = sprintf("null p = %s", signif(p_null, 2))),
            colour = MUTED, size = 2.6) +
  scale_colour_manual(values = KCOL, name = NULL) +
  scale_x_discrete(labels = function(x) gsub("\\+", " + ", x)) +
  coord_cartesian(ylim = c(0.7, 1.75)) +
  labs(title = "IL13, BATF3 and CD30 co-occur in the same cells at lineage-like strength; the activated-DC pair does not",
       subtitle = "IB1 source cells (CCL17+CCL22 >= 10, n = 16,590). Co-positivity (>= 1 count) observed / expected under independence within depth x brightness strata",
       x = NULL, y = "Observed / expected co-positivity",
       caption = paste("Point and bars: observed, 95% FOV cluster bootstrap. Grey: 1,000 random panel gene pairs / triples with each member's detection rate matched within 1.25x.",
                       "\nBATF3 is detected in 52% of sources, so only 13 panel genes match it and its null is thin. TNFRSF8 = CD30. Source: output/28_bright_sources/cooccurrence*.csv.")) +
  theme_viz(10) + theme(legend.position = "top")
sf("fig29_cooccurrence_vs_null.png", p2, 11.5, 5.4)

k3 <- fread(file.path(S28, "hrs3_k_of_3_by_tier.csv"))[n >= 50]
k3[, tier := factor(tier, levels = TIERS[-1])]
kl <- melt(k3, id.vars = c("samp", "tier", "n"), measure.vars = c("k0", "k1", "k2", "k3"),
           variable.name = "k", value.name = "share")
kl[, k := factor(sub("k", "", k), levels = c("3", "2", "1", "0"))]
kl[, tier := factor(tier, levels = TIERS[-1])]
KC <- c("3" = "#0b3d7a", "2" = "#4f78ad", "1" = "#a9bdd8", "0" = "#e9e7e1")
p3 <- ggplot(kl, aes(tier, share, fill = k)) +
  geom_col(width = .7) +
  geom_text(data = k3, aes(tier, 1.03, label = sprintf("n = %s", format(n, big.mark = ","))),
            inherit.aes = FALSE, size = 2.7, colour = MUTED) +
  facet_grid(~ samp, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = KC, name = "of IL13 / BATF3 / CD30 detected") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, .25)) +
  labs(title = "Share of source cells carrying 0-3 of the HRS-leaning genes, by brightness",
       x = "CCL17+CCL22 counts per cell", y = "Share of source cells",
       caption = "Tiers with < 50 cells omitted (RB1 has 9 cells >= 50, none >= 100). Detection rises with depth, so this is a description, not a call. Source: output/28_bright_sources/hrs3_k_of_3_by_tier.csv.") +
  theme_viz(10)
sf("fig29_k_of_3_by_tier.png", p3, 8, 4.6)

sc4 <- src[samp == "IB1"]
sc4[, roi_tier := fcase(chemo >= 100, "high (100+)", chemo >= 50, "medium (50-99)",
                        default = "low (10-49)")]
sc4[, roi_tier := factor(roi_tier, levels = c("low (10-49)", "medium (50-99)", "high (100+)"))]
cut20 <- cand[rank == 20 & grepl("_hrs$", group),
              .(roi_tier = factor(c(high = "high (100+)", medium = "medium (50-99)",
                                    low = "low (10-49)")[sub("_hrs", "", group)],
                                  levels = levels(sc4$roi_tier)), hrs_score)]
p4 <- ggplot(sc4, aes(roi_tier, hrs_score)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_violin(fill = "#dfe6ef", colour = NA, scale = "width", width = .8) +
  geom_boxplot(width = .12, outlier.shape = NA, colour = "#0b3d7a", fill = NA, linewidth = .35) +
  geom_point(data = cut20, colour = "#0b3d7a", shape = 95, size = 9) +
  labs(title = "HRS-leaning score by brightness tier, with the cut that selects each tier's top 20",
       subtitle = "Score = mean over IL13, BATF3, TNFRSF8 of log2(observed / expected), expected from depth-matched cells with < 10 chemokine counts",
       x = NULL, y = "HRS-leaning score (log2)",
       caption = "Bar: score of the 20th candidate in the *_hrs ROI group (before the renderer drops windows on blank morphology). Source: output/29_bright_cells/source_cells.csv.gz.") +
  theme_viz(10)
sf("fig29_hrs_score_by_tier.png", p4, 8, 4.6)

## ---- 5. are bright cells excluded from FL neighbourhoods? ----------------------
## Bright cells per distance band from the nearest FL-like (kappa-restricted)
## neighbourhood member -- step 16's d_fl_um, footprint <= 40 um, boundary 40-100 um.
## Control axis: distance to ANY called light-chain neighbourhood (d_nb_um), which
## separates "FL-specific" from "any B-cell cluster". FOV cluster bootstrap.
##
## DEPTH IS THE CONFOUND, and it runs toward ENRICHMENT: FL tissue is deeper (median
## depth0 ~500 in the FL bands vs ~290 beyond 800 um), and a count threshold is easier
## to reach in a deeper cell. So the plotted statistic is DEPTH-STANDARDISED: observed
## bright cells / expected, expected = sum over the band's cells of P(bright | the
## cell's depth0 decile), the probability taken over all IB1 cells here. The raw rate
## ratio is kept beside it. Depth-adjusting may over-correct -- bright cells are deep
## partly because they are large activated cells -- so the truth lies between the two.
## Rates are per defect-FREE cell, and the defect drop differs by region (fact 9).
N_BOOT <- if (SMOKE) 20L else 1000L
BANDS <- c(0, 40, 100, 200, 400, 800, Inf)
BAND_LAB <- c("<=40\n(FL footprint)", "40-100\n(boundary)", "100-200", "200-400", "400-800", ">800")
AXES <- c(d_fl_um = "Distance to an FL-like (kappa-restricted) neighbourhood",
          d_nb_um = "Control: distance to ANY light-chain neighbourhood")
THR <- c(n10 = ">= 10 counts", n50 = ">= 50 counts", n100 = ">= 100 counts")
CUT <- c(n10 = 10, n50 = 50, n100 = 100)
fl <- d[samp == "IB1" & !is.na(d_fl_um) & !is.na(d_nb_um)]
fl[, dq := cut(depth0, unique(quantile(depth0, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
for (k in names(THR)) {
  fl[, (k) := as.integer(chemo >= CUT[[k]])]
  fl[, (paste0("e_", k)) := mean(get(k)), by = dq]          # P(bright | depth decile)
}
long <- rbindlist(lapply(names(AXES), function(v)
  fl[, c(list(axis = v, band = cut(get(v), BANDS, labels = BAND_LAB, include.lowest = TRUE)),
         .SD), .SDcols = c("fov_id", "depth0", names(THR), paste0("e_", names(THR)))]))
SUMS <- c("cells", names(THR), paste0("e_", names(THR)))
fagg <- long[, c(list(cells = .N), lapply(.SD, sum)), by = .(axis, band, fov_id),
             .SDcols = c(names(THR), paste0("e_", names(THR)))]
band_stats <- function(A) {
  b <- A[, lapply(.SD, function(z) sum(z * w)), by = .(axis, band), .SDcols = SUMS]
  tot <- A[, lapply(.SD, function(z) sum(z * w)), by = axis, .SDcols = SUMS]
  b <- merge(b, tot, by = "axis", suffixes = c("", "_all"))
  rbindlist(lapply(names(THR), function(k)
    b[, .(axis, band, thr = k, n = get(k), cells, expected = get(paste0("e_", k)),
          raw = (get(k) / cells) / (get(paste0(k, "_all")) / cells_all),
          adj = get(k) / get(paste0("e_", k)))]))
}
fovs <- unique(fagg$fov_id)
R0 <- band_stats(fagg[, w := 1])
BB <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
  w <- data.table(fov_id = sample(fovs, replace = TRUE))[, .(w = .N), by = fov_id]
  band_stats(merge(fagg[, !"w"], w, by = "fov_id"))[, .(axis, band, thr, raw, adj)]
}))
ci <- BB[, .(raw_lo = quantile(raw, .025, na.rm = TRUE), raw_hi = quantile(raw, .975, na.rm = TRUE),
             adj_lo = quantile(adj, .025, na.rm = TRUE), adj_hi = quantile(adj, .975, na.rm = TRUE)),
         by = .(axis, band, thr)]
fx <- merge(R0, ci, by = c("axis", "band", "thr"))
dep <- long[, .(depth0_med = as.numeric(median(depth0))), by = .(axis, band)]
fx <- merge(fx, dep, by = c("axis", "band"))
num <- c("expected", "raw", "adj", "raw_lo", "raw_hi", "adj_lo", "adj_hi")
fx[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(fx[order(axis, thr, band)], file.path(OUT, "fl_exclusion.csv"))
cat("\n== bright cells by distance from FL: raw and depth-standardised ratios, FOV bootstrap\n")
print(fx[thr != "n100", .(axis, band, thr, cells, n, raw, adj, adj_lo, adj_hi, depth0_med)])

## The >= 100 tier stays in the CSV but is not plotted: ~390 cells over six bands and
## two axes leaves most bands with a handful or none, too sparse to read a trend from.
PLOT_THR <- c("n10", "n50")
fp <- fx[thr %in% PLOT_THR]
fp[, thr_lab := factor(THR[thr], levels = THR[PLOT_THR])]
fp[, axis_lab := factor(AXES[axis], levels = AXES)]
fp[, band := factor(band, levels = BAND_LAB)]
TCOL <- setNames(c("#7d9cc6", "#0b3d7a"), THR[PLOT_THR])
FLOOR <- 0.0625                    # a zero has no log; it is drawn here, open
fp[, `:=`(adj_plot = pmax(adj, FLOOR), lo_plot = pmax(adj_lo, FLOOR), hi_plot = pmax(adj_hi, FLOOR),
          raw_plot = pmax(raw, FLOOR), zero = n == 0)]
dg <- position_dodge(width = .35)
p5 <- ggplot(fp, aes(band, colour = thr_lab, group = thr_lab)) +
  geom_hline(yintercept = 1, colour = AXIS, linewidth = .4) +
  geom_line(aes(y = raw_plot), linewidth = .45, linetype = "22", alpha = .55, position = dg) +
  geom_errorbar(aes(ymin = lo_plot, ymax = hi_plot), width = .12, linewidth = .4, position = dg) +
  geom_line(aes(y = adj_plot), linewidth = .85, position = dg) +
  geom_point(aes(y = adj_plot, shape = zero), size = 1.9, fill = SURF, position = dg) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
  facet_wrap(~ axis_lab, nrow = 1) +
  scale_y_continuous(trans = "log2", breaks = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4),
                     labels = function(x) paste0(x, "x")) +
  scale_colour_manual(values = TCOL, name = "CCL17+CCL22 per cell") +
  labs(title = "Chemokine-bright cells are not excluded from FL neighbourhoods",
       subtitle = "IB1. Solid: observed / expected bright cells per distance band, expected from each cell's depth decile (1x = no preference). Dashed: raw rate ratio, not depth-adjusted",
       x = "Distance band (um)", y = "Depth-standardised ratio (log scale)",
       caption = paste("95% FOV cluster bootstrap on the solid line. OPEN points = no bright cell in that band (a zero, drawn at the floor).",
                       "\nFL footprint and boundary are step 16's (40 / 100 um from a kappa-restricted member). FL tissue is deeper, so the raw ratio overstates any enrichment there;",
                       "\ndepth-adjusting may over-correct (bright cells are deep partly because they are large). Rates are per defect-free cell, and the defect drop differs by region",
                       "\n(CLAUDE.md fact 9). The >= 100 tier is too sparse per band to plot; it is in output/29_bright_cells/fl_exclusion.csv.")) +
  theme_viz(10) + theme(strip.text = element_text(colour = INK2, face = "bold", hjust = 0))
sf("fig29_fl_exclusion.png", p5, 12.5, 5.3)

## ---- 6. chemokine count against its rank --------------------------------------------
## Every cell plotted has >= 10 counts, so the source list is the whole input. The x axis
## stops at RANK_MAX: by then IB1 is down near the 20s and the curve has flattened toward
## the 10-count floor; past it the line is just a long flat tail.
RANK_MAX <- 5000L
rk <- src[order(samp, -chemo)][, rank := seq_len(.N), by = samp][rank <= RANK_MAX]
SAMP_COLS <- c(IB1 = "#2a78d6", RB1 = "#eb6834")   # fixed sample hues (CLAUDE.md)
cuts <- rbindlist(lapply(c(50, 100), function(k)
  src[, .(cut = k, n = sum(chemo >= k)), by = samp]))
cut_lab <- cuts[, .(lab = paste(sprintf("%s: %s cells", samp, format(n, big.mark = ",")),
                                collapse = "   ")), by = cut]
cut_lab[, lab := sprintf(">= %d counts   %s", cut, lab)]
mark <- cuts[n > 0 & n <= RANK_MAX]                 # where each sample's curve crosses
p6 <- ggplot(rk, aes(rank, chemo, colour = samp)) +
  geom_hline(data = cut_lab, aes(yintercept = cut), colour = INK2, linewidth = .35, linetype = "22") +
  geom_text(data = cut_lab, aes(x = RANK_MAX, y = cut, label = lab), inherit.aes = FALSE,
            hjust = 1, vjust = -0.5, size = 3, colour = INK2) +
  geom_line(linewidth = .8) +
  geom_point(data = mark, aes(x = n, y = cut, colour = samp), size = 2.4) +
  geom_text(data = mark, aes(x = n, y = cut, label = format(n, big.mark = ","), colour = samp),
            vjust = 1.8, size = 3, show.legend = FALSE) +
  scale_y_continuous(trans = "log10", breaks = c(10, 20, 50, 100, 200, 500),
                     limits = c(10, NA)) +
  scale_x_continuous(breaks = seq(0, RANK_MAX, 1000), labels = scales::comma,
                     expand = expansion(mult = c(.01, .02))) +
  scale_colour_manual(values = SAMP_COLS, name = NULL) +
  labs(title = "CCL17+CCL22 counts per cell against their rank",
       subtitle = sprintf("Defect-free matched-tier cells with >= 10 counts, the top %s per sample. Points mark where each sample's curve crosses a cutoff",
                          format(RANK_MAX, big.mark = ",")),
       x = "Rank (brightest first)", y = "CCL17+CCL22 counts (log scale)",
       caption = sprintf("IB1 has %s cells >= 10 in all; RB1 %s, so RB1's line ends at its own %s. Source: output/29_bright_cells/source_cells.csv.gz.",
                         format(src[samp == "IB1", .N], big.mark = ","),
                         format(src[samp == "RB1", .N], big.mark = ","),
                         format(src[samp == "RB1", .N], big.mark = ","))) +
  theme_viz(10, xgrid = TRUE)
sf("fig29_chemokine_rank.png", p6, 9, 5)

cat("\ndone:", OUT, "\n")
