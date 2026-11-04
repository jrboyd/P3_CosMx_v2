## ============================================================================
## Step 31 -- Radial CD45 profiles around candidate HRS cells (IB1), segmentation-free.
##
## WHY. Step 30 showed per-cell segmented CD45 intensity cannot tell CD45+ from CD45-
## lineages (AUC 0.52): in packed lymphoid tissue a cell's mean is mostly its
## neighbours' membrane. A radial profile from the centroid reads the IMAGE directly and
## asks the shape question instead -- dim centre, bright ring of rosetting lymphocytes? --
## without trusting any boundary.
##
## Three stages, chained by scripts/31_submit_radial_cd45.sh:
##   select   (this script)  500 candidates + 500 area-matched controls, canvas px
##   profile  (31_radial_cd45_profiles.py)  CD45 + DNA mean per 5 um band, 0-100 um
##   analyze  (this script)  average profiles, class test, heatmaps
##
## CELLS (fixed before running)
##   candidates  seeded random 500 of the IB1 cells with CCL17+CCL22 >= 50 (step 29's
##               source list; >= 50 is ~300x Hodgkin-specific, CLAUDE.md steps 27-29)
##   controls    500 IB1 defect-free matched-tier cells with chemo < 10, sampled within
##               the candidates' Area deciles so the two groups have the same size
##               distribution. Any large cell has a DNA-filled, CD45-dimmer middle, so a
##               dim-centre profile means nothing without this control.
## NORMALISATION  each profile / the cell's own 50-100 um mean (IB1 staining varies by
##               region, CLAUDE.md facts 2/9), log2. Raw means are stored too.
## CLASS TEST    k-means on the normalised candidate profiles, k = 2..6 (nstart 25).
##               Classes are DETECTED at k only if mean silhouette exceeds the 95th pct of
##               the same k-means on N_NULL draws from ONE multivariate normal with the
##               data's own mean and covariance (a single smeary distribution), AND every
##               cluster holds >= 5% of cells. Otherwise: one distribution, reported so.
##
##   Rscript scripts/31_script_radial_cd45.R select|analyze [smoke]
## ============================================================================
suppressPackageStartupMessages(library(data.table))
args  <- commandArgs(trailingOnly = TRUE)
MODE  <- args[1]
SMOKE <- "smoke" %in% args
OUT   <- if (SMOKE) "output/31_radial_cd45_smoke" else "output/31_radial_cd45"
FIG   <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
N_EACH <- if (SMOKE) 12L else 500L
BAND_UM <- 5; R_MAX <- 100
NB <- R_MAX / BAND_UM
OUTER <- which(seq_len(NB) * BAND_UM > 50)        # bands covering 50-100 um
set.seed(31)

if (MODE == "select") {
  source("scripts/lib/cosmx_tiff.R")
  src <- fread("output/29_bright_cells/source_cells.csv.gz")[samp == "IB1"]
  cand <- src[chemo >= 50][sample(.N, min(N_EACH, .N))]
  cand[, group := "candidate"]
  ## controls: defect-free matched IB1 cells with chemo < 10, area-matched by decile
  reg <- fread(cmd = "zcat output/16_regions/matched/cell_regions.csv.gz",
               select = c("cell_id", "samp", "fov_id", "is_defect", "tier_matched",
                          "Area.um2", "chemo_ct"))[
    samp == "IB1" & tier_matched == TRUE & is_defect == FALSE & chemo_ct < 10]
  br <- unique(quantile(cand$Area.um2, 0:10 / 10))
  cand[, adec := cut(Area.um2, br, include.lowest = TRUE, labels = FALSE)]
  reg[, adec := cut(Area.um2, br, include.lowest = TRUE, labels = FALSE)]
  need <- cand[, .N, by = adec]
  ctrl <- rbindlist(lapply(seq_len(nrow(need)), function(i)
    reg[adec == need$adec[i]][sample(.N, need$N[i])]))
  ## canvas px for the controls from FOV-local centroids (the analytic transform)
  m <- fread("input/CosMx/export1/flatFiles/IB1/IB1_metadata_file.csv.gz",
             select = c("fov", "cell_ID", "CenterX_local_px", "CenterY_local_px"))
  m[, cell_id := sprintf("c_2_%d_%d", fov, cell_ID)]     # c_2_ = IB1 (CLAUDE.md)
  ctrl <- merge(ctrl, m, by = "cell_id")
  stopifnot(nrow(ctrl) == sum(need$N))
  xy <- cosmx_local_to_px("IB1", ctrl$fov, ctrl$CenterX_local_px, ctrl$CenterY_local_px)
  ctrl[, `:=`(X_px = xy$X_px, Y_px = xy$Y_px, group = "control", chemo = chemo_ct,
              hrs_score = NA_real_, k3 = NA_integer_)]
  cells <- rbind(cand[, .(cell_id, group, fov_id, X_px, Y_px, Area.um2, chemo, hrs_score, k3)],
                 ctrl[, .(cell_id, group, fov_id, X_px, Y_px, Area.um2, chemo, hrs_score, k3)])
  fwrite(cells, file.path(OUT, "cells.csv"))
  cat(sprintf("wrote %d candidates and %d area-matched controls; median area %.1f vs %.1f um^2\n",
              nrow(cand), nrow(ctrl), median(cand$Area.um2), median(ctrl$Area.um2)))
  quit(save = "no")
}

stopifnot(MODE == "analyze")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(cluster); library(MASS) })
N_NULL <- if (SMOKE) 5L else 50L
N_BOOT <- if (SMOKE) 20L else 1000L
K_RANGE <- 2:6; MIN_CLUSTER_FRAC <- 0.05

## ---- house style (canonical copy: 07_script_lightchain_neighborhoods.R:126-147) ----
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#d03b3b"; DIV_MID<-"#f0efec"; DIV_LO<-"#184f95"
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
GCOL <- c(candidate = "#0b3d7a", control = "#9a968d")

pr <- fread(file.path(OUT, "profiles.csv.gz"))
cells <- fread(file.path(OUT, "cells.csv"))
pr <- merge(cells, pr, by = c("cell_id", "group"))
cat(sprintf("profiled: %s; dropped (blank / off-canvas): %s\n",
            paste(pr[, .N, by = group][, sprintf("%s %d", group, N)], collapse = ", "),
            paste(pr[ok == FALSE, .N, by = .(group, reason)][, sprintf("%s %s %d", group, reason, N)],
                  collapse = ", ")))
pr <- pr[ok == TRUE]
CD <- sprintf("cd45_b%02d", seq_len(NB) - 1L)
DN <- sprintf("dna_b%02d", seq_len(NB) - 1L)
band_mid <- (seq_len(NB) - 0.5) * BAND_UM
band_lab <- sprintf("%d-%d", (seq_len(NB) - 1L) * BAND_UM, seq_len(NB) * BAND_UM)

raw <- as.matrix(pr[, ..CD])
outer <- rowMeans(raw[, OUTER, drop = FALSE])
rel <- log2(raw / outer)                                   # 0 = the cell's own 50-100 um level
colnames(rel) <- band_lab
dna <- as.matrix(pr[, ..DN]); dna_rel <- log2(dna / rowMeans(dna[, OUTER, drop = FALSE]))
## per-cell summaries that read the rosette question directly
pr[, `:=`(centre_rel = rel[, 1], ring_rel = rowMeans(rel[, 2:4, drop = FALSE]),   # 0-5 vs 5-20 um
          rosette_contrast = rowMeans(rel[, 2:4, drop = FALSE]) - rel[, 1],
          cd45_outer = outer)]
fwrite(cbind(pr[, .(cell_id, group, fov_id, Area.um2, chemo, hrs_score, k3, cd45_outer,
                    centre_rel, ring_rel, rosette_contrast)],
             as.data.table(round(rel, 4))), file.path(OUT, "profiles_relative.csv.gz"))

## ---- 1. average profiles, FOV cluster bootstrap ------------------------------------
avg <- rbindlist(lapply(c("candidate", "control"), function(g) {
  i <- which(pr$group == g); f <- pr$fov_id[i]; uf <- unique(f); idx <- split(i, f)
  m_rel <- colMeans(rel[i, ]); m_raw <- colMeans(raw[i, ]); m_dna <- colMeans(dna_rel[i, ])
  B <- t(replicate(N_BOOT, colMeans(rel[unlist(idx[as.character(sample(uf, replace = TRUE))]), , drop = FALSE])))
  data.table(group = g, band = band_lab, r_um = band_mid, rel = m_rel, raw = m_raw, dna_rel = m_dna,
             lo = apply(B, 2, quantile, .025), hi = apply(B, 2, quantile, .975), n = length(i))
}))
fwrite(avg, file.path(OUT, "average_profile.csv"))
## candidate minus control, per band, FOV bootstrap of each group independently
dif <- avg[group == "candidate", .(band, r_um, cand = rel)][avg[group == "control", .(band, ctrl = rel)], on = "band"]
dif[, d := cand - ctrl]
cat("\n== average log2 CD45 relative to own 50-100 um, by band\n")
print(dcast(avg, band + r_um ~ group, value.var = "rel")[order(r_um)][, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

## ---- 2. class test on candidate profiles -------------------------------------------
X <- rel[pr$group == "candidate", ]
Xc <- rel[pr$group == "control", ]
D <- dist(X)
sil_of <- function(M, k, d = dist(M)) {
  km <- kmeans(M, centers = k, nstart = 25, iter.max = 100)
  list(sil = mean(silhouette(km$cluster, d)[, 3]), cl = km$cluster,
       min_frac = min(tabulate(km$cluster, k)) / nrow(M))
}
mu <- colMeans(X); S <- cov(X)
kt <- rbindlist(lapply(K_RANGE, function(k) {
  o <- sil_of(X, k, D)
  nul <- replicate(N_NULL, { Z <- mvrnorm(nrow(X), mu, S); sil_of(Z, k)$sil })
  data.table(k = k, silhouette = o$sil, null_median = median(nul), null_p95 = quantile(nul, .95),
             min_cluster_frac = o$min_frac,
             detected = o$sil > quantile(nul, .95) & o$min_frac >= MIN_CLUSTER_FRAC)
}))
fwrite(kt, file.path(OUT, "class_test.csv"))
cat("\n== class test: k-means silhouette vs a single multivariate-normal null\n"); print(kt)
K_BEST <- if (any(kt$detected)) kt[detected == TRUE][which.max(silhouette - null_p95)]$k else NA_integer_
pc <- prcomp(X, center = TRUE, scale. = FALSE)
ve <- summary(pc)$importance[2, 1:3]
cat(sprintf("\nPCA of candidate profiles: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%% of variance\n",
            100 * ve[1], 100 * ve[2], 100 * ve[3]))
pr[group == "candidate", pc1 := pc$x[, 1]]
## orient PC1 so that positive = brighter ring relative to centre, for readability
if (cor(pr[group == "candidate"]$pc1, pr[group == "candidate"]$rosette_contrast) < 0)
  pr[group == "candidate", pc1 := -pc1]
if (!is.na(K_BEST)) {
  cl <- kmeans(X, centers = K_BEST, nstart = 50, iter.max = 200)$cluster
  pr[group == "candidate", cls := cl]
  prof_cl <- as.data.table(X)[, cls := cl][, lapply(.SD, mean), by = cls]
  fwrite(prof_cl, file.path(OUT, "class_profiles.csv"))
  cat(sprintf("\nCLASSES DETECTED at k = %d\n", K_BEST)); print(pr[group == "candidate", .N, by = cls])
} else cat("\nNO CLASSES DETECTED: candidate profiles form one smeary distribution\n")
fwrite(pr[, .(cell_id, group, pc1, cls = if ("cls" %in% names(pr)) cls else NA_integer_)],
       file.path(OUT, "cell_order.csv"))

## ---- 3. figures --------------------------------------------------------------------
f1a <- ggplot(avg, aes(r_um, rel, colour = group, fill = group)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) +
  geom_line(linewidth = .9) + geom_point(size = 1.3) +
  scale_colour_manual(values = GCOL, name = NULL) + scale_fill_manual(values = GCOL, name = NULL) +
  scale_x_continuous(breaks = seq(0, R_MAX, 10)) +
  labs(title = "Average CD45 profile around candidate HRS cells",
       subtitle = sprintf("log2 of each %g um band's mean CD45 over the cell's own 50-100 um mean; 95%% FOV cluster bootstrap", BAND_UM),
       x = "Distance from centroid (um)", y = "log2 CD45 vs own 50-100 um") + theme_viz(10, TRUE)
f1b <- ggplot(avg, aes(r_um, dna_rel, colour = group)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_line(linewidth = .9) + geom_point(size = 1.3) +
  scale_colour_manual(values = GCOL, guide = "none") +
  scale_x_continuous(breaks = seq(0, R_MAX, 10)) +
  labs(subtitle = "DNA, same normalisation -- marks where the nucleus sits",
       x = "Distance from centroid (um)", y = "log2 DNA vs own 50-100 um") + theme_viz(10, TRUE)
f1 <- f1a / f1b + plot_layout(heights = c(2, 1)) +
  plot_annotation(caption = sprintf(
    "Candidates: random %d of IB1 cells with CCL17+CCL22 >= 50. Controls: %d IB1 cells with < 10, matched on area decile. Image read directly at full resolution (no segmentation).\nSource: output/31_radial_cd45/average_profile.csv.",
    sum(pr$group == "candidate"), sum(pr$group == "control")),
    theme = theme(plot.caption = element_text(colour = MUTED, size = 7, hjust = 0), plot.background = element_rect(fill = SURF, colour = NA)))
sf("fig31_average_profile.png", f1, 9, 7)

## heatmap: rows = cells, columns = bands, fill = log2 relative CD45 (blue -> red, red = high)
lim <- max(abs(quantile(rbind(X, Xc), c(.02, .98))))
hm_of <- function(M, ord, title, anno = NULL) {
  M <- M[ord, , drop = FALSE]
  L <- data.table(row = rep(seq_len(nrow(M)), ncol(M)), band = factor(rep(colnames(M), each = nrow(M)), levels = band_lab),
                  v = pmax(pmin(as.vector(M), lim), -lim))
  p <- ggplot(L, aes(band, row, fill = v)) + geom_raster() +
    scale_fill_gradient2(low = DIV_LO, mid = DIV_MID, high = DIV_HI, midpoint = 0, limits = c(-lim, lim),
                         name = "log2 CD45 vs own 50-100 um") +
    scale_y_reverse(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) +
    labs(title = title, x = "Distance band (um)", y = sprintf("cells (n = %d)", nrow(M))) +
    theme_viz(9) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
                         axis.text.y = element_blank(), panel.grid = element_blank())
  if (!is.null(anno)) p <- p + geom_hline(yintercept = anno + .5, colour = INK, linewidth = .5)
  p
}
cand_i <- which(pr$group == "candidate"); ctrl_i <- which(pr$group == "control")
if (!is.na(K_BEST)) {
  ord_c <- order(pr$cls[cand_i], -pr$pc1[cand_i]); brk <- cumsum(tabulate(pr$cls[cand_i]))[-K_BEST]
  ttl <- sprintf("Candidates, %d classes (k-means), ordered by class then PC1", K_BEST)
} else {
  ord_c <- order(-pr$pc1[cand_i]); brk <- NULL
  ttl <- "Candidates, ordered by PC1 (no classes detected)"
}
ctrl_pc <- predict(pc, Xc)[, 1] * sign(cor(pr[group == "candidate"]$pc1, pc$x[, 1]))
h1 <- hm_of(X, ord_c, ttl, brk)
h2 <- hm_of(Xc, order(-ctrl_pc), "Area-matched controls (chemo < 10), same PC1 axis")
verdict <- if (is.na(K_BEST)) "No classes beyond a single smeary distribution" else sprintf("%d classes detected", K_BEST)
hm <- (h1 | h2) + plot_layout(guides = "collect") +
  plot_annotation(title = sprintf("Radial CD45 around candidate HRS cells: %s", verdict),
                  subtitle = sprintf("Rows = cells, columns = %g um distance bands from the centroid. Red = brighter than the cell's own 50-100 um surroundings, blue = dimmer", BAND_UM),
                  caption = sprintf("Class test: k-means k = %d-%d on candidate profiles, silhouette vs %d draws from one multivariate normal with the data's own mean and covariance; detected only above its 95th pct with every cluster >= 5%%. PC1 = %.0f%% of variance.\nColour clamped at +/- %.2f (2nd-98th pct). Source: output/31_radial_cd45/profiles_relative.csv.gz, class_test.csv.",
                                    min(K_RANGE), max(K_RANGE), N_NULL, 100 * ve[1], lim),
                  theme = theme(plot.title = element_text(colour = INK, face = "bold", size = 13),
                                plot.subtitle = element_text(colour = INK2, size = 9),
                                plot.caption = element_text(colour = MUTED, size = 7, hjust = 0),
                                plot.background = element_rect(fill = SURF, colour = NA))) &
  theme(legend.position = "bottom")
sf("fig31_heatmap.png", hm, 12, 9)

ks <- melt(kt, id.vars = "k", measure.vars = c("silhouette", "null_median", "null_p95"))
f3 <- ggplot(ks, aes(k, value, colour = variable, linetype = variable)) +
  geom_line(linewidth = .8) + geom_point(size = 1.6) +
  scale_colour_manual(values = c(silhouette = "#0b3d7a", null_median = MUTED, null_p95 = INK2),
                      labels = c("observed", "null median", "null 95th pct"), name = NULL) +
  scale_linetype_manual(values = c(silhouette = "solid", null_median = "22", null_p95 = "22"),
                        labels = c("observed", "null median", "null 95th pct"), name = NULL) +
  labs(title = "Class test: does k-means find more structure than one Gaussian would?",
       x = "k (clusters)", y = "Mean silhouette") + theme_viz(10, TRUE)
sf("fig31_class_test.png", f3, 7, 4.2)
cat("\ndone:", OUT, "\n")
