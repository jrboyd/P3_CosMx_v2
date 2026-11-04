## ============================================================================
## Step 33 -- Follow-up on step 31's five CD45-profile clusters: what each looks like on
## average (with the controls that share its shape), and whether the clusters -- defined
## on CD45 alone -- separate anything in the RNA.
##
## Step 31 clustered the 5 um radial CD45 profiles of 492 candidate HRS cells (k-means,
## k = 5, a borderline pass of its single-Gaussian null). Two shapes looked HRS-like:
## cluster 1 (dim centre, bright ring at 10-20 um) and cluster 5 (deeply dim centre,
## recovering by 25-35 um). This step asks the independent question: do the clusters
## differ in the step-28 gene framework?
##
## DESIGN (fixed before running)
##   profiles  per cluster: mean candidate profile with a FOV cluster bootstrap ribbon, and
##             the CONTROLS whose profile is nearest that cluster's centroid, so a shape
##             that ordinary cells also produce is visible as such.
##   excess    per gene and cluster: obs / exp counts, exp = depth0 x the gene's rate in
##             IB1 cells with CCL17+CCL22 < 10 in the same depth0 decile (step 28's
##             reference; depth0 excludes the chemokine counts). FOV cluster bootstrap.
##   test      per gene, heterogeneity across clusters = n-weighted variance of log
##             excess, against N_PERM permutations of the cluster labels among the
##             candidates; BH over genes. A permutation at the cell level ignores that
##             neighbouring candidates share tissue -- it is the right null for "do the
##             labels carry RNA information", not for a population claim.
##   gauge     CD3E / CD3D. A bright CD45 ring IS T cells and their transcripts leak into
##             the cell they surround (CLAUDE.md fact 6); IL13, LTA, TNFRSF8 are also made
##             by activated T cells. An HRS-leaning difference counts only where it exceeds
##             the gauge's in the same cluster.
##
##   sbatch --mem=32G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##       scripts/33_script_cluster_followup.R [smoke]
## ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
args  <- commandArgs(trailingOnly = TRUE)
SMOKE <- "smoke" %in% args
S31   <- "output/31_radial_cd45"
OUT   <- if (SMOKE) "output/33_cluster_followup_smoke" else "output/33_cluster_followup"
FIG   <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
FLAT   <- "input/CosMx/export1/flatFiles/IB1"
N_BOOT <- if (SMOKE) 20L else 1000L
N_PERM <- if (SMOKE) 50L else 1000L
set.seed(33)

SETS <- list(
  DC        = c("LAMP3", "CCR7", "CCL19", "IDO1", "CD1C", "CLEC9A", "XCR1", "CLEC10A",
                "FCER1A", "CD1A", "CD207", "MARCKSL1"),
  HRS       = c("TNFRSF8", "FUT4", "IL13", "LTA", "CSF2"),
  ambiguous = c("CD274", "PDCD1LG2", "CD40", "BATF3", "IRF4", "CD80", "CD86", "CD83", "JAK2"),
  macro     = c("CD68", "CD163", "C1QA", "CD14", "LYZ"),
  Bprog     = c("CD19", "CD79A", "MS4A1"),
  T         = c("CD3E", "CD3D"))
TEST  <- unlist(SETS, use.names = FALSE)
GENES <- c("CCL17", "CCL22", TEST)
SET_OF <- setNames(rep(names(SETS), lengths(SETS)), TEST)
## Short names for the five step-31 shapes, from output/31_radial_cd45/class_profiles.csv
CLS_NAME <- c(`1` = "1: dim centre,\nbright ring 10-20", `2` = "2: mildly dim\nthroughout",
              `3` = "3: bright centre", `4` = "4: deeply dim,\nbroad",
              `5` = "5: deeply dim centre,\nring 25-35")

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
GCOL <- c(candidate = "#0b3d7a", control = "#9a968d")
SET_COL <- c(HRS = "#0b3d7a", ambiguous = "#7a5aa6", DC = "#2a9d8f", other = "#9a968d")
OR_BREAKS <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)

## ---- 1. profiles by cluster, controls assigned to the nearest centroid ---------------
rel <- fread(file.path(S31, "profiles_relative.csv.gz"))
ord <- fread(file.path(S31, "cell_order.csv"))
cen <- fread(file.path(S31, "class_profiles.csv"))[order(cls)]
bands <- names(cen)[-1]
band_mid <- (seq_along(bands) - 0.5) * 5
rel <- merge(rel, ord[, .(cell_id, cls)], by = "cell_id")
X <- as.matrix(rel[, ..bands]); C <- as.matrix(cen[, ..bands])
rel[, nearest := apply(X, 1, function(x) which.min(colSums((t(C) - x)^2)))]
stopifnot(rel[group == "candidate", all(nearest == cls)])
rel[group == "control", cls := nearest]
## DNA profiles from the raw table (same normalisation as step 31)
raw <- fread(file.path(S31, "profiles.csv.gz"))[ok == TRUE]
DN <- grep("^dna_b", names(raw), value = TRUE)
dna <- as.matrix(raw[, ..DN]); outer <- rowMeans(dna[, (ncol(dna) / 2 + 1):ncol(dna)])
dna_rel <- log2(dna / outer); colnames(dna_rel) <- bands
dna_rel <- as.data.table(dna_rel)[, cell_id := raw$cell_id]
rel <- merge(rel, dna_rel, by = "cell_id", suffixes = c("", "_dna"))
Xd <- as.matrix(rel[, paste0(bands, "_dna"), with = FALSE])
X  <- as.matrix(rel[, ..bands])

fovs <- unique(rel$fov_id); idx <- split(seq_len(nrow(rel)), rel$fov_id)
prof <- rbindlist(lapply(1:5, function(k) rbindlist(lapply(c("candidate", "control"), function(g) {
  i <- which(rel$cls == k & rel$group == g)
  if (length(i) < 3) return(NULL)
  B <- t(replicate(N_BOOT, {
    rows <- unlist(idx[as.character(sample(fovs, replace = TRUE))], use.names = FALSE)
    rows <- rows[rel$cls[rows] == k & rel$group[rows] == g]
    if (length(rows) < 3) rep(NA_real_, length(bands)) else colMeans(X[rows, , drop = FALSE])
  }))
  data.table(cls = k, group = g, n = length(i), band = bands, r_um = band_mid,
             cd45 = colMeans(X[i, , drop = FALSE]), lo = apply(B, 2, quantile, .025, na.rm = TRUE),
             hi = apply(B, 2, quantile, .975, na.rm = TRUE), dna = colMeans(Xd[i, , drop = FALSE]))
}))))
fwrite(prof, file.path(OUT, "cluster_profiles.csv"))
share <- dcast(rel[, .N, by = .(cls, group)], cls ~ group, value.var = "N", fill = 0)
share[, `:=`(cand_pct = round(100 * candidate / sum(candidate), 1), ctrl_pct = round(100 * control / sum(control), 1))]
cat("== cells per cluster (controls by nearest centroid)\n"); print(share)

prof[, cls_lab := factor(CLS_NAME[as.character(cls)], levels = CLS_NAME)]
prof[, lab := sprintf("%s (n = %d)", group, n)]
pA <- ggplot(prof, aes(r_um, cd45, colour = group, fill = group)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) +
  geom_line(linewidth = .9) + geom_point(size = 1) +
  facet_wrap(~ cls_lab, nrow = 1) +
  scale_colour_manual(values = GCOL, name = NULL) + scale_fill_manual(values = GCOL, name = NULL) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  labs(title = "Average CD45 profile per step-31 cluster, with the controls that share each shape",
       subtitle = "log2 of each 5 um band's mean CD45 over the cell's own 50-100 um mean; ribbon = 95% FOV cluster bootstrap. Controls assigned to the nearest cluster centroid",
       x = NULL, y = "log2 CD45 vs own 50-100 um") +
  theme_viz(10, TRUE) + theme(strip.text = element_text(colour = INK2, face = "bold"))
pB <- ggplot(prof, aes(r_um, dna, colour = group)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_line(linewidth = .9) + geom_point(size = 1) +
  facet_wrap(~ cls_lab, nrow = 1) +
  scale_colour_manual(values = GCOL, guide = "none") + scale_x_continuous(breaks = seq(0, 100, 25)) +
  labs(subtitle = "DNA, same normalisation", x = "Distance from centroid (um)", y = "log2 DNA vs own 50-100 um") +
  theme_viz(10, TRUE) + theme(strip.text = element_blank())
ntxt <- share[, sprintf("cluster %d: %d candidates (%.0f%%) / %d controls (%.0f%%)", cls, candidate, cand_pct, control, ctrl_pct)]
fig1 <- pA / pB + plot_layout(heights = c(2, 1)) +
  plot_annotation(caption = paste(c(paste(ntxt, collapse = "; "),
                                    "Clusters were defined on the candidates' profiles alone (step 31, k-means k = 5, a borderline pass of a single-Gaussian null). Source: output/33_cluster_followup/cluster_profiles.csv."),
                                  collapse = "\n"),
                  theme = theme(plot.caption = element_text(colour = MUTED, size = 7, hjust = 0), plot.background = element_rect(fill = SURF, colour = NA)))
sf("fig33_cluster_profiles.png", fig1, 14, 6.5)

## ---- 2. gene excess per cluster, step-28 framework --------------------------------
reg <- fread(cmd = "zcat output/16_regions/matched/cell_regions.csv.gz",
             select = c("cell_id", "samp", "fov_id", "is_defect", "tier_matched", "nCount_RNA"))[
  samp == "IB1" & tier_matched == TRUE & is_defect == FALSE]
hdr <- names(fread(cmd = sprintf("zcat %s/IB1_exprMat_file.csv.gz | head -1", FLAT)))
stopifnot(all(GENES %in% hdr))
ex <- fread(cmd = sprintf("zcat %s/IB1_exprMat_file.csv.gz | cut -d, -f%s%s", FLAT,
                          paste(match(c("fov", "cell_ID", GENES), hdr), collapse = ","),
                          if (SMOKE) " | head -300001" else ""))
ex[, cell_id := sprintf("c_2_%d_%d", fov, cell_ID)]          # c_2_ = IB1 (CLAUDE.md)
d <- merge(reg, ex[, c("cell_id", GENES), with = FALSE], by = "cell_id")
d[, chemo := CCL17 + CCL22]; d[, depth0 := nCount_RNA - chemo]; d <- d[depth0 > 0]
d[, depth_q := cut(depth0, unique(quantile(depth0, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
ref <- d[chemo < 10, c(list(dep = sum(as.numeric(depth0))), lapply(.SD, function(z) sum(as.numeric(z)))),
         by = depth_q, .SDcols = TEST]
for (g in TEST) ref[, (paste0("r_", g)) := get(g) / dep]
cand <- merge(d, rel[group == "candidate", .(cell_id, cls)], by = "cell_id")
cand <- merge(cand, ref[, c("depth_q", paste0("r_", TEST)), with = FALSE], by = "depth_q")
for (g in TEST) cand[, (paste0("e_", g)) := depth0 * get(paste0("r_", g))]
cat(sprintf("\ncandidates with counts: %d of %d\n", nrow(cand), rel[group == "candidate", .N]))
if (!SMOKE) stopifnot(nrow(cand) == rel[group == "candidate", .N])

obs_m <- as.matrix(cand[, ..TEST]); exp_m <- as.matrix(cand[, paste0("e_", TEST), with = FALSE])
lr_by <- function(rows, labels) {              # 5 x genes matrix of log excess
  t(sapply(1:5, function(k) {
    i <- rows[labels == k]
    log((colSums(obs_m[i, , drop = FALSE]) + 0.5) / (colSums(exp_m[i, , drop = FALSE]) + 0.5))
  }))
}
L0 <- lr_by(seq_len(nrow(cand)), cand$cls)
cf <- unique(cand$fov_id); cidx <- split(seq_len(nrow(cand)), cand$fov_id)
BB <- lapply(seq_len(N_BOOT), function(b) {
  rows <- unlist(cidx[as.character(sample(cf, replace = TRUE))], use.names = FALSE)
  lr_by(rows, cand$cls[rows])
})
lo <- Reduce(function(a, b) a, list(apply(simplify2array(BB), 1:2, quantile, .025, na.rm = TRUE)))
hi <- apply(simplify2array(BB), 1:2, quantile, .975, na.rm = TRUE)
n_k <- tabulate(cand$cls, 5)
het <- function(L) { m <- colSums(L * n_k) / sum(n_k); colSums(n_k * sweep(L, 2, m)^2) / sum(n_k) }
H0 <- het(L0)
HP <- replicate(N_PERM, het(lr_by(seq_len(nrow(cand)), sample(cand$cls))))
p_het <- (1 + rowSums(HP >= H0)) / (1 + N_PERM)
q_het <- p.adjust(p_het, "BH")
exg <- rbindlist(lapply(1:5, function(k) data.table(cls = k, n = n_k[k], gene = TEST, set = SET_OF[TEST],
  obs = colSums(obs_m[cand$cls == k, , drop = FALSE]), exp = round(colSums(exp_m[cand$cls == k, , drop = FALSE]), 2),
  excess = round(exp(L0[k, ]), 3), lo = round(exp(lo[k, ]), 3), hi = round(exp(hi[k, ]), 3))))
exg <- merge(exg, data.table(gene = TEST, het_stat = round(H0, 4), p_het = signif(p_het, 3), q_het = signif(q_het, 3)), by = "gene")
fwrite(exg[order(set, gene, cls)], file.path(OUT, "gene_excess_by_cluster.csv"))
## the gauge: each HRS/ambiguous gene's spread across clusters vs the T gauge's
gauge <- rowMeans(L0[, c("CD3E", "CD3D")])   # one value per CLUSTER
net <- rbindlist(lapply(1:5, function(k) data.table(cls = k, gene = c(SETS$HRS, SETS$ambiguous),
  net_vs_gauge = round(exp(L0[k, c(SETS$HRS, SETS$ambiguous)] - gauge[k]), 3))))
fwrite(net, file.path(OUT, "net_vs_T_gauge.csv"))
desc <- merge(cand[, .(n = .N, chemo_med = as.numeric(median(chemo)), depth0_med = as.numeric(median(depth0))), by = cls],
              rel[group == "candidate", .(hrs_score_med = round(median(hrs_score), 2), area_med = round(median(Area.um2), 1)), by = cls], by = "cls")[order(cls)]
fwrite(desc, file.path(OUT, "cluster_descriptives.csv"))
cat("\n== candidates per cluster\n"); print(desc)
cat("\n== excess by cluster (obs/exp vs depth-matched < 10 cells)\n")
w <- dcast(exg, set + gene + q_het ~ cls, value.var = "excess")
print(w[order(set, q_het)])
cat("\nheterogeneity across clusters, q < 0.1:\n"); print(exg[q_het < 0.1, .(gene, set, p_het, q_het)][!duplicated(gene)][order(q_het)])

## ---- 3. figure in the fig29 style: x = cluster -------------------------------------
exg[, panel := factor(fcase(set == "HRS", "HRS-leaning (pre-declared)",
                            set == "ambiguous", "HRS or activated APC (shared)",
                            set == "DC", "Dendritic-cell lineage",
                            default = "Macrophage / B / T (T = leakage gauge)"),
                      levels = c("HRS-leaning (pre-declared)", "HRS or activated APC (shared)",
                                 "Dendritic-cell lineage", "Macrophage / B / T (T = leakage gauge)"))]
exg[, col := fifelse(set %in% names(SET_COL), set, "other")]
SHOW <- c("IL13", "TNFRSF8", "LTA", "FUT4", "BATF3", "IRF4", "CD80", "CD83", "CD274", "PDCD1LG2",
          "LAMP3", "XCR1", "CLEC9A", "CD1C", "IDO1", "CCR7", "LYZ", "CD68", "CD79A", "CD3E", "CD3D")
exg[, shown := gene %in% SHOW]
exg[, glab := fifelse(q_het < 0.05, paste0(gene, " *"), gene)]
lab <- exg[cls == 5 & shown == TRUE]
p3 <- ggplot(exg, aes(cls, excess, group = gene)) +
  geom_hline(yintercept = 1, colour = AXIS, linewidth = .4) +
  geom_line(data = exg[shown == FALSE], colour = "#d8d6cf", linewidth = .35) +
  geom_errorbar(data = exg[shown == TRUE], aes(ymin = lo, ymax = hi, colour = col), width = .08, linewidth = .3, alpha = .5) +
  geom_line(data = exg[shown == TRUE], aes(colour = col), linewidth = .8) +
  geom_point(data = exg[shown == TRUE], aes(colour = col), size = 1.4) +
  ggrepel::geom_text_repel(data = lab, aes(label = glab, colour = col), size = 2.8, nudge_x = .35,
                           direction = "y", hjust = 0, segment.size = .2, min.segment.length = 0, seed = 1) +
  scale_y_continuous(trans = "log2", breaks = OR_BREAKS, labels = function(x) paste0(x, "x")) +
  scale_x_continuous(breaks = 1:5, limits = c(0.8, 6.4), labels = sub("\\n", " ", CLS_NAME)) +
  scale_colour_manual(values = SET_COL, guide = "none") +
  facet_wrap(~ panel, nrow = 1) +
  labs(title = "Gene excess by CD45-profile cluster, same framework as the brightness-tier figure",
       subtitle = sprintf("The %d candidate cells (CCL17+CCL22 >= 50) in their five step-31 clusters. Excess = observed / expected counts, expected from < 10 cells in the same depth decile", nrow(cand)),
       x = "Step-31 CD45 profile cluster", y = "Excess over depth-matched reference (log scale)",
       caption = paste("Bars: 95% FOV cluster bootstrap. * = heterogeneity across clusters q < 0.05 (n-weighted variance of log excess vs 1,000 permutations of the cluster labels, BH over genes).",
                       "\nCD3E / CD3D are the leakage gauge: a bright CD45 ring is T cells, whose transcripts leak into the cell they surround. Read an HRS-leaning difference only where it exceeds the gauge's.",
                       "\nSource: output/33_cluster_followup/gene_excess_by_cluster.csv.")) +
  theme_viz(10) + theme(strip.text = element_text(colour = INK2, face = "bold", hjust = 0),
                        axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5))
sf("fig33_gene_excess_by_cluster.png", p3, 14, 6)
cat("\ndone:", OUT, "\n")
