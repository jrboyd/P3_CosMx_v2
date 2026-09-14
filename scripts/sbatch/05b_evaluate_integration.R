# Step 5 evaluation: did integration buy confident B-cell identification in IB1
# without dissolving rare populations?
#
# Three questions, each with a metric that survives the fact that the variants
# produce different numbers of clusters:
#
# 1. MIXING -- per-cluster IB1 share; clusters (and cells) within 1.5x of the
#    whole-merge expectation. Unintegrated baseline: 0 of 7 at res 0.4.
#
# 2. RARE-POPULATION PRESERVATION -- the concern. Cluster-count-independent:
#    ARI/NMI between each clustering and the pre-existing InSituType / marker
#    lineage partitions, plus per-label BEST-F1 (the single cluster that best
#    recovers that label). A rare population that survives integration keeps a
#    high best-F1; one that is integrated away smears across clusters and drops.
#
# 3. PAYOFF vs SAFETY -- integration FORCES mixing, so co-clustering alone cannot
#    distinguish "IB1 B cells correctly aligned to RB1 B cells" from "IB1 junk
#    forcibly aligned". The discriminator: inside clusters where RB1 cells are
#    strongB-enriched, are the IB1 cells strongB-ENRICHED relative to IB1's own
#    3.9% baseline (real alignment), or at baseline / cluster-i dominated (junk)?
#
# Usage: sbatch --mem=200G scripts/sbatch/run_rscript.sh scripts/sbatch/05b_evaluate_integration.R

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table); setDTthreads(1)
})
OUT <- "output/05_integrate"; FIG <- file.path(OUT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
S_IB1 <- "#2a78d6"; S_RB1 <- "#eb6834"
GOOD <- "#0ca30c"; CRIT <- "#d03b3b"                     # status palette
SEQ <- c("#cde2fb","#9ec5f4","#6da7ec","#3987e5","#256abf","#184f95","#0d366b")

theme_viz <- function(base = 11, xgrid = FALSE) {
  th <- theme_minimal(base_size = base) +
    theme(plot.background = element_rect(fill = SURF, colour = NA),
          panel.background = element_rect(fill = SURF, colour = NA),
          panel.grid.major.y = element_line(colour = GRID, linewidth = 0.25),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.line.x = element_line(colour = AXIS, linewidth = 0.3), axis.ticks = element_blank(),
          axis.text = element_text(colour = MUTED, size = base - 2),
          axis.title = element_text(colour = INK2, size = base - 1),
          plot.title = element_text(colour = INK, face = "bold", size = base + 2),
          plot.subtitle = element_text(colour = INK2, size = base - 1),
          plot.caption = element_text(colour = MUTED, size = base - 3, hjust = 0),
          strip.text = element_text(colour = INK, face = "bold", size = base - 1),
          legend.text = element_text(colour = INK2, size = base - 2),
          legend.title = element_text(colour = INK2, size = base - 2),
          legend.position = "top", legend.justification = "left",
          legend.key.height = unit(9, "pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x = element_line(colour = GRID, linewidth = 0.25))
  th
}
sf <- function(n, p, w, h) { ggsave(file.path(FIG, n), p, width = w, height = h, dpi = 150, bg = SURF)
                             cat("wrote", file.path(FIG, n), "\n") }
sc <- function(n, d) { fwrite(d, file.path(FIG, n)); cat("wrote", file.path(FIG, n), "\n") }

## chance-corrected partition agreement, from the contingency table
ari_nmi <- function(a, b) {
  t <- table(a, b); n <- sum(t)
  ai <- rowSums(t); bj <- colSums(t)
  cmb <- function(x) sum(x * (x - 1) / 2)
  idx <- cmb(as.vector(t)); ea <- cmb(ai); eb <- cmb(bj); tot <- n * (n - 1) / 2
  exp_idx <- ea * eb / tot
  ari <- (idx - exp_idx) / ((ea + eb) / 2 - exp_idx)
  p <- t / n; pa <- ai / n; pb <- bj / n
  mi <- sum(ifelse(p > 0, p * log(p / outer(pa, pb)), 0))
  ha <- -sum(pa * log(pa)); hb <- -sum(pb * log(pb))
  c(ARI = ari, NMI = 2 * mi / (ha + hb))
}
## best single-cluster recovery of each label
best_f1 <- function(label, clust) {
  t <- table(label, clust); nl <- rowSums(t); nc <- colSums(t)
  f1 <- 2 * t / outer(nl, nc, "+")
  data.table(label = rownames(t), n_label = as.vector(nl),
             best_f1 = round(apply(f1, 1, max), 3),
             best_cluster = colnames(t)[apply(f1, 1, which.max)])
}

md <- as.data.table(readRDS(file.path(OUT, "merged_integrated.Rds"))@meta.data,
                    keep.rownames = "cell")
md[, samp := factor(merge_sample, levels = c("IB1", "RB1"))]
md[, is_i := inSitu_clusters == "i"]
EXP <- md[, mean(samp == "IB1")]
BASE_SB <- md[, .(base = mean(strongB, na.rm = TRUE)), by = samp]
cat("cells:", nrow(md), "| expected IB1 share:", round(100*EXP, 2), "%\n")
print(BASE_SB)

VARIANTS <- c("none (res 0.4)" = "pearson_clusters_low",
              "none (res 0.8)" = "pearson_clusters_high",
              "RPCA (res 0.4)" = "int_rpca_low_full",
              "RPCA (res 0.8)" = "int_rpca_high_full",
              "CCA (res 0.4)"  = "int_cca_low_full",
              "CCA (res 0.8)"  = "int_cca_high_full")
VARIANTS <- VARIANTS[VARIANTS %in% colnames(md)]
cat("variants found:", paste(names(VARIANTS), collapse = " | "), "\n")
stopifnot(length(VARIANTS) >= 3)
md[, (VARIANTS) := lapply(.SD, as.character), .SDcols = VARIANTS]

## ---- 1. MIXING -------------------------------------------------------------
comp <- rbindlist(lapply(names(VARIANTS), function(v) {
  g <- VARIANTS[[v]]
  d <- md[, .(n = .N, n_IB1 = sum(samp == "IB1"),
              pct_IB1 = round(100*mean(samp == "IB1"), 2),
              pct_strongB_IB1 = round(100*mean(strongB[samp=="IB1"], na.rm=TRUE), 2),
              pct_strongB_RB1 = round(100*mean(strongB[samp=="RB1"], na.rm=TRUE), 2),
              pct_i_IB1 = round(100*mean(is_i[samp=="IB1"]), 2)),
          by = .(cluster = get(g))]
  d[, `:=`(variant = v, log2_enrich_IB1 = round(log2((pct_IB1/100 + 1e-6)/EXP), 3))]
  d[, within := abs(log2_enrich_IB1) < log2(1.5)][]
}), use.names = TRUE)
sc("05_cluster_composition_by_method.csv", comp[order(variant, -pct_IB1)])
mix <- comp[, .(clusters = .N, within_1.5x = sum(within),
                pct_cells_within = round(100*sum(n[within])/sum(n), 1),
                median_abs_log2 = round(median(rep(abs(log2_enrich_IB1), n)), 3)),
            by = variant]
cat("\n=== MIXING ===\n"); print(mix)
sc("05_mixing_summary.csv", mix)

## ---- 2. RARE-POPULATION PRESERVATION ---------------------------------------
agree <- rbindlist(lapply(names(VARIANTS), function(v) {
  g <- md[[VARIANTS[[v]]]]
  rbindlist(lapply(c("inSitu_clusters", "lineage"), function(lb) {
    a <- ari_nmi(md[[lb]], g)
    data.table(variant = v, reference = lb, n_clusters = uniqueN(g),
               ARI = round(a[["ARI"]], 4), NMI = round(a[["NMI"]], 4))
  }))
}))
cat("\n=== PARTITION AGREEMENT (chance-corrected) ===\n"); print(agree)
sc("05_partition_agreement.csv", agree)

rare <- rbindlist(lapply(names(VARIANTS), function(v) {
  g <- md[[VARIANTS[[v]]]]
  rbindlist(lapply(c("inSitu_clusters", "lineage"), function(lb) {
    b <- best_f1(md[[lb]], g); b[, `:=`(variant = v, reference = lb)]
    b[, pct_of_cells := round(100*n_label/nrow(md), 2)][]
  }))
}))
rare[, is_rare := pct_of_cells < 5]
sc("05_rare_preservation.csv", rare[order(reference, label, variant)])
cat("\n=== RARE LABELS (<5% of cells): best-F1 by variant ===\n")
print(dcast(rare[is_rare == TRUE], reference + label + pct_of_cells ~ variant,
            value.var = "best_f1"))

## ---- 3. PAYOFF vs SAFETY ---------------------------------------------------
bb <- BASE_SB[samp == "IB1", base]; rb <- BASE_SB[samp == "RB1", base]
pay <- comp[, {
  bcl <- pct_strongB_RB1 >= 200*rb          # RB1 strongB >= 2x RB1 baseline (in %)
  .(rb1_B_clusters = sum(bcl),
    IB1_cells_in_them = sum(n_IB1[bcl]),
    IB1_strongB_pct = round(weighted.mean(pct_strongB_IB1[bcl], n_IB1[bcl]), 2),
    IB1_strongB_enrich = round(weighted.mean(pct_strongB_IB1[bcl], n_IB1[bcl])/(100*bb), 2),
    IB1_cluster_i_pct = round(weighted.mean(pct_i_IB1[bcl], n_IB1[bcl]), 2))
}, by = variant]
pay[, IB1_baseline_strongB_pct := round(100*bb, 2)]
pay[, IB1_baseline_cluster_i_pct := round(100*md[samp == "IB1", mean(is_i)], 2)]
cat("\n=== PAYOFF: IB1 cells inside RB1-strongB-enriched clusters ===\n"); print(pay)
sc("05_bcell_recovery.csv", pay)

## ---- figures ---------------------------------------------------------------
comp[, variant := factor(variant, levels = names(VARIANTS))]
comp[, clu := factor(paste(variant, cluster),
       levels = comp[order(variant, pct_IB1), paste(variant, cluster)])]
p1 <- ggplot(comp, aes(pct_IB1, clu, colour = ifelse(log2_enrich_IB1 > 0, "IB1", "RB1"))) +
  geom_segment(aes(x = 100*EXP, xend = pct_IB1, yend = clu), linewidth = 0.4) +
  geom_point(aes(size = n)) +
  geom_vline(xintercept = 100*EXP, colour = INK, linetype = "22", linewidth = 0.4) +
  facet_grid(variant ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(IB1 = S_IB1, RB1 = S_RB1), name = "enriched in") +
  scale_size_area(max_size = 6, labels = function(x) format(x, big.mark = ","), name = "cells") +
  scale_x_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  scale_y_discrete(labels = function(v) sub("^.*\\) ", "", v)) +
  labs(title = "Does integration make the samples co-cluster?",
       subtitle = sprintf("IB1 share per cluster; dashed line = whole-merge expectation (%.1f%%).", 100*EXP),
       x = "IB1 share of cluster", y = "cluster",
       caption = "Data: 05_cluster_composition_by_method.csv / 05_mixing_summary.csv") +
  theme_viz(xgrid = TRUE) +
  theme(panel.grid.major.y = element_blank(), strip.text.y = element_text(angle = 0))
sf("fig_05_mixing_by_method.png", p1, 9, 11)

rr <- rare[is_rare == TRUE]
rr[, variant := factor(variant, levels = names(VARIANTS))]
p2 <- ggplot(rr, aes(best_f1, reorder(label, best_f1))) +
  geom_col(fill = SEQ[4], width = 0.7) +
  facet_grid(reference ~ variant, scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Are rare populations integrated away?",
       subtitle = paste0("Best-F1 recovery of each rare label (<5% of cells) by any single cluster.\n",
                         "A population that survives keeps a high best-F1; one that smears drops toward 0."),
       x = "best-F1 by a single cluster", y = NULL,
       caption = "Data: 05_rare_preservation.csv") +
  theme_viz(xgrid = TRUE) + theme(panel.grid.major.y = element_blank())
sf("fig_05_rare_population_preservation.png", p2, 12, 8)

payl <- melt(pay[, .(variant, `IB1 strongB in RB1-B clusters` = IB1_strongB_pct,
                     `IB1 cluster-i in RB1-B clusters` = IB1_cluster_i_pct)],
             id.vars = "variant", variable.name = "measure", value.name = "pct")
bl <- data.table(measure = factor(c("IB1 strongB in RB1-B clusters",
                                    "IB1 cluster-i in RB1-B clusters"),
                 levels = levels(payl$measure)),
                 y = c(100*bb, 100*md[samp == "IB1", mean(is_i)]))
payl[, variant := factor(variant, levels = names(VARIANTS))]
p3 <- ggplot(payl, aes(variant, pct)) +
  geom_col(fill = SEQ[4], width = 0.68) +
  geom_hline(data = bl, aes(yintercept = y), colour = CRIT, linetype = "22", linewidth = 0.45) +
  facet_wrap(~measure, ncol = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(c(0, 0.08)), labels = function(x) paste0(x, "%")) +
  labs(title = "Can we call B cells in IB1 by co-clustering with RB1?",
       subtitle = paste0("Inside clusters where RB1 is strongB-enriched (>=2x its baseline): what are the ",
         "IB1 cells there?\nDashed line = IB1's whole-tissue baseline. Above it on strongB = real ",
         "alignment; at/below, or high cluster-i = forced mixing."),
       x = NULL, y = NULL, caption = "Data: 05_bcell_recovery.csv") +
  theme_viz() + theme(axis.text.x = element_text(angle = 20, hjust = 1, colour = INK2))
sf("fig_05_bcell_recovery_IB1.png", p3, 9, 7)
cat("\nOK\n")
