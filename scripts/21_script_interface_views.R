#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 21d -- figures for the FL/Hodgkin interface.
##
## Three panels-sets, one per question the step asks:
##   fig_21_axis_selection    why distance-to-hotspot was rejected and field level accepted
##   fig_21_interface_response the FL-specific field response, with its detection floor drawn
##   fig_21_arrangement        the two nulls, and the diagnostic that chooses between them
##
## COLOUR follows the project's fixed decisions: FL orange #eb6834, Hodgkin blue #0b3d7a,
## and the diverging scale runs blue -> red with red the UP end. All panels are within IB1,
## so the IB1/RB1 sample hues do not apply here.
##
## AXES: everything here is in um or cp10k, not canvas px and not slide mm, so THE AXIS RULE
## does not apply -- no scale reverse anywhere. (The rule: canvas px reverses both axes,
## slide mm reverses neither.)
##
## Usage: Rscript scripts/21_script_interface_views.R [tier]
## ---------------------------------------------------------------------------------------
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
setDTthreads(1)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
OUT  <- file.path("output/21_interface", TIER)
FIG  <- file.path(OUT, "figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

FL_COL <- "#eb6834"; HK_COL <- "#0b3d7a"; UP <- "#c0392b"; DN <- "#2a78d6"
INK <- "#1c1b19"; INK2 <- "#4a4844"; MUTED <- "#898781"; SURF <- "#fcfcfb"
theme_v <- theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#eceae4", linewidth = .3),
        axis.text = element_text(colour = MUTED, size = 8),
        axis.title = element_text(colour = INK2, size = 9),
        plot.title = element_text(colour = INK, face = "bold", size = 10),
        plot.subtitle = element_text(colour = INK2, size = 8.5),
        strip.text = element_text(colour = INK2, face = "bold", size = 9),
        legend.position = "bottom", legend.text = element_text(size = 8))
rd <- function(f) { p <- file.path(OUT, f); if (file.exists(p)) fread(p) else NULL }

GV  <- rd("gate_verdict.csv"); FLA <- rd("axisA_fl_fovs.csv")
AV  <- rd("arrangement_verdict.csv"); ND <- rd("arrangement_null_draws.csv")
ADJ <- rd("arrangement_null_adjudication.csv"); T1 <- rd("interface_tier1.csv")

## =========================== 1. axis selection =========================================
if (!is.null(FLA) && !is.null(GV)) {
  FLA[, expo := chemo_cp10k]
  pa <- ggplot(FLA, aes(d_hot_um, expo)) +
    geom_hline(yintercept = 11.38, linetype = "22", colour = MUTED, linewidth = .35) +
    geom_point(aes(size = fl_bins), colour = HK_COL, alpha = .75) +
    scale_size_area(max_size = 5, name = "FL bins in FOV") +
    scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
    annotate("label", x = max(FLA$d_hot_um), y = min(FLA$expo), hjust = 1, vjust = 0,
             size = 2.9, colour = INK, fill = SURF, linewidth = .25,
             label = sprintf("rho = %.3f, p = %.2f\nn = %d FL FOVs",
                             GV$axisA_rho_exposure_vs_dist, GV$axisA_p, nrow(FLA))) +
    labs(title = "A. Distance to a hotspot does not order field exposure",
         subtitle = "each point is one FL-carrying FOV; dashed line is the field cut (11.4 cp10k)",
         x = "distance to nearest hotspot bin (um, log)",
         y = "CCL17+CCL22 in the FOV's FL bins (cp10k, log)") + theme_v

  pc <- rbindlist(list(
    data.table(axis = "AXIS A\ndistance to hotspot", stratum = c("near","far"),
               val = c(GV$axisA_chemo_near, GV$axisA_chemo_far)),
    data.table(axis = "AXIS B\nown field level", stratum = c("low","high"),
               val = FLA[, {s <- median(expo)
                 c(sum(expo[expo<=s]*cells[expo<=s])/sum(cells[expo<=s]),
                   sum(expo[expo>s]*cells[expo>s])/sum(cells[expo>s]))}])))
  pc[, stratum := factor(stratum, c("near","far","low","high"))]
  pb <- ggplot(pc, aes(stratum, val, fill = axis)) +
    geom_col(width = .62) +
    geom_text(aes(label = sprintf("%.1f", val)), vjust = -0.4, size = 3, colour = INK2) +
    facet_wrap(~axis, scales = "free_x") +
    scale_fill_manual(values = c(MUTED, HK_COL), guide = "none") +
    expand_limits(y = max(pc$val)*1.18) +
    labs(title = "B. The positive control each axis has to pass",
         subtitle = "CCL17+CCL22 must be higher in the stratum the axis calls more exposed.\nAxis A fails (far > near); axis B passes by construction and is the axis used.",
         x = NULL, y = "CCL17+CCL22 (cp10k)") + theme_v
  ggsave(file.path(FIG, "fig_21_axis_selection.png"), pa + pb + plot_layout(widths = c(1, 1)),
         width = 10, height = 4.3, dpi = 150, bg = SURF)
  cat("wrote", file.path(FIG, "fig_21_axis_selection.png"), "\n")
}

## =========================== 2. the FL-specific field response ==========================
if (!is.null(T1) && nrow(T1)) {
  D <- T1[testable == TRUE & is.finite(delta)]
  D[, lab := factor(gene, levels = unique(gene[order(celltype, delta)]))]
  D[, sig := ifelse(is.finite(q_delta) & q_delta < 0.05 & concordant, "q < 0.05", "ns")]
  ## The MDE band is the honest half of any null: |delta| inside it could not have been
  ## resolved, so its absence is not evidence of absence.
  mde <- D[, .(m = median(mde_log2, na.rm = TRUE)), by = celltype]
  pr <- ggplot(D, aes(delta, lab, colour = sig)) +
    geom_rect(data = mde, inherit.aes = FALSE,
              aes(xmin = -m, xmax = m, ymin = -Inf, ymax = Inf),
              fill = "#eceae4", alpha = .55) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = .35) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0, linewidth = .45) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = c("ns" = MUTED, "q < 0.05" = UP), name = NULL) +
    facet_wrap(~celltype, scales = "free_y", nrow = 1) +
    labs(title = "The field response specific to FL tissue, within cell type",
         subtitle = paste0("delta = (slope on field exposure inside FL bins) - (the same slope in IB1 background bins).\n",
                           "Grey band is the median minimum detectable effect: a delta inside it could not have been resolved.\n",
                           "Bars are 95% intervals from a stratified FOV cluster bootstrap. Pre-declared gene sets only."),
         x = "delta (log2 expression per log2 field exposure)", y = NULL) +
    theme_v + theme(axis.text.y = element_text(size = 6.5))
  ggsave(file.path(FIG, "fig_21_interface_response.png"), pr,
         width = 11, height = 7.2, dpi = 150, bg = SURF)
  cat("wrote", file.path(FIG, "fig_21_interface_response.png"), "\n")
}

## =========================== 3. arrangement =============================================
if (!is.null(ND) && !is.null(AV)) {
  LV <- c("R_relocation","Rd_relocation_defect_matched","P_permutation")
  LB <- c("NULL R - uniform relocation\n(ignores the defect rim)",
          "NULL Rd - defect-matched\n(matches the hotspots' own tissue)",
          "NULL P - label permutation\n(destroys clustering; rejected)")
  ND <- ND[null %in% LV]
  ND[, null := factor(null, LV, LB)]
  ## obs$null must be the SAME FACTOR, not a character vector: a character in one layer
  ## coerces the combined facet variable and ggplot then re-sorts the panels
  ## alphabetically, which put the secondary null first.
  obs <- data.table(null = factor(LB, levels = LB), o = AV$obs_h2f_um)
  pl <- data.table(null = factor(LB, LB),
                   lab = c(sprintf("p = %.2f", AV$nullR_h2f_p),
                           sprintf("p = %.3f", AV$nullRd_h2f_p),
                           sprintf("p = %.3f (rejected)", AV$nullP_h2f_p)))
  p1 <- ggplot(ND, aes(h2f_um)) +
    geom_histogram(bins = 46, fill = "#d7d5cd", colour = NA) +
    geom_vline(data = obs, aes(xintercept = o), colour = FL_COL, linewidth = .8) +
    geom_text(data = pl, aes(x = Inf, y = Inf, label = lab), hjust = 1.08, vjust = 1.6,
              inherit.aes = FALSE, size = 3, colour = INK) +
    facet_wrap(~null, scales = "free_y", nrow = 1) +
    labs(title = "A. Hotspot to nearest FL focus, observed against three nulls",
         subtitle = sprintf(paste0("orange line = observed %d um. The two relocation nulls ",
                            "bracket the answer: matching the hotspots' own\ndefect burden ",
                            "moves the null median %d -> %d um and the observed value from ",
                            "ordinary to marginal."),
                            AV$obs_h2f_um, AV$nullR_h2f_med, AV$nullRd_h2f_med),
         x = "median hotspot -> nearest FL distance (um)", y = "null draws") + theme_v

  A <- copy(ADJ); A[, pattern := factor(pattern, rev(pattern))]
  p2 <- ggplot(A, aes(mnn_self_um, pattern,
                      fill = ifelse(grepl("OBSERVED", pattern), FL_COL, "#d7d5cd"))) +
    geom_col(width = .6) + scale_fill_identity() +
    geom_text(aes(label = sprintf("%d um  (%.2fx)", mnn_self_um, ratio_to_observed)),
              hjust = -0.06, size = 3, colour = INK2) +
    expand_limits(x = max(A$mnn_self_um)*1.45) +
    labs(title = "B. Which null holds the observed clustering fixed?",
         subtitle = sprintf(paste0("Mean nearest-neighbour spacing WITHIN the hotspot pattern. A valid null reproduces it.\n",
                           "Both relocation nulls do (1.03x, 1.02x) and are the admissible pair. NULL P inflates it 1.37x:\n",
                           "it spreads the hotspots out, which lengthens every distance, and that is where its nominal\n",
                           "p = %.3f comes from - not from arrangement."), AV$nullP_h2f_p),
         x = "within-pattern nearest-neighbour spacing (um)", y = NULL) + theme_v
  ggsave(file.path(FIG, "fig_21_arrangement.png"), p1 / p2 + plot_layout(heights = c(1, .85)),
         width = 12, height = 7.4, dpi = 150, bg = SURF)
  cat("wrote", file.path(FIG, "fig_21_arrangement.png"), "\n")
}
cat("OK\n")
