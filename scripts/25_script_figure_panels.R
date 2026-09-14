#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 25b -- the manuscript figure panels tagged [make] in reports/manuscript_plan.html.
##
## WHAT THIS IS. Plotting, not analysis. Every number is read from a saved result table and
## nothing is recomputed here, so a panel cannot drift from the result it depicts. Each panel
## is written as its own PNG named for its position in the review's figure list, so the
## figures can be assembled without re-running anything.
##
## CONVENTIONS OBSERVED (CLAUDE.md "Colour and channel conventions"):
##   - the two samples have fixed hues: IB1 #2a78d6, RB1 #eb6834
##   - diverging scales run blue -> red with RED the UP/positive end
##   - the disease hues (FL orange, Hodgkin blue) are reserved; where a sample hue and a
##     disease hue coincide (RB1 and FL are both #eb6834) the panels here are per-sample or
##     per-disease, never both at once, so there is no clash inside a panel
##   - `Hodgkin chemokine niche` is a REGION, not a lineage, and any panel showing it must
##     render it distinctly and say so
##
## NO AXIS-REVERSAL ANYWHERE. None of these panels is in canvas-px or slide-mm space -- they
## are all statistical summaries -- so THE AXIS RULE does not apply. Spatial panels live in
## their own steps' view scripts.
##
## Usage: Rscript scripts/25_script_figure_panels.R
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
setDTthreads(1); options(width = 200)

OUT <- "output/25_figures"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
IF21 <- "output/21_interface/matched"; PD17 <- "output/17_pd1/matched"
MISSING <- character()
rd <- function(p) { if (!file.exists(p)) { MISSING <<- c(MISSING, p); return(NULL) }; fread(p) }
sv <- function(nm, p, w, h) { f <- file.path(OUT, paste0("fig_25_", nm, ".png"))
  ggsave(f, p, width = w, height = h, dpi = 200, bg = "white"); cat("wrote", f, "\n") }

SAMP <- c(IB1 = "#2a78d6", RB1 = "#eb6834")
## An explicit log2 axis. ggplot's default breaks on trans="log2" rendered as
## "7.6293954.8628123.045..." on the odds-ratio forests -- unreadable and run together. Always
## give breaks for a log axis rather than trusting the default.
OR_BREAKS <- c(0.25, 0.5, 1, 2, 4, 8)
scale_or <- function(lim = NULL) scale_x_continuous(
  trans = "log2", breaks = OR_BREAKS, labels = function(x) format(x, drop0trailing = TRUE),
  limits = lim)
UP <- "#c0392b"; DN <- "#2a78d6"; HODG <- "#0b3d7a"; FLO <- "#eb6834"
INK <- "#1c1b19"; INK2 <- "#4a4844"; MUTED <- "#898781"; SURF <- "#fcfcfb"
th <- theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#eceae4", linewidth = .3),
        axis.text = element_text(colour = INK2, size = 8),
        axis.title = element_text(colour = INK2, size = 9),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = INK2, size = 8.2),
        plot.caption = element_text(colour = MUTED, size = 7.2, hjust = 0),
        strip.text = element_text(colour = INK2, face = "bold", size = 8.5),
        legend.position = "bottom", legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

## ===== Fig 1d -- composition with FOV-level spread (closes G8) =========================
C <- rd(file.path(OUT, "composition_fov_summary.csv"))
if (!is.null(C)) {
  C[, lab := ifelse(is_region, paste0(celltype_cluster, " *"), celltype_cluster)]
  ord <- C[samp == "IB1"][order(share_pooled), lab]
  C[, lab := factor(lab, levels = unique(c(ord, setdiff(lab, ord))))]
  p <- ggplot(C, aes(100*share_pooled, lab, colour = samp)) +
    geom_linerange(aes(xmin = 100*share_q1, xmax = 100*share_q3), linewidth = 2.4, alpha = .30) +
    geom_linerange(aes(xmin = 100*share_p10, xmax = 100*share_p90), linewidth = .5, alpha = .55) +
    geom_point(aes(x = 100*share_med), shape = 124, size = 3.4) +
    geom_point(size = 2.3) +
    facet_wrap(~samp, nrow = 1) +
    scale_colour_manual(values = SAMP, guide = "none") +
    labs(title = "Composition, with the spread across fields of view",
         subtitle = paste0("Dot = pooled share; tick = FOV median; thick bar = IQR; thin bar = 10th-90th percentile.\n",
           "High-confidence cells only. * marks the Hodgkin chemokine niche, which is a REGION label, not a lineage."),
         x = "share of high-confidence labelled cells in the FOV (%)", y = NULL,
         caption = paste0("The spread is heterogeneity across the section, not a confidence interval: n = 1 patient, so no population mean exists.\n",
           "Note IB1's myeloid and B IQRs overlap almost completely, so its two largest fractions are not separable across the section;\n",
           "RB1's B-cell fraction is by contrast tight. Source: output/25_figures/composition_fov_summary.csv")) + th
  sv("1d_composition_fov", p, 9.6, 4.6)
}

## ===== Fig 2d -- caller agreement by lineage ===========================================
A <- rd("output/24_arch_relabel/matched/label_agreement.csv")
if (!is.null(A)) {
  AG <- A[, .(cells = sum(N), agree = sum(N[lab16 == lab19])), by = lab19][
          , pct := 100*agree/cells][order(pct)]
  AG[, lab19 := factor(lab19, levels = lab19)]
  p <- ggplot(AG, aes(pct, lab19)) +
    geom_vline(xintercept = 70.2, linetype = "22", colour = MUTED, linewidth = .4) +
    geom_segment(aes(x = 0, xend = pct, yend = lab19), colour = "#d7d5cd", linewidth = 1.6) +
    geom_point(aes(size = cells, colour = pct < 50)) +
    scale_colour_manual(values = c(`FALSE` = HODG, `TRUE` = UP), guide = "none") +
    scale_size_area(max_size = 6, name = "cells named by both", labels = scales::comma) +
    annotate("text", x = 70.2, y = .6, label = "overall 70.2%", hjust = -.06, size = 2.7,
             colour = MUTED) +
    scale_x_continuous(limits = c(0, 100)) +
    labs(title = "The two cell-type callers agree on 70.2% of cells",
         subtitle = "Agreement between the step-16 argmax caller and the step-19 cluster caller, on the 998,544 cells both name.",
         x = "agreement (%)", y = NULL,
         caption = paste0("Red marks lineages where the callers agree less than half the time. NK and plasma cell are effectively a coin flip,\n",
           "so no claim should rest on either identity from one caller alone. Source: output/24_arch_relabel/matched/label_agreement.csv")) + th
  sv("2d_caller_agreement", p, 7.4, 4.4)
}

## ===== Fig 4b -- BATF3 against distance to the nearest cDC =============================
B <- rd("output/16_checkpoint/batf3_spillover.csv")
if (!is.null(B)) {
  M <- melt(B[, .(dband, cp10k = ratio_cp10k, detection = ratio_det)],
            id.vars = "dband", variable.name = "measure", value.name = "ratio")
  M[, dband := factor(dband, levels = B$dband)]
  p <- ggplot(M, aes(dband, ratio, colour = measure, group = measure)) +
    geom_hline(yintercept = 1, colour = MUTED, linewidth = .4) +
    geom_line(linewidth = .8) + geom_point(size = 2.6) +
    scale_colour_manual(values = c(cp10k = HODG, detection = "#1baf7a"), name = NULL) +
    labs(title = "BATF3 enrichment is HIGHEST farthest from any cDC",
         subtitle = "Hotspot : background ratio for BATF3, by distance from the nearest cDC. Spillover predicts the opposite.",
         x = "distance to nearest cDC (um)", y = "hotspot : background ratio",
         caption = paste0("A BATF3-high dendritic cell bleeding transcripts into its neighbours would raise BATF3 most CLOSE to a cDC and decay away.\n",
           "The ratio rises instead, to 5.14x beyond 120 um, so this is a regional field and not a segmentation artefact.\n",
           "Source: output/16_checkpoint/batf3_spillover.csv")) + th
  sv("4b_batf3_vs_cdc_distance", p, 7.2, 4.4)
}

## ===== Fig 4c -- the CCR4 gradient, and what does NOT come with it =====================
G <- rd("output/16_gradient/matched/ccr4_axis.csv")
if (!is.null(G)) {
  K <- G[gene %in% c("CCR4","FOXP3","IL2RA","CTLA4","PDCD1","IKZF2")]
  K[, set := factor(fifelse(grepl("^FL", set), "inside FL foci", "all non-hotspot bins"),
                    levels = c("all non-hotspot bins", "inside FL foci"))]
  K[, gene := factor(gene, levels = c("CCR4","CTLA4","PDCD1","IKZF2","FOXP3","IL2RA"))]
  K[, sig := q_field < 0.05]
  p <- ggplot(K, aes(rho_field, gene, colour = set, shape = sig)) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = .4) +
    geom_point(size = 3, position = position_dodge(width = .55)) +
    scale_colour_manual(values = c(`all non-hotspot bins` = MUTED, `inside FL foci` = FLO),
                        name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "q < 0.05") +
    labs(title = "A CCR4 gradient without a Treg signature",
         subtitle = "Spearman rho of each gene against the CCL17+CCL22 field level, per 100 um bin.",
         x = "rho with field level", y = NULL,
         caption = paste0("CCR4 -- the receptor for CCL17 and CCL22 -- was declared as the readout BEFORE any genome-wide scan, so it is not a post-hoc pick.\n",
           "It is ~3.5x more strongly field-associated inside FL foci (0.415) than tissue-wide (0.120), while FOXP3 is flat and IL2RA runs negative.\n",
           "This is chemoattraction of CCR4+ cells or CCR4 upregulation, NOT Treg recruitment. Source: 16_gradient/matched/ccr4_axis.csv")) + th
  sv("4c_ccr4_axis", p, 7.6, 4.4)
}

## ===== Fig 5a -- module contrasts, with the detectable effect drawn ====================
MC <- rd(file.path(PD17, "module_contrasts.csv")); PB <- rd(file.path(PD17, "power_bounds.csv"))
if (!is.null(MC) && !is.null(PB)) {
  X <- MC[unit == "FOVxregion"][PB[, .(contrast, module, mde_log2)], on = .(contrast, module),
                                nomatch = 0L]
  X[, sig := q < 0.05 & concordant]
  ## Two of these contrasts are CROSS-SAMPLE, and CLAUDE.md fact 5 forecloses exactly those
  ## for the offset-affected genes: the between-slide reference-zone difference reaches OR
  ## 5.33, larger than any within-slide regional effect. Rendering them beside the
  ## within-sample contrasts without saying so invites a reader to take a red point in the
  ## "IB1 FL vs RB1 foll" column at face value. They are kept (dropping them would hide that
  ## they were computed) and labelled.
  X[, cross := grepl("IB1", contrast) & grepl("RB1", contrast)]
  X[, contrast := factor(ifelse(cross, paste0("[cross-sample] ", contrast), contrast))]
  X[, contrast := factor(contrast, levels = c(sort(unique(contrast[!cross])),
                                              sort(unique(contrast[cross]))))]
  X[, module := factor(module, levels = X[contrast == levels(contrast)[1]][order(log2fc), module])]
  p <- ggplot(X, aes(log2fc, module)) +
    geom_rect(aes(xmin = -mde_log2, xmax = mde_log2, ymin = -Inf, ymax = Inf),
              fill = "#eceae4", alpha = .5) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = .4) +
    geom_point(aes(colour = sig, shape = cross), size = 2.6) +
    scale_colour_manual(values = c(`FALSE` = MUTED, `TRUE` = UP),
                        labels = c("not significant", "q < 0.05, concordant"), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                       labels = c("within sample", "CROSS-SAMPLE: see the offset bound"),
                       name = NULL) +
    facet_wrap(~contrast, nrow = 1) +
    labs(title = "Module scores by region, with each contrast's detection limit",
         subtitle = "Grey band is the minimum detectable effect: a value inside it could not have been resolved.",
         x = "log2 fold change vs the same sample's reference region", y = NULL,
         caption = paste0("IB1's hotspots carry the PD-1 axis and the 9p24.1 module while the interferon-gamma response does NOT rise -- ligand-side\n",
           "without interferon. IB1's FL foci carry neither, and their MDE of 1.14-fold means an effect of the hotspot's size would have been seen.\n",
           "TRIANGLES ARE CROSS-SAMPLE and must not be read at face value: the between-slide offset on these genes reaches OR 5.33 in the two\n",
           "samples' REFERENCE zones, larger than any within-slide regional effect, so those columns are bounded by the design (CLAUDE.md fact 5).\n",
           "Sources: 17_pd1/matched/{module_contrasts,power_bounds}.csv")) + th
  sv("5a_module_contrasts", p, 10.4, 4.6)
}

## ===== Fig 5b/5c -- per-cell binomial odds ratios ======================================
BM <- rd(file.path(PD17, "binomial_model.csv"))
forest <- function(ctr, ttl, sub, cap, both_ref) {
  D <- BM[contrast == ctr & (if (both_ref) TRUE else refdef == "far_background")]
  if (!nrow(D)) return(NULL)
  D[, stratum := factor(stratum, levels = c("all","B","T","myeloid"))]
  D[, gene := factor(gene, levels = D[stratum == "all" & refdef == D$refdef[1]][order(or), gene])]
  D[, sig := q < 0.05]
  p <- ggplot(D, aes(or, gene, colour = refdef, shape = sig)) +
    geom_vline(xintercept = 1, colour = MUTED, linewidth = .4) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .45,
                   position = position_dodge(width = .6)) +
    geom_point(size = 1.9, position = position_dodge(width = .6)) +
    scale_or() +
    ## Clip rather than let one very wide interval (PDCD1 in B) compress every other row.
    ## Limits are given SORTED, so this only clips -- it never flips an axis.
    coord_cartesian(xlim = c(0.28, 6)) +
    scale_colour_manual(values = c(all_background = MUTED, far_background = HODG),
                        labels = c("all background", "far background (>800um)"), name = "reference") +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                       labels = c("not significant", "q < 0.05"), name = NULL) +
    facet_wrap(~stratum, nrow = 1) +
    labs(title = ttl, subtitle = sub,
         x = "odds ratio for positivity (log2 scale; intervals beyond the axis are clipped)",
         y = NULL, caption = cap) + th
  p
}
p <- forest("cHL_hotspot vs background",
  "IB1's Hodgkin hotspots: the ligands rise, the receptor does not",
  "Per-cell binomial model, positivity >= 2 counts, region and depth decile as covariates; intervals resample FOVs.",
  paste0("Both reference definitions are shown because 40% of IB1's `background` sits within 800 um of a hotspot and carries field signal.\n",
    "Cleaning the reference STRENGTHENS every ligand result (CD274 1.74 -> 1.83, PDCD1LG2 1.81 -> 1.92), so the asymmetry is not an artefact\n",
    "of where the reference was drawn. Source: 17_pd1/matched/binomial_model.csv"), TRUE)
if (!is.null(p)) sv("5b_forest_ib1_hotspot", p, 10.8, 4.8)

p <- forest("FL_follicular vs interfollicular",
  "RB1's follicles: the receptor rises, in the expected cells",
  "Per-cell binomial model, RB1 follicular vs interfollicular. Same construction as panel b.",
  paste0("PDCD1 in T cells is 1.82x more likely inside follicles (q = 4.6e-14) with CTLA4, CXCL13, TOX and TIGIT alongside -- the recognised\n",
    "follicular PD-1 pattern of FL, recovered without being sought, since the Tfh caller deliberately excludes PDCD1. CD274 falls in follicles.\n",
    "Source: 17_pd1/matched/binomial_model.csv"), FALSE)
if (!is.null(p)) sv("5c_forest_rb1_follicular", p, 10.8, 4.8)

## ===== Fig 5d -- the difference-in-differences, beside its reference contrast ==========
DD <- rd(file.path(PD17, "did_within_celltype.csv"))
if (!is.null(DD)) {
  M <- melt(DD[, .(stratum, gene, ref_class,
                   `disease contrast` = or_focus_raw, `reference contrast` = or_reference,
                   `difference-in-differences` = or_did)],
            id.vars = c("stratum","gene","ref_class"), variable.name = "which", value.name = "or")
  M[, gene := factor(gene, levels = DD[stratum == "T"][order(or_focus_raw), gene])]
  M[, offset := grepl("LARGE", ref_class)]
  p <- ggplot(M, aes(or, gene, colour = which, shape = offset)) +
    geom_vline(xintercept = 1, colour = MUTED, linewidth = .4) +
    geom_point(size = 2.2, position = position_dodge(width = .6)) +
    scale_or() + coord_cartesian(xlim = c(0.3, 6)) +
    scale_colour_manual(values = c(`disease contrast` = UP, `reference contrast` = MUTED,
                                   `difference-in-differences` = HODG), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                       labels = c("reference NULL: power-limited", "reference LARGE: offset"),
                       name = NULL) +
    facet_wrap(~stratum, nrow = 1) +
    labs(title = "Cross-sample contrasts, shown against the contrast between the two sections' reference zones",
         subtitle = "Where the reference contrast is as large as the disease contrast, the apparent difference is the section, not the biology.",
         x = "odds ratio, IB1 vs RB1 (log2 scale)", y = NULL,
         caption = paste0("CXCL13 deflates from 5.21 to 0.98 and LAG3 from 4.64 to 1.45 -- retractions. But where the reference contrast is NULL (open symbols)\n",
           "a null difference-in-differences is a power failure, not a retraction: its detectable effect is only 1.8-2.5 fold. CD274 is the live exception,\n",
           "with a null reference and a residual 1.76-fold. Source: 17_pd1/matched/did_within_celltype.csv")) + th
  sv("5d_did_vs_reference", p, 11.0, 5.0)
}

## ===== Fig 5e -- 9p24.1 dosage, and why it is not dosage ===============================
## Two panels rather than one. A first version printed the cytoband percentile as text inside
## the score panel at x = Inf; it landed on top of the fibroblast points. The percentile is a
## different quantity on a different scale, so it gets its own panel.
##
## COLOUR: deliberately NOT red/blue here. Red-to-blue is this project's DIVERGING convention
## (red = up), and these two series are categorical gene subsets, not directions -- rendering
## them red and blue invites the reader to read a direction that is not there. Grey is the
## circular subset (immune genes, which interferon alone can raise) and purple the
## load-bearing one (genes no cytokine story explains).
DS <- rd("output/23_dosage/matched/dosage_9p24.csv")
if (!is.null(DS)) {
  IMM <- "#8f8b84"; CLEAN <- "#7b4ea3"
  DS[, ctr := factor(sub(" vs .*$", "", contrast),
                     levels = c("IB1 hotspot","IB1 FL focus","RB1 follicular"))]
  DS[, celltype := factor(celltype, levels = c("B","fibro","myeloid","T"))]
  M <- melt(DS[, .(ctr, celltype, q_clean,
                   `immune / IFN-inducible (3 genes)` = score_immune,
                   `not immune-functional (5 genes)` = score_clean)],
            id.vars = c("ctr","celltype","q_clean"), variable.name = "subset",
            value.name = "score")
  M[, sig := q_clean < 0.05 & grepl("^not", subset)]
  pA <- ggplot(M, aes(score, celltype, colour = subset, shape = sig)) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = .4) +
    geom_point(size = 3, position = position_dodge(width = .55)) +
    scale_colour_manual(values = c(IMM, CLEAN), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
    facet_wrap(~ctr, nrow = 1) +
    labs(title = "9p24.1 rises in IB1's hotspots -- but this is not evidence of copy gain",
         subtitle = paste0("Mean log2 fold change of the band's genes vs the same sample's reference region. ",
                           "Filled = q < 0.05 on the non-immune subset."),
         x = "mean log2 fold change of the 9p24.1 gene set", y = NULL) + th
  pB <- ggplot(unique(DS[, .(ctr, celltype, pct_of_bands)]),
               aes(pct_of_bands, celltype)) +
    geom_vline(xintercept = c(50, 95), linetype = "22", colour = MUTED, linewidth = .35) +
    geom_segment(aes(x = 0, xend = pct_of_bands, yend = celltype),
                 colour = "#d7d5cd", linewidth = 1.4) +
    geom_point(aes(colour = pct_of_bands >= 95), size = 2.6) +
    scale_colour_manual(values = c(`FALSE` = MUTED, `TRUE` = CLEAN), guide = "none") +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 95)) +
    facet_wrap(~ctr, nrow = 1) +
    labs(subtitle = "Percentile of the 9p24.1 band among the 440+ cytobands measured on the same panel, in the same cells",
         x = "percentile among all cytobands", y = NULL,
         caption = paste0("Copy gain is indifferent to gene function and raises the whole interval; interferon raises only CD274/PDCD1LG2/JAK2. The non-immune\n",
           "subset does rise in IB1's T and fibroblast compartments, at the 99th percentile of all cytobands. But the SAME statistic is significant in RB1,\n",
           "which has no Hodgkin component, and the signal spans three of four lineages when a clonal event must be confined to one.\n",
           "Source: output/23_dosage/matched/dosage_9p24.csv")) + th
  sv("5e_9p24_dosage", pA / pB + patchwork::plot_layout(heights = c(1, .85)), 10.0, 6.4)
}

## ===== Fig 6b -- architecture under both callers, at matched n =========================
CE <- rd("output/24_arch_relabel/matched/control_effect_size.csv")
if (!is.null(CE)) {
  CE[, cond := factor(condition, c("A_full","B_matched"),
                      c("all labelled cells", "matched cell counts"))]
  CE[, caller := factor(labels, c("lab16","lab19"),
                        c("step-16 argmax caller", "step-19 cluster caller"))]
  p <- ggplot(CE, aes(med_abs_l2, caller, colour = samp)) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .6) +
    geom_point(size = 3) +
    scale_colour_manual(values = SAMP, name = NULL) +
    facet_grid(samp ~ cond, switch = "y") +
    labs(title = "The better labels see more structure, and it is not a power artefact",
         subtitle = paste0("Median |log2 enrichment| over all cell-type pairs: how far the arrangement sits from chance.\n",
           "Right column down-samples every lineage in every region to min(n16, n19), so power is equalised by construction."),
         x = "median |log2 enrichment| from a composition- and density-preserving permutation", y = NULL,
         caption = paste0("Matching does not close the gap. RB1 goes 0.059 -> 0.128 between callers at matched n, i.e. the step-19 labels put its enrichments\n",
           "2.15x further from chance; the same pair is further from chance under them in 80.3% of comparisons (p = 2.9e-39). The earlier reading that\n",
           "RB1 'sits at chance' measured the caller's silence, not the tissue. Both nodes are non-random; IB1 about twice as strongly.\n",
           "Source: output/24_arch_relabel/matched/control_effect_size.csv")) +
    th + theme(legend.position = "none")
  sv("6b_architecture_both_callers", p, 8.6, 4.6)
}

## ===== Fig 7d -- the distance-band profile that separates biology from leakage =========
## USES enrich_fov, NOT the analytic `enrich` in lr_proximity_profile.csv. The two are
## different quantities -- the analytic one includes cross-FOV cell pairs, enrich_fov uses
## within-FOV pairs only and is what the bootstrap resamples and what step 22 reports. A first
## version of this panel plotted the analytic column and disagreed with the review's own
## numbers: PDCD1LG2-PDCD1 came out NA / 1.65 / 1.13 / 1.03 instead of 2.28 / 1.54 / 1.14 /
## 1.01, dropping the 0-10 band entirely to the expected-count floor and losing the very pair
## the panel exists to show.
B7 <- rd("output/22_lr/matched/lr_proximity_tier1.csv")
TG <- rd("output/22_lr/matched/snapshots/snapshot_targets.csv")
if (!is.null(B7) && !is.null(TG)) {
  K <- merge(B7, TG[, .(samp, bin_region, ligand, receptor, status)],
             by = c("samp","bin_region","ligand","receptor"))
  K[, pair := paste0(ligand, "-", receptor)]
  K[, band := factor(band, levels = c("0-10","10-20","20-40","40-80"))]
  K[, lab := fifelse(band == "40-80", pair, NA_character_)]
  p <- ggplot(K, aes(band, enrich_fov, group = pair, colour = status)) +
    geom_hline(yintercept = 1, colour = MUTED, linewidth = .4) +
    geom_linerange(aes(ymin = lo, ymax = hi), alpha = .30, linewidth = .8) +
    geom_line(linewidth = .8) + geom_point(size = 2.2) +
    ggrepel::geom_text_repel(aes(label = lab), size = 2.6, na.rm = TRUE, hjust = 0,
                             direction = "y", nudge_x = .38, segment.size = .2,
                             min.segment.length = 0, seed = 1, xlim = c(4.05, NA)) +
    scale_colour_manual(values = c(likely = HODG, `near-likely` = UP), name = NULL) +
    scale_y_continuous(trans = "log2", breaks = c(0.8, 1, 1.25, 1.6, 2, 2.5)) +
    scale_x_discrete(expand = expansion(mult = c(.06, .34))) +
    labs(title = "The distance profile separates a neighbourhood association from leakage",
         subtitle = paste0("Co-expression-in-proximity enrichment across distance bands, for the eight pairs step 22 implicated.\n",
           "Shaded bars are 95% intervals from the FOV cluster bootstrap."),
         x = "centre-to-centre distance band (um)",
         y = "enrichment over expectation (log2 scale)",
         caption = paste0("Leakage is a segmentation-boundary effect and must decay within a cell diameter or two; a real neighbourhood association persists.\n",
           "PDCD1LG2-PDCD1, the largest value in the whole table, runs 2.28 -> 1.54 -> 1.14 -> 1.01: the leakage shape, on 10 candidate pairs all at\n",
           "the 2-count floor. CXCL13-CXCR5 holds across all four bands. Never quote a single band.\n",
           "Source: output/22_lr/matched/lr_proximity_tier1.csv (enrich_fov)")) + th
  sv("7d_lr_distance_profile", p, 8.6, 5.4)
}


## ===== S3 -- the measurement floor that bounds everything ==============================
## Three views of the same limit, because it constrains different claims in different ways.
AU <- rd("output/16_checkpoint/depth_auc.csv"); PB2 <- rd(file.path(PD17, "power_bounds.csv"))
if (!is.null(AU) && !is.null(PB2)) {
  A2 <- AU[thr == 2]
  A2[, kind := fifelse(gene %in% c("IGKC","CD3E","CD68","MS4A1"), "positive control",
               fifelse(gene %in% c("PDCD1","CD274","PDCD1LG2","TNFRSF8"), "gene under suspicion",
                       "other"))]
  A2 <- A2[kind != "other"]
  A2[, gene := factor(gene, levels = A2[, .(m = mean(auc_depth)), by = gene][order(m), gene])]
  s1 <- ggplot(A2, aes(auc_depth, gene, colour = samp, shape = kind)) +
    geom_vline(xintercept = 0.5, colour = MUTED, linewidth = .4) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = SAMP, name = NULL) +
    scale_shape_manual(values = c(`positive control` = 16, `gene under suspicion` = 17),
                       name = NULL) +
    labs(title = "S3a. Depth predicts positivity at AUC ~0.70 for EVERY sparse gene",
         subtitle = "Including IGKC, the transcript the FL clone is defined by. Positivity is >= 2 counts.",
         x = "AUC of sequencing depth for predicting positivity", y = NULL,
         caption = paste0("An AUC near 0.70 against depth is a property of thresholding a sparse count on this platform, NOT a signature of noise.\n",
           "Step 10 rejected CD30 partly on its 0.716, which is indistinguishable from IGKC's. Source: 16_checkpoint/depth_auc.csv")) + th
  P2 <- PB2[order(mde_fold)]
  P2[, lab := paste0(contrast, " | ", module)]
  P2[, lab := factor(lab, levels = lab)]
  s2 <- ggplot(P2, aes(mde_fold, lab, colour = abs(observed_log2) > mde_log2)) +
    geom_segment(aes(x = 1, xend = mde_fold, yend = lab), colour = "#d7d5cd", linewidth = 1.2) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(`FALSE` = MUTED, `TRUE` = UP),
                        labels = c("observed effect below the floor", "observed effect above it"),
                        name = NULL) +
    labs(title = "S3b. What each module contrast could have resolved",
         x = "minimum detectable effect (fold)", y = NULL,
         caption = paste0("Read beside every null: a contrast on 20 FOVs resolves 1.14-fold at best, one on 346 resolves 1.01-fold.\n",
           "Source: 17_pd1/matched/power_bounds.csv")) +
    th + theme(axis.text.y = element_text(size = 6))
  sv("S3_measurement_floor", s1 / s2 + patchwork::plot_layout(heights = c(.8, 1)), 9.4, 9.0)
}

## ===== S10 -- the over-fixation defect is spatially confounded =========================
DB <- rd(file.path(IF21, "defect_burden_by_region.csv"))
DE <- rd(file.path(IF21, "defect_edge_geometry.csv"))
DV <- rd(file.path(IF21, "defect_qc_verdict.csv"))
if (!is.null(DB) && !is.null(DE) && !is.null(DV)) {
  DB[, bin_region := factor(bin_region, levels = DB[order(med_bin_defect_frac), bin_region])]
  t1 <- ggplot(DB, aes(med_bin_defect_frac, bin_region)) +
    geom_linerange(aes(xmin = q1_, xmax = q3_), linewidth = 2.2, colour = "#d7d5cd") +
    geom_point(aes(colour = bin_region == "cHL_hotspot"), size = 3) +
    scale_colour_manual(values = c(`FALSE` = MUTED, `TRUE` = HODG), guide = "none") +
    labs(title = "S10a. The defect does not fall evenly on the regions being contrasted",
         subtitle = "Per-bin fraction of cells removed by the per-cell over-fixation drop. Bar = IQR.",
         x = "defect fraction of the bin's cells", y = NULL) + th
  DE2 <- melt(DE[, .(bin_region, `defect fraction` = med_defect,
                     `distance to tissue edge (um)` = med_edge_um / 400,
                     `distance to nearest FL focus (um)` = med_d_fl_um / 1600)],
              id.vars = "bin_region", variable.name = "measure", value.name = "v")
  t2 <- ggplot(DE2, aes(v, bin_region, colour = measure)) +
    geom_point(size = 3) +
    scale_colour_manual(values = c(MUTED, HODG, FLO), name = NULL) +
    labs(title = "S10b. And it is a rim, which is why it tracks distance to the FL",
         subtitle = paste0("Rescaled to a common axis: edge distance / 400 um, FL distance / 1600 um.\n",
           sprintf("rho(defect, edge distance) = %.3f; rho(defect, distance to nearest FL) = %+.3f",
                   DV$rho_defect_edge[1], DV$rho_defect_dist_to_fl[1])),
         x = "rescaled value", y = NULL,
         caption = paste0("Because the defect is a peripheral rim it rises with distance from the FL foci, so it is confounded with the outcome of any\n",
           "distance-based contrast. Consequence: carry per-bin defect fraction as a COVARIATE, not just the per-cell drop. And IB1's tissue\n",
           sprintf("'pieces' are %d before the drop and %d after -- coverage, not anatomy.\n", DV$pieces_pre[1], DV$pieces_post[1]),
           "Source: output/21_interface/matched/defect_{burden_by_region,edge_geometry,qc_verdict}.csv")) + th
  sv("S10_defect_confounding", t1 / t2, 8.8, 7.2)
}


## ===== S2 -- the QC asymmetry, which every cross-sample number inherits ================
AT <- rd("output/02_spatial_QC/qc_spatial_attrition.csv")
if (!is.null(AT)) {
  M <- melt(AT[!rule %in% c("total","passing_all","passing_spatial_only")],
            id.vars = "rule", variable.name = "samp", value.name = "cells")
  M[, samp := as.character(samp)]
  M[, rule := factor(rule, levels = M[, .(s = sum(cells)), by = rule][order(s), rule])]
  p <- ggplot(M, aes(cells, rule, fill = samp)) +
    geom_col(position = position_dodge(width = .7), width = .62) +
    scale_fill_manual(values = SAMP, name = NULL) +
    scale_x_continuous(labels = scales::comma) +
    labs(title = "S2. The spatial QC is asymmetric between the two samples",
         subtitle = paste0("Cells removed by each rule. IB1 loses 2,013,225 to a staining rule that RB1 was never evaluated for\n",
           "(qc_spatial_staining == 'not_evaluated', fail_staining = 0)."),
         x = "cells removed", y = NULL,
         caption = paste0("Any cross-sample claim on the qc_passing_all pools (IB1 400,111 vs RB1 1,399,189) inherits this asymmetry. The symmetric tier used\n",
           "throughout this project is qc_spatial_fovpx == 'ok' & qcCellsPassed, which gives IB1 1,877,697 / RB1 1,399,189 -- IB1 the LARGER pool.\n",
           "Source: output/02_spatial_QC/qc_spatial_attrition.csv")) + th
  sv("S2_qc_attrition", p, 8.4, 4.4)
}

## ===== S4 -- why the cross-sample FL comparison is not identifiable ====================
BB <- rd("output/16_regions/matched/bins_all_scales.csv.gz")
if (!is.null(BB)) {
  Z <- BB[bin_um == 100 & bin_low_info == FALSE & !is.na(fb_a)]
  Q <- Z[, {qq <- quantile(fb_a, 0:4/4)
            .(qt = paste0("Q", 1:4),
              qmean = vapply(1:4, function(i) mean(fb_a[fb_a >= qq[i] & fb_a <= qq[i+1]]), 0))},
         by = samp]
  ib_top <- Q[samp == "IB1" & qt == "Q4", qmean]; rb_bot <- Q[samp == "RB1" & qt == "Q1", qmean]
  p <- ggplot(Z, aes(fb_a, fill = samp, colour = samp)) +
    geom_density(alpha = .28, linewidth = .5) +
    geom_vline(data = Q, aes(xintercept = qmean, colour = samp), linetype = "22", linewidth = .35) +
    annotate("rect", xmin = ib_top, xmax = rb_bot, ymin = -Inf, ymax = Inf,
             fill = "#b8860b", alpha = .10) +
    annotate("text", x = (ib_top + rb_bot)/2, y = Inf, vjust = 1.6, size = 2.8, colour = "#8a6608",
             label = sprintf("no overlap:\nIB1 Q4 mean %.3f < RB1 Q1 mean %.3f", ib_top, rb_bot)) +
    scale_fill_manual(values = SAMP, name = NULL) + scale_colour_manual(values = SAMP, name = NULL) +
    labs(title = "S4. The two nodes' tissue does not overlap in B-cell composition",
         subtitle = paste0("Per-100 um-bin B share among lineage-assigned cells. Dashed lines are each sample's quartile means.\n",
           "IB1's MOST B-dense quartile is less B-rich than RB1's LEAST B-dense quartile."),
         x = "B share among lineage-assigned cells in the bin", y = "density",
         caption = paste0("This is why the direct cross-sample FL contrast is not identifiable: matching 100 um bins on B share, depth, cell count and\n",
           "typability matched 2 of 89. The difference-in-differences meant to rescue it was withdrawn, because the two reference tissues differ too.\n",
           "Step 16 quotes 0.286 vs 0.484 on its matched-bin subset; the quartile means here are over all non-low-information bins.\n",
           "Source: output/16_regions/matched/bins_all_scales.csv.gz")) + th
  sv("S4_composition_disjoint", p, 8.4, 4.8)
}

if (length(MISSING)) { cat("\nMISSING INPUTS:\n"); cat(paste0("  ", unique(MISSING), collapse = "\n"), "\n") } else
  cat("\nall inputs present\n")
writeLines(format(Sys.time()), file.path(OUT, "figure_panels.done"))
cat("\nOK\n")
