## ============================================================================
## 07_script_lightchain_size_diagnostics.R -- how does a neighbourhood's
## kappa:lambda ratio depend on its SIZE, in cells and in area?
##
## Written to support choosing MIN_CELLS. The motivating observation: at r40um the
## neighbourhoods that tilt lambda are markedly smaller -- median 24 called cells and
## 0.013 mm2, against 47 cells and 0.034 mm2 for those above the expectation.
##
## Two very different things could produce that, and they call for opposite responses:
##   1. SAMPLING NOISE. With n=24 and an expected 62% kappa, chance alone spans a wide
##      range of ratios. But noise is SYMMETRIC -- it would throw up small
##      extreme-kappa neighbourhoods as often as small extreme-lambda ones.
##   2. A STRUCTURAL SIZE-RATIO LINK. This sample is kappa-dominant, so lambda cells
##      are sparser; a community that is mostly lambda is therefore smaller by
##      construction. If so, raising MIN_CELLS preferentially deletes lambda-tilted
##      neighbourhoods and makes the sample look MORE kappa-skewed as an artefact of
##      the cutoff, not as a finding.
##
## The funnel plot is what tells these apart: it draws the binomial sampling envelope
## at the null, so a point outside the envelope is more extreme than chance explains at
## that n. The cutoff table then shows what each candidate MIN_CELLS would retain AND
## what it would remove, split by tilt direction, so the selection effect is visible
## before it is applied.
##
## Reads only the step-07 csvs -- login-node cheap, no Seurat object, no TIFF.
##
## Usage: Rscript scripts/07_script_lightchain_size_diagnostics.R [in_dir]
## ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); setDTthreads(1)
})

args   <- commandArgs(trailingOnly = TRUE)
IN_DIR <- if (length(args) >= 1) args[1] else "output/07_lightchain/IB1"
FIG    <- file.path(IN_DIR, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

RADII  <- c(20, 30, 40, 60)
FOCUS_R <- 40          # the radius under consideration
NULL_K <- 0.62
EXP_L2 <- log2(NULL_K / (1 - NULL_K))
CUTOFFS <- c(20, 25, 30, 40, 50, 60, 80)

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
DIV_HI<-"#184f95"; DIV_LO<-"#d03b3b"
theme_v <- theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = GRID, linewidth = .25),
        plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = INK2, size = 9),
        plot.caption = element_text(colour = MUTED, size = 8, hjust = 0),
        strip.text = element_text(colour = INK2, face = "bold"))
sf <- function(n, p, w, h) {
  ggsave(file.path(FIG, n), p, width = w, height = h, dpi = 200, bg = SURF)
  cat("wrote ", file.path(FIG, n), "\n", sep = "")
}
.l2 <- function(k, l) log2((k + 0.5) / (l + 0.5))

## ---- 1. assemble neighbourhoods with size and area -------------------------
## Area = convex hull of the member cells, in mm2 (shoelace on the hull). The table's
## own span_mm is a max radial distance, so it does not describe a non-circular focus.
hull_area <- function(x, y) {
  if (length(x) < 3L) return(NA_real_)
  h <- grDevices::chull(x, y); hx <- x[h]; hy <- y[h]
  nx <- c(hx[-1], hx[1]); ny <- c(hy[-1], hy[1])
  abs(sum(hx * ny - nx * hy)) / 2
}

A <- rbindlist(lapply(RADII, function(R) {
  ## NOT named `nb`: inside lc[...] the COLUMN nb would shadow it and `nb$nb` fails
  nbt <- fread(file.path(IN_DIR, sprintf("lightchain_neighborhoods_r%dum.csv", R)))
  lc  <- fread(file.path(IN_DIR, sprintf("lightchain_cells_r%dum.csv", R)))
  keep_nb <- nbt$nb
  ar <- lc[nb %in% keep_nb, .(area_mm2 = hull_area(x_slide_mm, y_slide_mm)), by = nb]
  nbt <- merge(nbt, ar, by = "nb", all.x = TRUE)
  nbt[, `:=`(radius_um = R, l2 = .l2(n_kappa, n_lambda))]
  nbt[, `:=`(tilt = ifelse(l2 < EXP_L2, "below expectation (lambda-ward)",
                           "above expectation (kappa-ward)"),
             dens = n_called / area_mm2)]
  nbt
}), fill = TRUE)
A[, radius_lab := factor(sprintf("r%dum", radius_um),
                         levels = sprintf("r%dum", RADII))]
cat("neighbourhoods per radius:\n"); print(A[, .N, by = radius_lab])

## ---- 2. binomial envelope: what could chance alone produce? ----------------
env <- rbindlist(lapply(seq(5, max(A$n_called), by = 1), function(n) {
  qs <- qbinom(c(0.025, 0.975, 0.005, 0.995), n, NULL_K)
  data.table(n_called = n,
             lo95 = .l2(qs[1], n - qs[1]), hi95 = .l2(qs[2], n - qs[2]),
             lo99 = .l2(qs[3], n - qs[3]), hi99 = .l2(qs[4], n - qs[4]))
}))
A[, out95 := l2 < env$lo95[match(n_called, env$n_called)] |
              l2 > env$hi95[match(n_called, env$n_called)]]
cat("\noutside the 95% binomial envelope, by radius and tilt:\n")
print(dcast(A[, .N, by = .(radius_lab, tilt, out95)], radius_lab + tilt ~ out95,
            value.var = "N", fill = 0))

## ---- 3. figure: funnel, ratio vs number of cells ---------------------------
p1 <- ggplot(A, aes(n_called, l2)) +
  geom_ribbon(data = env, aes(x = n_called, ymin = lo99, ymax = hi99),
              inherit.aes = FALSE, fill = GRID, alpha = .55) +
  geom_ribbon(data = env, aes(x = n_called, ymin = lo95, ymax = hi95),
              inherit.aes = FALSE, fill = MUTED, alpha = .30) +
  geom_hline(yintercept = EXP_L2, colour = INK, linetype = "22", linewidth = .4) +
  geom_vline(xintercept = CUTOFFS, colour = AXIS <- "#c3c2b7", linewidth = .2) +
  geom_point(aes(colour = out95, size = area_mm2), alpha = .85) +
  scale_colour_manual(values = c(`TRUE` = DIV_HI, `FALSE` = MUTED),
                      labels = c(`TRUE` = "beyond chance", `FALSE` = "within chance"),
                      name = NULL) +
  scale_size_area(max_size = 4, name = "area mm2") +
  scale_x_log10() + facet_wrap(~ radius_lab, nrow = 1) +
  labs(title = "Ratio vs neighbourhood size, against the binomial sampling envelope",
       subtitle = paste0("Bands = what chance alone produces at each n under the polytypic null ",
                         "(dark 95%, light 99%). Dashed = expected ", sprintf("%.2f", EXP_L2),
                         ".\nA point inside the band is not more skewed than its size explains. ",
                         "Vertical lines mark the candidate MIN_CELLS cutoffs."),
       x = "called cells in the neighbourhood (log scale)", y = "log2 kappa:lambda",
       caption = "Data: lightchain_neighborhoods_r*um.csv") + theme_v
sf("fig_07_size_funnel_cells.png", p1, 15, 5)

## ---- 4. figure: ratio vs area ----------------------------------------------
p2 <- ggplot(A, aes(area_mm2, l2)) +
  geom_hline(yintercept = EXP_L2, colour = INK, linetype = "22", linewidth = .4) +
  geom_point(aes(colour = tilt, size = n_called), alpha = .8) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
              colour = INK2, fill = GRID, linewidth = .5) +
  scale_colour_manual(values = c("above expectation (kappa-ward)" = DIV_HI,
                                 "below expectation (lambda-ward)" = DIV_LO), name = NULL) +
  scale_size_area(max_size = 4, name = "called cells") +
  scale_x_log10() + facet_wrap(~ radius_lab, nrow = 1) +
  labs(title = "Ratio vs neighbourhood area",
       subtitle = "Hull area of the member cells. If the trend is flat, area is not driving the ratio.",
       x = "convex hull area, mm2 (log scale)", y = "log2 kappa:lambda",
       caption = "Data: lightchain_neighborhoods_r*um.csv") + theme_v
sf("fig_07_size_ratio_area.png", p2, 15, 5)

## ---- 5. figure: are the lambda-ward ones systematically smaller? -----------
mlt <- melt(A[, .(radius_lab, tilt, `called cells` = n_called,
                  `area mm2` = area_mm2, `cells per mm2` = dens)],
            id.vars = c("radius_lab", "tilt"))
p3 <- ggplot(mlt, aes(tilt, value, fill = tilt)) +
  geom_violin(scale = "width", colour = NA, alpha = .8) +
  geom_boxplot(width = .15, outlier.size = .5, fill = SURF, colour = INK2, linewidth = .3) +
  scale_fill_manual(values = c("above expectation (kappa-ward)" = DIV_HI,
                              "below expectation (lambda-ward)" = DIV_LO), guide = "none") +
  scale_y_log10() +
  facet_grid(variable ~ radius_lab, scales = "free_y") +
  labs(title = "Size of the neighbourhoods that tilt lambda vs those that tilt kappa",
       subtitle = paste0("If the lambda-ward ones are consistently smaller, a higher MIN_CELLS ",
                         "removes them preferentially\nand the sample will look more kappa-skewed ",
                         "as a consequence of the cutoff."),
       x = NULL, y = NULL, caption = "Data: lightchain_neighborhoods_r*um.csv") +
  theme_v + theme(axis.text.x = element_blank())
sf("fig_07_size_tilt_vs_size.png", p3, 13, 8)

## ---- 6. what each candidate cutoff would do -------------------------------
cut_tab <- rbindlist(lapply(CUTOFFS, function(cut) A[, {
  keep <- n_called >= cut
  .(min_cells = cut, kept = sum(keep), dropped = sum(!keep),
    kept_lambda_ward = sum(keep & tilt %like% "below"),
    dropped_lambda_ward = sum(!keep & tilt %like% "below"),
    pct_kept_lambda_ward = round(100 * mean(tilt[keep] %like% "below"), 1),
    med_l2_kept = round(median(l2[keep]), 2),
    fl_like_kept = sum(keep & class == "kappa-restricted (FL-like)" & !artifact_suspect),
    beyond_chance_kept = sum(keep & out95))
}, by = radius_lab]))
cat("\n=== what each candidate MIN_CELLS would retain\n")
print(cut_tab[radius_lab == sprintf("r%dum", FOCUS_R)])
fwrite(cut_tab, file.path(IN_DIR, "size_cutoff_tradeoff.csv"))
cat("wrote ", file.path(IN_DIR, "size_cutoff_tradeoff.csv"), "\n", sep = "")

p4 <- ggplot(cut_tab, aes(min_cells)) +
  geom_line(aes(y = kept, colour = "kept"), linewidth = .7) +
  geom_point(aes(y = kept, colour = "kept")) +
  geom_line(aes(y = kept_lambda_ward, colour = "kept, lambda-ward"), linewidth = .7) +
  geom_point(aes(y = kept_lambda_ward, colour = "kept, lambda-ward")) +
  geom_line(aes(y = beyond_chance_kept, colour = "kept, beyond chance"), linewidth = .7) +
  geom_point(aes(y = beyond_chance_kept, colour = "kept, beyond chance")) +
  scale_colour_manual(values = c(kept = INK2, `kept, lambda-ward` = DIV_LO,
                                 `kept, beyond chance` = DIV_HI), name = NULL) +
  facet_wrap(~ radius_lab, nrow = 1, scales = "free_y") +
  labs(title = "Cost of each candidate MIN_CELLS",
       subtitle = paste0("`beyond chance` = outside the 95% binomial envelope, i.e. more skewed ",
                         "than its size explains --\nthe neighbourhoods worth keeping. Watch how ",
                         "fast the lambda-ward ones disappear."),
       x = "MIN_CELLS", y = "neighbourhoods",
       caption = "Data: size_cutoff_tradeoff.csv") + theme_v
sf("fig_07_size_cutoff_tradeoff.png", p4, 15, 4.6)

## ---- 7. the summary the decision actually rests on ------------------------
cat("\n=== r", FOCUS_R, "um focus: correlation of ratio with size\n", sep = "")
f <- A[radius_um == FOCUS_R]
cat(sprintf("  Spearman rho(log2 ratio, called cells) = %.3f\n",
            cor(f$l2, f$n_called, method = "spearman")))
cat(sprintf("  Spearman rho(log2 ratio, area mm2)     = %.3f\n",
            cor(f$l2, f$area_mm2, method = "spearman", use = "complete.obs")))
cat(sprintf("  lambda-ward: n=%d, median %d cells, %.4f mm2\n",
            f[tilt %like% "below", .N], f[tilt %like% "below", median(n_called)],
            f[tilt %like% "below", median(area_mm2, na.rm = TRUE)]))
cat(sprintf("  kappa-ward : n=%d, median %d cells, %.4f mm2\n",
            f[tilt %like% "above", .N], f[tilt %like% "above", median(n_called)],
            f[tilt %like% "above", median(area_mm2, na.rm = TRUE)]))
cat(sprintf("  outside the 95%% envelope: %d of %d (%d lambda-ward, %d kappa-ward)\n",
            sum(f$out95), nrow(f), sum(f$out95 & f$tilt %like% "below"),
            sum(f$out95 & f$tilt %like% "above")))
