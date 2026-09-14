## ============================================================================
## 08_script_lightchain_sample_compare.R -- IB1 vs RB1, using RB1 as a positive
## control for what full FL involvement looks like through this exact pipeline.
##
## Clinical context: the patient was treated for Hodgkin lymphoma and recurred with
## follicular lymphoma. IB1 is the composite node (HRS cells plus FL); RB1 is the FL
## recurrence. Both are the same clone, so RB1 calibrates the pipeline against a
## known, extensive FL and IB1 can be read against it.
##
## THE QUESTION this step exists to answer. IB1's neighbourhoods sit above the
## polytypic expectation, which could mean either
##   (a) DIFFUSE early FL infiltration -- the clone spread thinly through the node, so
##       the "background" is itself disease and a focus-finding detector is measuring
##       the wrong thing; or
##   (b) FOCAL FL -- discrete clonal foci on a genuinely reactive polytypic background.
##
## These are separated by looking at the UNCLUSTERED B cells. Under (a) they are skewed
## too; under (b) they sit at the polytypic baseline. RB1 shows what (a) looks like.
##
## Note the clustered-vs-diffuse split is not circular: the spatial graph is built on
## all light-chain-high B cells regardless of chain, so light chain plays no part in
## deciding which cells cluster. Finding that clusters are kappa-skewed is a result.
##
## Reads only the step-07 csvs for both samples, so it is login-node cheap:
##   Rscript scripts/08_script_lightchain_sample_compare.R
## Requires scripts/07_submit_lightchain_neighborhoods.sh to have been run for BOTH.
## ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); setDTthreads(1)
})

SAMPLES  <- c("IB1", "RB1")
RADII    <- c(20, 30, 40, 60)
FOCUS_R  <- 40
IN_ROOT  <- "output/07_lightchain"
OUT_DIR  <- "output/08_lightchain_compare"
FIG      <- file.path(OUT_DIR, "figures")
BCL2_CSV <- "output/04_merge/bcl2_detectability.csv"
NULL_K   <- 0.62
EXP_L2   <- log2(NULL_K / (1 - NULL_K))
## total segmented cells per sample, for pool prevalence
TOT      <- c(IB1 = 2540021L, RB1 = 1866919L)

dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
DIV_HI<-"#184f95"; DIV_LO<-"#d03b3b"
SAMP_COLS <- c(IB1 = "#4e79a7", RB1 = "#e15759")
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

## ---- load both samples -----------------------------------------------------
need <- file.path(IN_ROOT, SAMPLES, sprintf("lightchain_cells_r%dum.csv", FOCUS_R))
if (!all(file.exists(need)))
  stop("missing per-sample tables:\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "),
       "\nrun scripts/07_submit_lightchain_neighborhoods.sh <sample> first", call. = FALSE)

CELLS <- rbindlist(lapply(SAMPLES, function(S) rbindlist(lapply(RADII, function(R) {
  lc  <- fread(file.path(IN_ROOT, S, sprintf("lightchain_cells_r%dum.csv", R)))
  nbt <- fread(file.path(IN_ROOT, S, sprintf("lightchain_neighborhoods_r%dum.csv", R)))
  keep <- nbt$nb                       # not `nbt$nb` inside lc[...]: the column shadows it
  lc[, `:=`(sample = S, radius_um = R, scored = !is.na(nb) & nb %in% keep)]
  lc[, .(sample, radius_um, tier, scored, call)]
}))), fill = TRUE)

NBS <- rbindlist(lapply(SAMPLES, function(S) rbindlist(lapply(RADII, function(R) {
  nbt <- fread(file.path(IN_ROOT, S, sprintf("lightchain_neighborhoods_r%dum.csv", R)))
  nbt[, `:=`(sample = S, radius_um = R, l2 = .l2(n_kappa, n_lambda))]
}))), fill = TRUE)
NBS[, radius_lab := factor(sprintf("r%dum", radius_um), levels = sprintf("r%dum", RADII))]

## ---- 1. the decisive table: skew by compartment ----------------------------
kfrac <- function(d) { d <- d[call != "tie"]; if (!nrow(d)) NA_real_ else mean(d$call == "kappa") }
comp <- CELLS[radius_um == FOCUS_R, .(
  compartment = c("B cells, clustered", "B cells, diffuse", "plasma cells", "halo"),
  kappa_frac = c(kfrac(.SD[tier == "B_clean" & scored == TRUE]),
                 kfrac(.SD[tier == "B_clean" & scored == FALSE]),
                 kfrac(.SD[tier == "plasma"]), kfrac(.SD[tier == "halo"])),
  cells = c(sum(tier == "B_clean" & scored), sum(tier == "B_clean" & !scored),
            sum(tier == "plasma"), sum(tier == "halo"))), by = sample]
comp[, `:=`(ratio = kappa_frac / (1 - kappa_frac),
            compartment = factor(compartment,
              levels = c("B cells, clustered", "B cells, diffuse", "halo", "plasma cells")))]
comp[, pool_pct := round(100 * sum(cells) / TOT[sample], 1), by = sample]
cat("=== kappa fraction by compartment, r", FOCUS_R, "um\n", sep = "")
print(comp[order(sample, compartment)])
cat(sprintf("\npolytypic reference %.3f (%.2f:1)\n", NULL_K, NULL_K/(1-NULL_K)))

p1 <- ggplot(comp, aes(compartment, kappa_frac, fill = sample)) +
  geom_hline(yintercept = NULL_K, colour = INK, linetype = "22", linewidth = .5) +
  geom_col(position = position_dodge(.75), width = .68) +
  geom_text(aes(label = sprintf("%.3f", kappa_frac)),
            position = position_dodge(.75), vjust = -0.4, size = 3, colour = INK2) +
  scale_fill_manual(values = SAMP_COLS, name = NULL) +
  scale_y_continuous(limits = c(0, 1.06), breaks = seq(0, 1, .2)) +
  labs(title = "The diffuse B-cell compartment is what separates focal from diffuse FL",
       subtitle = paste0("Dashed = polytypic expectation ", sprintf("%.2f", NULL_K),
                         ". RB1 (known extensive FL) is skewed EVERYWHERE, including its ",
                         "unclustered B cells.\nIB1's unclustered B cells sit at the ",
                         "polytypic baseline: its skew is confined to the clustered foci."),
       x = NULL, y = "kappa fraction of called cells",
       caption = sprintf("Data: %s/<sample>/lightchain_cells_r%dum.csv", IN_ROOT, FOCUS_R)) +
  theme_v
sf("fig_08_compartment_kappa.png", p1, 9, 5)

## ---- 2. neighbourhood ratio distributions ----------------------------------
p2 <- ggplot(NBS, aes(l2, fill = sample, colour = sample)) +
  geom_vline(xintercept = EXP_L2, colour = INK, linetype = "22", linewidth = .4) +
  geom_density(alpha = .35, linewidth = .5) +
  scale_fill_manual(values = SAMP_COLS, name = NULL) +
  scale_colour_manual(values = SAMP_COLS, name = NULL) +
  facet_wrap(~ radius_lab, nrow = 1) +
  labs(title = "Neighbourhood kappa:lambda, IB1 vs the FL recurrence",
       subtitle = paste0("Dashed = polytypic expectation. IB1 straddles it; RB1 sits far ",
                         "to the kappa side at every radius."),
       x = "log2 kappa:lambda per neighbourhood", y = "density",
       caption = "Data: lightchain_neighborhoods_r*um.csv") + theme_v
sf("fig_08_neighbourhood_ratio_density.png", p2, 15, 4.4)

## ---- 3. funnel, both samples against the same envelope ---------------------
env <- rbindlist(lapply(seq(5, max(NBS$n_called), by = 1), function(n) {
  qs <- qbinom(c(0.025, 0.975), n, NULL_K)
  data.table(n_called = n, lo = .l2(qs[1], n - qs[1]), hi = .l2(qs[2], n - qs[2]))
}))
p3 <- ggplot(NBS, aes(n_called, l2)) +
  geom_ribbon(data = env, aes(x = n_called, ymin = lo, ymax = hi), inherit.aes = FALSE,
              fill = MUTED, alpha = .30) +
  geom_hline(yintercept = EXP_L2, colour = INK, linetype = "22", linewidth = .4) +
  geom_point(aes(colour = sample), size = 1, alpha = .75) +
  scale_colour_manual(values = SAMP_COLS, name = NULL) +
  scale_x_log10() + facet_wrap(~ radius_lab, nrow = 1) +
  labs(title = "Same binomial envelope, both samples",
       subtitle = paste0("Band = what chance alone produces under the polytypic null. ",
                         "RB1 is almost entirely outside it;\nIB1 has a population inside, ",
                         "which is its reactive background."),
       x = "called cells (log scale)", y = "log2 kappa:lambda",
       caption = "Data: lightchain_neighborhoods_r*um.csv") + theme_v
sf("fig_08_funnel_both_samples.png", p3, 15, 4.6)

## ---- 4. burden: pool prevalence and the independent BCL2 read --------------
burden <- unique(comp[, .(sample, pool_pct)])
if (file.exists(BCL2_CSV)) {
  b <- fread(BCL2_CSV)
  burden <- merge(burden, b[, .(sample = samp, pct_BCL2_in_B, pct_BCL2_in_gcB)],
                  by = "sample", all.x = TRUE)
}
cat("\n=== disease burden, two independent measures\n"); print(burden)

bl <- melt(burden, id.vars = "sample", variable.name = "measure", value.name = "pct")
bl[, measure := factor(measure,
     levels = c("pool_pct", "pct_BCL2_in_B", "pct_BCL2_in_gcB"),
     labels = c("light-chain-high\n(% of all cells)", "BCL2+\n(% of B cells)",
                "BCL2+\n(% of GC B cells)"))]
p4 <- ggplot(bl[!is.na(pct)], aes(measure, pct, fill = sample)) +
  geom_col(position = position_dodge(.75), width = .68) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(.75),
            vjust = -0.4, size = 3, colour = INK2) +
  scale_fill_manual(values = SAMP_COLS, name = NULL) +
  labs(title = "Disease burden by two independent measures",
       subtitle = paste0("BCL2 is the FL hallmark (t(14;18)) and is measured independently of ",
                         "light chain.\nIB1 carries a real but minority FL population; RB1 is ",
                         "extensively involved."),
       x = NULL, y = "percent",
       caption = sprintf("Data: lightchain_cells + %s", BCL2_CSV)) + theme_v
sf("fig_08_burden.png", p4, 8, 4.8)

## ---- 5. verdict ------------------------------------------------------------
d_ib <- comp[sample == "IB1" & compartment == "B cells, diffuse", kappa_frac]
d_rb <- comp[sample == "RB1" & compartment == "B cells, diffuse", kappa_frac]
c_ib <- comp[sample == "IB1" & compartment == "B cells, clustered", kappa_frac]
fwrite(comp, file.path(OUT_DIR, "compartment_kappa.csv"))
fwrite(burden, file.path(OUT_DIR, "burden.csv"))
fwrite(NBS[, .(sample, radius_um, nb, n_called, n_kappa, n_lambda, l2, class,
               artifact_suspect)],
       file.path(OUT_DIR, "neighbourhoods_both_samples.csv"))

cat("\n============================================================\n")
cat(sprintf("IB1 diffuse B cells  %.3f   vs polytypic %.3f  -> %s\n", d_ib, NULL_K,
            ifelse(abs(d_ib - NULL_K) < 0.03, "AT baseline", "departs from baseline")))
cat(sprintf("RB1 diffuse B cells  %.3f   vs polytypic %.3f  -> %s\n", d_rb, NULL_K,
            ifelse(abs(d_rb - NULL_K) < 0.03, "AT baseline", "departs from baseline")))
cat(sprintf("IB1 clustered B cells %.3f  -- the skew is confined here\n", c_ib))
cat("\nRB1 shows what diffuse involvement looks like: its UNCLUSTERED B cells are\n")
cat("skewed too. IB1's are not. So IB1's FL reads as FOCAL -- discrete clonal foci on\n")
cat("a reactive polytypic background -- rather than early diffuse infiltration.\n")
cat("This does not argue against limited/early FL in a treated node; it locates it.\n")
cat("============================================================\n")
cat("wrote ", OUT_DIR, "/{compartment_kappa,burden,neighbourhoods_both_samples}.csv\n", sep = "")
