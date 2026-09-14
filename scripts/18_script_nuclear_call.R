## ============================================================================
## 18_script_nuclear_call.R -- does restricting argmax to NUCLEAR transcripts fix it?
##
## THE CAUSAL VERSION OF STEP 18. The density maps show that marker transcripts sit
## displaced toward the marker-rich neighbour, more so where argmax and the cluster
## disagree. That is a correlation. This asks the question that has an actionable answer:
## if the argmax call is recomputed from only the transcripts LEAST likely to have bled in,
## does it agree with the cluster labels more?
##
## THE CONTROL IS THE WHOLE EXPERIMENT. Nuclear transcripts are 79% of the total, so a
## nuclear-only call sees less depth and could easily agree LESS for reasons that have
## nothing to do with leakage. Comparing `nuclear` against `all` would therefore measure
## depth, not compartment. Every arm is instead compared against a RANDOM subsample of the
## same cell matched to its exact per-cell size:
##
##   nuc         Nuclear only                    vs  rand_nuc      random, n = n_nuclear
##   nonnuc      Cytoplasm + Membrane only       vs  rand_nonnuc   random, n = n_non-nuclear
##
## and the two arms predict OPPOSITE signs, which is what makes this a test rather than a
## hope. If bleed is driving miscalls, `nuc` beats its random control and `nonnuc` loses to
## its own. If both arms match their controls, compartment carries no information about
## miscalling and the leakage signal, though real, does not reach the call.
##
## THE RANDOM SUBSAMPLE IS EXACT, NOT APPROXIMATE. Drawing k of a cell's N transcripts
## without replacement is a multivariate hypergeometric, taken here as a sequence of
## rhyper() draws over the marker genes with the unmarked transcripts as the remainder.
## rhyper is vectorised over its parameters, so all ~100k cells are drawn at once per gene
## -- 35 vectorised calls rather than 3.5M scalar ones, and no need to rank 102.7M rows.
##
## THE BACKGROUND FLOOR IS RECOMPUTED PER ARM. assign_lineage()'s raw floor is the per-cell
## negprobe rate, and it is what stops a dim lineage being minted out of spillover. An arm
## holding 79% of the transcripts must be judged against 79% of the background, so the
## negprobes are compartment-split alongside the markers rather than carried over whole.
##
## Region is the unit: every interval here is a bootstrap over FOVs.
##
## Usage: sbatch --mem=64G --time=2:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/18_script_nuclear_call.R IB1
## ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(ggplot2); setDTthreads(2)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/lineage.R")

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"
OUT    <- "output/18_leakage"
FIG    <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
TOT_F  <- file.path(OUT, sprintf("%s_nuc_totals.csv.gz",  SAMPLE))
MRK_F  <- file.path(OUT, sprintf("%s_nuc_markers.csv.gz", SAMPLE))
FOC_F  <- file.path(OUT, sprintf("%s_focal_cells.csv.gz", SAMPLE))
NEG_F  <- "output/16_regions/negprobe_counts_all.Rds"
for (f in c(TOT_F, MRK_F, FOC_F, NEG_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}
SEED <- 1L; BOOT <- 400L
FAMS  <- names(LINEAGE_SETS_COARSE)
GENES <- unique(unlist(LINEAGE_SETS_COARSE, use.names = FALSE))

cat("============================================================\n")
cat("step 18 -- nuclear-restricted argmax recall: ", SAMPLE, "\n", sep = "")
cat("============================================================\n")

foc <- fread(FOC_F)[scoreable == TRUE]
tot <- fread(TOT_F)
mrk <- fread(MRK_F)
cat(sprintf("focal scoreable cells %d | totals %d | marker rows %d\n",
            nrow(foc), nrow(tot), nrow(mrk)))

## The awk aggregator carries its own hardcoded copy of the gene list. If the two ever
## drift, everything below is silently computed on the wrong markers.
extra <- setdiff(unique(mrk$target), c(GENES, "NEGPROBE"))
if (length(extra)) stop("aggregator emitted genes not in LINEAGE_SETS_COARSE: ",
                        paste(extra, collapse = ", "), call. = FALSE)
absent <- setdiff(GENES, unique(mrk$target))
if (length(absent)) cat("NOTE: no transcripts anywhere for: ", paste(absent, collapse = ", "), "\n", sep = "")

## n_neg from the project's own negprobe matrix, so the floor is the one step 16 used.
N_NEG <- nrow(readRDS(NEG_F)); invisible(gc())
cat(sprintf("negative probes: %d\n", N_NEG))

cells <- intersect(foc$cell_id, tot$cell_id)
## Both tables must be put in `cells` order, not merely keyed. setkey() SORTS, so keying
## foc and then reading foc$A positionally against a matrix built in `cells` order compares
## two different cells per row -- which is what made the regression check read 19.3%.
foc   <- foc[cell_id %in% cells]; setkey(foc, cell_id); foc <- foc[J(cells)]
tot   <- tot[cell_id %in% cells]; setkey(tot, cell_id); tot <- tot[J(cells)]
stopifnot(identical(foc$cell_id, cells), identical(tot$cell_id, cells))
tot[, n_all := n_nuc + n_cyt + n_mem]
cat(sprintf("cells scored: %d | median transcripts/cell %d (nuclear %d)\n",
            length(cells), as.integer(median(tot$n_all)), as.integer(median(tot$n_nuc))))

## ---- count matrices per compartment ----------------------------------------
CATS <- c(GENES, "NEGPROBE")
widen <- function(d) {
  w <- dcast(d[, .(n = sum(n)), by = .(cell_id, target)], cell_id ~ target,
             value.var = "n", fill = 0L)
  for (m in setdiff(CATS, names(w))) w[[m]] <- 0L
  setkey(w, cell_id)
  w <- w[J(cells)]
  for (m in CATS) set(w, which(is.na(w[[m]])), m, 0L)
  as.matrix(w[, ..CATS])          # cells x categories
}
C_all <- widen(mrk)
C_nuc <- widen(mrk[comp == "Nuclear"])
C_non <- widen(mrk[comp != "Nuclear"])
stopifnot(identical(dim(C_all), dim(C_nuc)), nrow(C_all) == length(cells))
cat(sprintf("\nmarker transcripts: all %d | nuclear %d (%.1f%%) | non-nuclear %d\n",
            sum(C_all), sum(C_nuc), 100*sum(C_nuc)/sum(C_all), sum(C_non)))

## ---- exact depth-matched random subsamples ---------------------------------
## Sequential multivariate hypergeometric: draw k of the cell's N transcripts without
## replacement, category by category, with the unmarked transcripts as the implicit
## remainder. Vectorised across cells.
set.seed(SEED)
mv_hyper <- function(cnt, Ntot, k) {
  out  <- matrix(0L, nrow(cnt), ncol(cnt), dimnames = dimnames(cnt))
  remN <- as.numeric(Ntot); remk <- as.numeric(k)
  for (j in seq_len(ncol(cnt))) {
    m <- as.numeric(cnt[, j])
    x <- rhyper(length(m), m = m, n = pmax(remN - m, 0), k = pmin(remk, remN))
    out[, j] <- as.integer(x); remN <- remN - m; remk <- remk - x
  }
  out
}
C_rnuc <- mv_hyper(C_all, tot$n_all, tot$n_nuc)
C_rnon <- mv_hyper(C_all, tot$n_all, tot$n_cyt + tot$n_mem)

## A SECOND CONTROL, MATCHED ON MARKER COUNT RATHER THAN DEPTH -- because the first one is
## not a fair fight. Marker genes are only 64.5% nuclear while transcripts overall are
## 78.3%: the panel's lineage markers are cytoplasm- and membrane-enriched. So a nuclear
## arm holding 78% of a cell's transcripts holds barely 64% of its marker signal, and it
## would lose to a depth-matched random draw even if leakage did not exist at all. That
## comparison answers "is nuclear-only calling better in practice", which is worth knowing,
## but it cannot isolate leakage.
##
## The marker-matched control fixes the amount of marker evidence and varies only WHICH
## transcripts supply it: k is the cell's NUCLEAR marker count, drawn from the cell's marker
## transcripts of every compartment. Negprobes are carried over from the matching arm
## unchanged so the background floor is identical and only composition differs.
MK   <- GENES
Mtot <- rowSums(C_all[, MK, drop = FALSE])
C_rnucM <- cbind(mv_hyper(C_all[, MK, drop = FALSE], Mtot, rowSums(C_nuc[, MK, drop = FALSE])),
                 NEGPROBE = C_nuc[, "NEGPROBE"])
C_rnonM <- cbind(mv_hyper(C_all[, MK, drop = FALSE], Mtot, rowSums(C_non[, MK, drop = FALSE])),
                 NEGPROBE = C_non[, "NEGPROBE"])
cat(sprintf("\nmarker transcripts per arm: nuc %d | rand_nuc %d (depth-matched) | rand_nuc_mk %d (marker-matched)\n",
            sum(C_nuc[, MK]), sum(C_rnuc[, MK]), sum(C_rnucM[, MK])))
cat(sprintf("                            non %d | rand_non %d (depth-matched) | rand_non_mk %d (marker-matched)\n",
            sum(C_non[, MK]), sum(C_rnon[, MK]), sum(C_rnonM[, MK])))
cat(sprintf("markers are %.1f%% nuclear against %.1f%% for transcripts overall -- the panel's lineage markers are cytoplasm-enriched\n",
            100*sum(C_nuc[, MK])/sum(C_all[, MK]), 100*sum(tot$n_nuc)/sum(tot$n_all)))

## ---- recall the lineage on each arm ----------------------------------------
call_arm <- function(C) {
  M <- t(C[, GENES, drop = FALSE])            # genes x cells
  r <- do.call(assign_lineage,
               c(list(cts = M, sets = LINEAGE_SETS_COARSE,
                      neg_total = C[, "NEGPROBE"], n_neg = N_NEG), LINEAGE_CALL_ARGS))
  r$label
}
ARMS <- list(all = C_all, nuc = C_nuc, rand_nuc = C_rnuc, rand_nuc_mk = C_rnucM,
             nonnuc = C_non, rand_nonnuc = C_rnon, rand_nonnuc_mk = C_rnonM)
lab <- lapply(ARMS, call_arm)

## REGRESSION CHECK. The `all` arm is rebuilt from transcripts, while the stored call came
## from the vendor expression matrix over all 4.4M cells; the fold-over-mean scaling basis
## therefore differs and exact agreement is not expected. A LOW number here would mean the
## transcript reconstruction is not the same data, and nothing below would be interpretable.
cat(sprintf("\nregression check -- `all` arm vs the stored lineage16 call: %.1f%% of cells\n",
            100*mean(lab$all == foc$A)))

## ---- score every arm against the cluster label -----------------------------
D <- rbindlist(lapply(names(lab), function(a)
  data.table(arm = a, cell_id = cells, lab = lab[[a]])))
D <- foc[, .(cell_id, fov_id, A, arb)][D, on = .(cell_id)]
D[, `:=`(called = lab %in% FAMS)]
D[, ok := called & lab == arb]
sc(sprintf("%s_nuclear_call_percell.csv.gz", SAMPLE), D)

summ <- function(d) d[, .(cells = .N, call_rate = mean(called),
                          agree_called = sum(ok)/pmax(sum(called),1),
                          agree_all = mean(ok)), by = arm]
S <- summ(D)
cat("\n=== per arm, all scoreable focal cells ===\n")
print(S[, .(arm, cells, call_rate = round(100*call_rate,2),
            agree_called = round(100*agree_called,2), agree_all = round(100*agree_all,2))])

## ---- the two contrasts, bootstrapped over FOVs -----------------------------
PAIRS <- list(c("nuc","rand_nuc"), c("nuc","rand_nuc_mk"),
              c("nonnuc","rand_nonnuc"), c("nonnuc","rand_nonnuc_mk"))
fovs  <- unique(D$fov_id)
stat_of <- function(d) {
  s <- summ(d); setkey(s, arm)
  vapply(PAIRS, function(p) c(s[p[1], agree_called] - s[p[2], agree_called],
                              s[p[1], agree_all]    - s[p[2], agree_all]), numeric(2))
}
obs <- stat_of(D)
set.seed(SEED + 1L)
bs <- vapply(seq_len(BOOT), function(i) {
  f <- data.table(fov_id = sample(fovs, length(fovs), replace = TRUE))
  stat_of(f[D, on = .(fov_id), nomatch = NULL, allow.cartesian = TRUE])
}, matrix(0, 2, length(PAIRS)))
res <- rbindlist(lapply(seq_along(PAIRS), function(j) rbindlist(lapply(1:2, function(i)
  data.table(contrast = paste(PAIRS[[j]], collapse = " - "),
             metric = c("agree_called","agree_all")[i],
             diff_pp = 100*obs[i,j],
             lo = 100*quantile(bs[i,j,], .025, na.rm = TRUE),
             hi = 100*quantile(bs[i,j,], .975, na.rm = TRUE),
             p_boot = 2*min(mean(bs[i,j,] <= 0, na.rm = TRUE),
                            mean(bs[i,j,] >= 0, na.rm = TRUE)))))))
cat("\n=== compartment arm MINUS its depth-matched random control (percentage points) ===\n")
print(res[, .(contrast, metric, diff_pp = round(diff_pp,3),
              ci = sprintf("[%.3f, %.3f]", lo, hi), p_boot = round(p_boot,4))])
sc(sprintf("%s_nuclear_call_contrasts.csv.gz", SAMPLE), res)

## ---- per argmax call, so it lines up with the density-map facets ------------
byA <- D[, .(cells = .N, agree_called = 100*sum(ok)/pmax(sum(called),1)), by = .(A, arm)]
byA <- dcast(byA, A + cells ~ arm, value.var = "agree_called")
byA[, `:=`(nuc_vs_depth = nuc - rand_nuc, nuc_vs_marker = nuc - rand_nuc_mk,
           non_vs_marker = nonnuc - rand_nonnuc_mk)]
cat("\n=== agreement with the cluster label (%), by argmax call ===\n")
print(byA[, lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)][order(-cells)])
sc(sprintf("%s_nuclear_call_byA.csv.gz", SAMPLE), byA)

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
pd <- melt(byA, id.vars = c("A","cells"),
           measure.vars = c("all","nuc","rand_nuc","rand_nuc_mk","nonnuc","rand_nonnuc","rand_nonnuc_mk"),
           variable.name = "arm", value.name = "agree")
pd[, arm := factor(arm, levels = c("all","nuc","rand_nuc","rand_nuc_mk","nonnuc","rand_nonnuc","rand_nonnuc_mk"))]
pd[, A := factor(A, levels = byA[order(-cells), A])]
p <- ggplot(pd, aes(arm, agree, fill = arm)) +
  geom_col(width = .72) +
  geom_text(aes(label = sprintf("%.1f", agree)), vjust = -0.4, size = 2.5, colour = INK2) +
  facet_wrap(~ A, nrow = 1) +
  scale_fill_manual(values = c(all = "#9a9892", nuc = "#184f95", rand_nuc = "#9db8d8",
                               rand_nuc_mk = "#c8d8ea", nonnuc = "#d03b3b",
                               rand_nonnuc = "#e8a9a9", rand_nonnuc_mk = "#f3d4d4"), guide = "none") +
  labs(title = sprintf("%s -- does calling argmax from nuclear transcripts only agree with the clusters better?", SAMPLE),
       subtitle = paste0("Each compartment arm sits beside TWO random subsamples of the same cells: one matched on total transcript count (rand_*),\n",
                         "one matched on marker-transcript count (rand_*_mk). The panel's markers are cytoplasm-enriched, so only the second is a fair fight."),
       x = NULL, y = "agreement with pearson_clusters_low (%, among cells called)") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = GRID, linewidth = .25),
        axis.text.x = element_text(colour = MUTED, size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = MUTED), axis.title = element_text(colour = INK2),
        strip.text = element_text(colour = INK2, face = "bold"),
        plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = INK2, size = 9))
sf(sprintf("%s_nuclear_call.png", SAMPLE), p, 15, 4.6)
cat("\ndone\n")
