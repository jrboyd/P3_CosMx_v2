#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 22b -- LIGAND-RECEPTOR CO-EXPRESSION IN PROXIMITY.  The G4 gap.
##
## WHY THIS GAP MATTERED. The whole project narrative rests on two receptor-ligand axes --
## the chemokine axis (CCL17/CCL22 -> CCR4) and the checkpoint axis (PD-L1/PD-L2 -> PD-1) --
## and both had only ever been measured as separate single-gene gradients. No cell-cell
## communication analysis of any kind existed in the repo.
##
## WHAT IS MEASURED, AND WHAT THE WORDS MEAN. For a curated ligand-receptor pair (L, R) and
## a distance band, the rate at which a cell within that band of an L-positive cell is
## R-positive. **This is CO-EXPRESSION IN PROXIMITY: not signalling, not engagement, not
## communication, and this step will not use those words for its results.** The distinction
## is not pedantry. Step 17 built three nulls for transcript-level PD-1/PD-L1 proximity, all
## three failed their own positive control, and the readout was withdrawn under a
## pre-committed stopping rule. Spatial engagement of a receptor-ligand axis is UNMEASURED in
## this project and nothing here changes that.
##
## THE NULL IS ANALYTIC, NOT SIMULATED, AND THAT IS WHAT MAKES THE STEP POSSIBLE.
## The first version of this script permuted positivity labels 300 times per pair per band
## per region: 2,315 pairs x 5 bands x 7 regions x 300 permutations over million-element
## neighbour lists, i.e. days of compute. It is unnecessary, because the permutation null has
## a closed form.
##
## Write the band's neighbour graph as a sparse adjacency A (source cell -> target cell) and
## per-cell positivity as a sparse logical matrix P (gene x cell). Then for every pair at once
##      OBS = P A P^T / sum(A)
## is one sparse triple product. And if positivity is permuted independently within strata,
## the expected rate depends on the strata only, so with S the stratum indicator (stratum x
## cell) and Pbar = P S^T / rowSums(S) the per-stratum positivity rates,
##      EXP = Pbar (S A S^T) Pbar^T / sum(A)
## where S A S^T is a tiny stratum-by-stratum matrix. Both are exact, and the whole 1,125-gene
## square comes out in seconds.
##
## THE STRATA ARE DEPTH DECILE x CELL TYPE, and both terms are load-bearing. Depth, because
## positivity tracks sequencing depth at AUC ~0.70 for every gene on this platform, positive
## controls included. Cell type, because cell types cluster spatially, so a pair whose ligand
## and receptor are expressed by the same lineage would show enrichment from lineage
## clustering alone. Stratifying on both means the observed arrangement of depths and
## lineages is held fixed and only which cell within a stratum is positive varies.
##
## THE LEAKAGE DISCRIMINATOR, which the design exists to survive. Step 18 measured marker
## transcripts displaced toward marker-rich neighbours (+0.127 log2 in myeloid, a floor). A
## statistic of the form "ligand here, receptor next door" is EXACTLY what leakage fabricates.
## So the statistic is computed across 0-10, 10-20, 20-40, 40-80 and 80-160 um. Leakage is a
## segmentation-boundary effect and must decay within a cell diameter or two; a genuine
## neighbourhood association persists. This is step 16's BATF3-vs-distance-to-cDC idiom, which
## found the BATF3 ratio HIGHEST farthest from any cDC and so not spillover. A pair enriched
## only in the first band is reported as leakage-shaped, never as biology.
##
## TWO TIERS, AS IN STEP 21. Tier 1 is the set of axes this project already declared, so it
## carries the power and cannot be a post-hoc pick from 2,315. Tier 2 is the full curated
## scan, corrected across everything, run so the discovery surface has been looked at.
## Inference is a FOV cluster bootstrap, restricted to Tier 1 plus screened Tier 2 hits so
## the per-FOV aggregates fit in memory.
##
## Cell types are STEP 19's, which also closes part of gap G1. The `Hodgkin chemokine niche`
## label is excluded from cell-type roles: CLAUDE.md is explicit that it is a region, not a
## lineage. Within-sample only -- no IB1-vs-RB1 subtraction anywhere.
##
## Usage: sbatch --mem=128G --time=4:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/22_script_lr_proximity.R [tier] [smoke]
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(Matrix); library(RANN); library(data.table)})
setDTthreads(1); options(width = 200)

AG    <- commandArgs(TRUE)
TIER  <- if (length(AG) >= 1) AG[1] else "matched"
SMOKE <- length(AG) >= 2 && AG[2] == "smoke"
LR   <- "output/22_lr"
OUT  <- if (SMOKE) file.path(LR, TIER, "smoke") else file.path(LR, TIER)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

## Bands stop at 80 um, and K_NN is 200, both for the same measured reason: with a k-nearest
## cap the far bands are only sampled where tissue is sparse, which biases them. In IB1's FL
## foci the 60th neighbour sits at 48 um median and the 200th at 106 um, so a 80-160 band
## under a 60-neighbour cap was populated almost entirely by cells at the edge of the focus.
## 0-10 vs 40-80 is enough to separate a segmentation-boundary effect from a neighbourhood
## one, which is all the discriminator needs. Saturation is measured per region per band and
## a truncated band is flagged, not silently reported.
BANDS     <- c(0, 10, 20, 40, 80)
POS_MIN   <- 2L          # >=2 counts, the project's positivity definition
K_NN      <- 200L
N_DEC     <- 10L
MIN_POS   <- 50L         # positive cells needed, per side, per region
MIN_PAIRS <- 500L
MIN_EXP   <- 10       # expected co-positive pairs; below this obs/exp is unstable
MAX_SAT   <- 0.20     # a band is 'truncated' if >20% of cells cannot reach its upper edge
SCREEN_E  <- 1.15        # Tier 2 screen before the (costlier) bootstrap
NBOOT     <- 2000L
Q_CUT     <- 0.05
SEED      <- 22L
NICHE     <- "Hodgkin chemokine niche"
if (SMOKE) { NBOOT <- 50L; K_NN <- 60L }

## Tier 1: declared here, before any scan.
T1 <- unique(data.table(
  ligand   = c("CCL17","CCL22","CCL19","CCL21","CXCL13","CXCL9","CXCL10","CXCL11",
               "CD274","PDCD1LG2","CD40LG","ICOSLG","IL6"),
  receptor = c("CCR4","CCR4","CCR7","CCR7","CXCR5","CXCR3","CXCR3","CXCR3",
               "PDCD1","PDCD1","CD40","ICOS","IL6R")))

cat("========================================================================\n")
if (SMOKE) cat("*** SMOKE MODE -> ", OUT, " ***\n", sep = "")
cat("step 22b -- LR co-expression in proximity, tier =", TIER, "\n")
cat("========================================================================\n")

M  <- readRDS(file.path(LR, "lr_counts.Rds"))
CE <- fread(file.path(LR, "lr_cells.csv.gz"))
PR <- unique(fread(file.path(LR, "lr_pairs_used.csv"))[, .(ligand, receptor)])
setkey(CE, cell_id); CE <- CE[colnames(M)]
stopifnot(identical(CE$cell_id, colnames(M)))

RG <- fread(file.path("output/16_regions", TIER, "cell_regions.csv.gz"),
            select = c("cell_id","bin_region"))
CE <- merge(CE, RG, by = "cell_id", all.x = TRUE, sort = FALSE)
setkey(CE, cell_id); CE <- CE[colnames(M)]
CE[, ct  := fifelse(celltype_cluster == NICHE, NA_character_, celltype_cluster)]
CE[, dep := cut(nCount_RNA, quantile(nCount_RNA, seq(0,1,length.out = N_DEC+1), na.rm = TRUE),
                include.lowest = TRUE, labels = FALSE)]
CE[, ok := !is.na(bin_region) & bin_region != "unbinned" & !is.na(ct) & !is.na(dep)]
cat(sprintf("\ncells usable: %s of %s\n", format(CE[, sum(ok)], big.mark=","),
            format(nrow(CE), big.mark=",")))
print(CE[ok == TRUE, .(cells=.N, fovs=uniqueN(fov_id)), by=.(samp,bin_region)][order(samp,-cells)])

PR[, tier := fifelse(paste(ligand,receptor) %in% paste(T1$ligand,T1$receptor),
                     "T1_declared","T2_database")]
PR <- PR[ligand %in% rownames(M) & receptor %in% rownames(M)]
cat(sprintf("\npairs: Tier 1 %d | Tier 2 %d\n",
            PR[tier=="T1_declared", .N], PR[tier=="T2_database", .N]))
if (!PR[tier=="T1_declared", .N]) stop("no Tier 1 pair is on the panel", call. = FALSE)

POS <- as(as(M >= POS_MIN, "dMatrix"), "generalMatrix")     # gene x cell, 1/0
GEN <- rownames(POS)
gi  <- setNames(seq_along(GEN), GEN)

REG <- CE[ok == TRUE, .N, by = .(samp, bin_region)][N >= 5000L][order(samp, -N)]
if (SMOKE) REG <- REG[1:min(2, .N)]
set.seed(SEED)
stage1 <- list(); fovagg <- list()

for (r in seq_len(nrow(REG))) {
  SM <- REG$samp[r]; BR <- REG$bin_region[r]
  z  <- CE[ok == TRUE & samp == SM & bin_region == BR]
  ii <- match(z$cell_id, colnames(POS))
  P  <- POS[, ii, drop = FALSE]
  n  <- nrow(z)
  xy <- as.matrix(z[, .(x_slide_mm, y_slide_mm)]) * 1000
  nn <- RANN::nn2(xy, xy, k = min(K_NN, n - 1L) + 1L)
  NI <- nn$nn.idx[, -1, drop = FALSE]; ND <- nn$nn.dists[, -1, drop = FALSE]
  ## k-NN saturation: the fraction of cells whose furthest retrieved neighbour falls short
  ## of a band's upper edge. Those cells contribute nothing to that band, and if they are
  ## many the band describes sparse tissue rather than the region.
  dmax <- ND[, ncol(ND)]
  sat  <- vapply(BANDS[-1], function(h) mean(dmax < h), 0)
  names(sat) <- sprintf("<%d", BANDS[-1])
  cat("    kNN saturation (frac of cells unable to reach the band edge): ",
      paste(sprintf("%s %.2f", names(sat), sat), collapse = "  "), "\n", sep = "")

  ## strata = depth decile x cell type; S is (stratum x cell)
  z[, stratum := paste(dep, ct, sep = "|")]
  sl <- sort(unique(z$stratum)); si <- setNames(seq_along(sl), sl)
  S  <- sparseMatrix(i = si[z$stratum], j = seq_len(n), x = 1, dims = c(length(sl), n))
  Pbar <- (P %*% t(S)) %*% Diagonal(x = 1/pmax(Matrix::rowSums(S), 1))   # gene x stratum
  cat(sprintf("\n--- %s / %s: %s cells, %d FOVs, %d strata\n", SM, BR,
              format(n, big.mark=","), uniqueN(z$fov_id), length(sl)))

  for (bi in seq_len(length(BANDS) - 1L)) {
    lo <- BANDS[bi]; hi <- BANDS[bi + 1L]
    sel <- ND > lo & ND <= hi
    if (sum(sel) < MIN_PAIRS) next
    src <- row(NI)[sel]; tgt <- NI[sel]
    A   <- sparseMatrix(i = src, j = tgt, x = 1, dims = c(n, n))
    npair <- sum(A)
    ## OBS and EXP for every pair at once
    PA  <- P %*% A                                  # gene x cell
    OBS <- as.matrix(PA %*% t(P)) / npair           # gene x gene
    SAS <- as.matrix(S %*% A %*% t(S))              # stratum x stratum (small)
    EXP <- as.matrix(Pbar %*% SAS %*% t(Pbar)) / npair
    npL <- Matrix::rowSums(P)
    keep <- PR[npL[gi[ligand]] >= MIN_POS & npL[gi[receptor]] >= MIN_POS]
    if (!nrow(keep)) next
    o <- OBS[cbind(gi[keep$ligand], gi[keep$receptor])]
    e <- EXP[cbind(gi[keep$ligand], gi[keep$receptor])]
    stage1[[length(stage1)+1]] <- data.table(
      samp = SM, bin_region = BR, band = sprintf("%d-%d", lo, hi), band_lo = lo,
      ligand = keep$ligand, receptor = keep$receptor, tier = keep$tier,
      n_pairs = npair, n_Lpos = npL[gi[keep$ligand]], n_Rpos = npL[gi[keep$receptor]],
      obs = o, exp = e, exp_count = e * npair,
      sat = sat[[sprintf("<%d", hi)]],
      truncated = sat[[sprintf("<%d", hi)]] > MAX_SAT,
      enrich = fifelse(e * npair >= MIN_EXP, o / pmax(e, 1e-12), NA_real_))
    ## per-FOV aggregates for the bootstrap, on a reduced gene set
    G1 <- unique(c(PR[tier == "T1_declared", ligand], PR[tier == "T1_declared", receptor]))
    keepg <- gi[intersect(G1, GEN)]
    ## Per-FOV aggregates, built from the flat (src, tgt) vectors rather than by subsetting
    ## the big sparse adjacency twice per FOV -- that subsetting dominated the runtime.
    fsrc <- z$fov_id[src]; ftgt <- z$fov_id[tgt]
    within <- fsrc == ftgt
    Pk <- P[keepg, , drop = FALSE]
    for (fv in unique(fsrc[within])) {
      m2 <- within & fsrc == fv
      npf <- sum(m2)
      if (npf < 20L) next
      s2 <- src[m2]; t2 <- tgt[m2]
      jj <- which(z$fov_id == fv)
      Af <- sparseMatrix(i = match(s2, jj), j = match(t2, jj), x = 1,
                         dims = c(length(jj), length(jj)))
      Pf <- Pk[, jj, drop = FALSE]; Sf <- S[, jj, drop = FALSE]
      Pbf <- (Pf %*% t(Sf)) %*% Diagonal(x = 1/pmax(Matrix::rowSums(Sf), 1))
      fovagg[[length(fovagg)+1]] <- list(
        samp = SM, bin_region = BR, band = sprintf("%d-%d", lo, hi), fov = fv,
        npair = npf, OBS = as.matrix(Pf %*% Af %*% t(Pf)),
        EXP = as.matrix(Pbf %*% as.matrix(Sf %*% Af %*% t(Sf)) %*% t(Pbf)),
        genes = names(keepg))
    }
  }
}

R <- rbindlist(stage1)
if (!nrow(R)) { cat("\nnothing testable\n"); quit(status = 0) }
cat(sprintf("\n\nstage 1 complete: %d (region x band x pair) rows\n", nrow(R)))

## ---- stage 2: FOV cluster bootstrap on Tier 1 -------------------------------------------
boot_one <- function(SM, BR, BD, Lg, Rc) {
  L <- Filter(function(x) x$samp == SM && x$bin_region == BR && x$band == BD, fovagg)
  if (length(L) < 6) return(NULL)
  gn <- L[[1]]$genes
  if (!(Lg %in% gn) || !(Rc %in% gn)) return(NULL)
  iL <- match(Lg, gn); iR <- match(Rc, gn)
  no <- vapply(L, function(x) x$OBS[iL, iR], 0)
  ne <- vapply(L, function(x) x$EXP[iL, iR], 0)
  np <- vapply(L, function(x) x$npair, 0)
  bs <- replicate(NBOOT, {
    k <- sample(seq_along(L), length(L), TRUE)
    (sum(no[k]) / sum(np[k])) / max(sum(ne[k]) / sum(np[k]), 1e-12)
  })
  bs <- bs[is.finite(bs)]
  if (!length(bs)) return(NULL)
  data.table(fovs = length(L),
             enrich_fov = (sum(no)/sum(np)) / max(sum(ne)/sum(np), 1e-12),
             lo = quantile(bs, .025), hi = quantile(bs, .975),
             p_boot = 2*min((1+sum(bs <= 1))/(1+length(bs)), (1+sum(bs >= 1))/(1+length(bs))))
}
T1R <- R[tier == "T1_declared"]
bt <- rbindlist(lapply(seq_len(nrow(T1R)), function(i) {
  o <- boot_one(T1R$samp[i], T1R$bin_region[i], T1R$band[i], T1R$ligand[i], T1R$receptor[i])
  if (is.null(o)) return(NULL)
  cbind(T1R[i, .(samp, bin_region, band, band_lo, ligand, receptor, n_pairs, enrich)], o)
}), fill = TRUE)
if (nrow(bt)) bt[, q := p.adjust(p_boot, "BH"), by = .(samp, bin_region, band)]

cat("\n=== TIER 1 (declared axes), analytic enrichment and FOV cluster bootstrap\n")
cat("    `enrich` uses ALL band pairs including cross-FOV ones; `enrich_fov` uses only\n")
cat("    WITHIN-FOV pairs, which is what the bootstrap resamples. They are different\n")
cat("    quantities and the inference belongs to enrich_fov.\n")
if (nrow(bt)) print(bt[order(samp, bin_region, ligand, band_lo),
   .(samp, bin_region, ligand, receptor, band, n_pairs, fovs,
     enrich = round(enrich,3), enrich_fov = round(enrich_fov,3),
     ci = sprintf("[%.2f, %.2f]", lo, hi), p = signif(p_boot,3), q = signif(q,3))]) else
   cat("  no Tier 1 pair had enough FOVs to bootstrap\n")

## ---- the leakage discriminator ----------------------------------------------------------
prof <- dcast(R, samp + bin_region + ligand + receptor + tier ~ band, value.var = "enrich")
bcols <- setdiff(names(prof), c("samp","bin_region","ligand","receptor","tier"))
setnames(prof, bcols, paste0("e_", bcols))
shape <- function(d) {
  e1 <- d[["e_0-10"]]; e4 <- d[["e_40-80"]]
  fifelse(!is.finite(e1) | !is.finite(e4), "incomplete (a band lacked expected pairs)",
   fifelse(e4 >= 1.10, "persists to 40-80um: NOT leakage-shaped",
    fifelse(e1 > e4 * 1.5, "first band only: LEAKAGE-SHAPED", "intermediate")))
}
prof[, shape := shape(.SD)]
cat("\n=== band truncation (a truncated band is reported but not used for the shape call)\n")
print(R[, .(rows = .N, med_sat = round(median(sat),3), truncated = sum(truncated)),
        by = .(samp, bin_region, band)][order(samp, bin_region, band)])
cat("\n=== rows dropped for too few expected co-positive pairs (obs/exp unstable)\n")
print(R[, .(rows = .N, no_enrich = sum(is.na(enrich)),
            pct = round(100*mean(is.na(enrich)),1)), by = tier])

cat("\n=== leakage-shape classification, all pairs with a full profile\n")
print(prof[, .N, by = .(tier, shape)][order(tier, -N)])
cat("\n=== Tier 1 profiles across the distance bands\n")
print(prof[tier == "T1_declared"][order(samp, bin_region, ligand),
   c("samp","bin_region","ligand","receptor", paste0("e_", bcols), "shape"), with = FALSE])

## Tier 2: screen then report, corrected across everything
R2 <- R[tier == "T2_database"]
R2[, q_screen := NA_real_]
top2 <- merge(R2[!is.na(enrich) & enrich >= SCREEN_E & !truncated], prof[, .(samp,bin_region,ligand,receptor,shape)],
              by = c("samp","bin_region","ligand","receptor"), all.x = TRUE)
cat(sprintf("\n=== TIER 2: %d rows, %d at enrichment >= %.2f, of which not leakage-shaped: %d\n",
            nrow(R2), nrow(top2), SCREEN_E,
            top2[shape == "persists to 40-80um: NOT leakage-shaped", .N]))
if (nrow(top2)) print(top2[shape == "persists to 40-80um: NOT leakage-shaped"][
  order(-enrich)][1:min(25,.N), .(samp, bin_region, ligand, receptor, band,
  n_pairs, enrich = round(enrich,3))])
cat("\n  Tier 2 carries no bootstrap and therefore no p-value: it is a screen, not a test.\n")

sc("lr_proximity_stage1.csv.gz", R)
sc("lr_proximity_profile.csv", prof)
if (nrow(bt)) sc("lr_proximity_tier1.csv", bt)
if (nrow(top2)) sc("lr_proximity_tier2_screen.csv", top2)
writeLines(format(Sys.time()), file.path(OUT, "lr_proximity.done"))
cat("\nOK\n")
