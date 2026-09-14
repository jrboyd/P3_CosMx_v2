#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 24 -- DOES THE ARCHITECTURE RESULT SURVIVE THE BETTER CELL LABELS?  The G1 gap.
##
## THE GAP. Step 19 built the project's best cell identity and nothing consumed it: every
## cell-type-conditioned result in steps 16 and 17 still reads
## output/16_celltypes/cell_celltypes.csv.gz. That matters most for step 16's co-occurrence
## result -- "IB1's node is spatially structured; RB1's is not" -- because the review called
## it the one place where "differently organised" is demonstrated rather than inferred, and
## because the substrate it ran on is largely unnamed: the step-16 argmax leaves 2.18M of
## 4.4M cells in `ambiguous`, `all_markers_0` or `below_background`, against step 19's 2.05M
## named cells at the matched tier.
##
## WHAT IS AND IS NOT RE-TESTABLE, AND THIS IS NOT A LIMITATION OF THIS SCRIPT.
## Step 16's headline numbers are about **Tfh** (self-aggregation 2.74x, Tfh-B 1.62x,
## Tfh-macrophage 0.47x). **Step 19 has no Tfh label, on purpose**: cluster-based T-cell
## subsetting failed three ways and is a closed avenue, so T stays a single label and Tfh
## comes only from the per-cell argmax subsets -- which is what step 16 already does. So the
## Tfh-specific findings CANNOT be re-derived from step-19 labels, and neither can step 16's
## Tfh composition result. Any claim that they were "confirmed" here would be false.
##
## What IS testable is the part that does not depend on a T subset: the coarse architecture
## contrast. Step 16 reports lymphoid/myeloid segregation running 0.48-0.57x throughout IB1
## against 0.93-1.09 (chance) everywhere in RB1. That is a statement about B, T, myeloid,
## NK, stromal, endothelial and plasma cells, all seven of which exist in both vocabularies.
##
## THE DESIGN. The identical permutation statistic -- mean count of type-B neighbours within
## R_UM of a type-A cell, labels permuted within local-density deciles so composition and
## geometry are held fixed -- run TWICE on the SAME cells and the SAME regions, once per
## label source, and compared pair by pair. Agreement is the result; disagreement localises
## which labels carry the effect.
##
## TWO STEP-19 LABELS ARE HANDLED SPECIALLY.
##   `Hodgkin chemokine niche` is EXCLUDED. CLAUDE.md is explicit that it is a REGION, not a
##   cell type (is_niche = TRUE), and is barred from cell-type statements. Including it would
##   put 82,228 cells into a co-occurrence table as though they were a lineage.
##   `pDC` has no step-16 counterpart, so it is reported for step-19 only, never compared.
##
## Cells are the matched tier, non-defect, step-19-labelled. Region comes from step 16 so the
## two runs are compared on identical tissue.
##
## Usage: sbatch --mem=64G --time=2:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/24_script_arch_relabel.R [tier]
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(RANN); library(data.table)})
setDTthreads(1); options(width = 200)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
REG_F <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT16  <- "output/16_celltypes/cell_celltypes.csv.gz"
CT19  <- "output/19_celltypes/all_cells_celltypes.csv.gz"
OUT   <- file.path("output/24_arch_relabel", TIER)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }
for (f in c(REG_F, CT16, CT19)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)

## Constants lifted verbatim from 16_script_cooccurrence.R so the two runs are comparable.
R_UM   <- 20; K_MAX <- 40L; N_PERM <- 200L; MIN_A <- 100L; N_DEC <- 10L; SEED <- 24L
NICHE_LABEL <- "Hodgkin chemokine niche"

## step-19 -> the shared coarse vocabulary. Subsets collapse the same way step 16's coarse
## lineage does, so "T" means the same population in both columns.
MAP19 <- c("B cell"="B", "T cell"="T", "Myeloid"="myeloid", "Stromal/FDC"="fibro",
           "Endothelial"="endo", "NK cell"="NK", "Plasma cell"="PC",
           "Plasmablast"="PC", "pDC"="pDC")
MAP16 <- c("B"="B", "T_unresolved"="T", "Tfh"="T", "Treg"="T", "CD8_T"="T", "CD4_conv"="T",
           "Monocyte"="myeloid", "Macrophage"="myeloid", "cDC"="myeloid",
           "myeloid_unresolved"="myeloid", "fibro"="fibro", "endo"="endo",
           "NK"="NK", "PC"="PC")
SHARED <- c("B","T","myeloid","fibro","endo","NK","PC")

cat("========================================================================\n")
cat("step 24 -- architecture under step-16 vs step-19 labels, tier =", TIER, "\n")
cat("========================================================================\n")

md <- fread(REG_F, select = c("cell_id","samp","fov_id","x_slide_mm","y_slide_mm",
                              "bin_region","is_defect","tier_matched"))
md <- md[tier_matched == TRUE & is_defect == FALSE & bin_region != "unbinned"]
c16 <- fread(CT16, select = c("cell_id","celltype"))
c19 <- fread(CT19, select = c("cell_id","celltype_cluster","confidence"))
md <- merge(merge(md, c16, by = "cell_id", all.x = TRUE), c19, by = "cell_id", all.x = TRUE)

md[, lab16 := MAP16[celltype]]
md[, lab19 := fifelse(celltype_cluster == NICHE_LABEL, NA_character_, MAP19[celltype_cluster])]
cat(sprintf("\ncells: %s | named by step 16: %s | by step 19: %s | by both: %s\n",
  format(nrow(md), big.mark=","), format(md[!is.na(lab16), .N], big.mark=","),
  format(md[!is.na(lab19), .N], big.mark=","),
  format(md[!is.na(lab16) & !is.na(lab19), .N], big.mark=",")))
cat(sprintf("step-19 `%s` cells excluded as a region, not a lineage: %s\n",
  NICHE_LABEL, format(md[celltype_cluster == NICHE_LABEL, .N], big.mark=",")))

cat("\n=== agreement between the two label sources, on cells named by both\n")
xt <- md[!is.na(lab16) & !is.na(lab19) & lab19 %in% SHARED & lab16 %in% SHARED,
         .N, by = .(lab16, lab19)]
tot <- xt[, sum(N)]
cat(sprintf("  overall agreement: %.1f%% of %s cells\n",
            100*xt[lab16 == lab19, sum(N)]/tot, format(tot, big.mark=",")))
ag <- xt[, .(cells = sum(N), agree_pct = round(100*sum(N[lab16 == lab19])/sum(N), 1)), by = lab19]
print(ag[order(-cells)]); sc("label_agreement.csv", xt)

## ---- the permutation statistic, verbatim in structure from 16_script_cooccurrence.R -----
run_one <- function(z, labcol) {
  z <- z[!is.na(get(labcol))]
  if (nrow(z) < 3*MIN_A) return(NULL)
  xy <- as.matrix(z[, .(x_slide_mm, y_slide_mm)]) * 1000     # mm -> um
  k  <- min(K_MAX, nrow(z) - 1L)
  nn <- RANN::nn2(xy, xy, k = k + 1L)
  idx <- nn$nn.idx[, -1, drop = FALSE]; dst <- nn$nn.dists[, -1, drop = FALSE]
  keep <- dst <= R_UM
  ## density stratum = distance to the 20th neighbour, in deciles
  dcol <- min(20L, ncol(dst))
  dens <- dst[, dcol]
  strat <- cut(dens, quantile(dens, seq(0, 1, length.out = N_DEC + 1L)),
               include.lowest = TRUE, labels = FALSE)
  lab0 <- z[[labcol]]
  types <- intersect(c(SHARED, "pDC"), unique(lab0))
  ti <- setNames(seq_along(types), types)
  stat <- function(lab) {
    li <- ti[lab]
    m <- matrix(0, length(types), length(types), dimnames = list(types, types))
    for (j in seq_len(ncol(idx))) {
      ok <- keep[, j]
      if (!any(ok)) next
      a <- li[ok]; b <- li[idx[ok, j]]
      good <- !is.na(a) & !is.na(b)
      if (!any(good)) next
      m <- m + as.matrix(Matrix::sparseMatrix(i = a[good], j = b[good], x = 1,
             dims = c(length(types), length(types))))
    }
    na <- as.numeric(table(factor(lab, levels = types)))
    sweep(m, 1, pmax(na, 1), "/")
  }
  obs <- stat(lab0)
  set.seed(SEED)
  perms <- lapply(seq_len(N_PERM), function(p) {
    lab <- lab0
    for (s in unique(strat[!is.na(strat)])) {
      ii <- which(strat == s); lab[ii] <- sample(lab0[ii])
    }
    stat(lab)
  })
  nu <- Reduce(`+`, perms) / N_PERM
  ge <- Reduce(`+`, lapply(perms, function(m) m >= obs))
  na <- as.numeric(table(factor(lab0, levels = types)))
  data.table(a = rep(types, times = length(types)),
             b = rep(types, each = length(types)),
             n_a = rep(na, times = length(types)),
             n_b = rep(na, each = length(types)),
             obs = as.vector(obs), null = as.vector(nu),
             p_emp = (1 + as.vector(ge)) / (1 + N_PERM))[
             , enrich := obs / pmax(null, 1e-9)][]
}

REGIONS <- md[, .N, by = .(samp, bin_region)][N >= 3*MIN_A][order(samp, -N)]
cat("\nregions run:\n"); print(REGIONS)

res <- rbindlist(lapply(seq_len(nrow(REGIONS)), function(i) {
  SM <- REGIONS$samp[i]; BR <- REGIONS$bin_region[i]
  z <- md[samp == SM & bin_region == BR]
  rbindlist(lapply(c("lab16","lab19"), function(L) {
    o <- run_one(z, L)
    if (is.null(o)) return(NULL)
    cbind(samp = SM, bin_region = BR, labels = L, o)
  }), fill = TRUE)
}), fill = TRUE)

res <- res[n_a >= MIN_A & n_b >= MIN_A]
res[, q := p.adjust(p_emp, "BH"), by = .(samp, bin_region, labels)]
sc("cooccurrence_both_labels.csv", res)

## ---- do the two label sources agree about architecture? ---------------------------------
W <- dcast(res[a %in% SHARED & b %in% SHARED],
           samp + bin_region + a + b ~ labels, value.var = c("enrich","q"))
W <- W[!is.na(enrich_lab16) & !is.na(enrich_lab19)]
cat(sprintf("\n=== %d (region x pair) cells testable under BOTH label sources\n", nrow(W)))
cat(sprintf("  Spearman rho of enrichment, step16 vs step19: %.3f\n",
            cor(W$enrich_lab16, W$enrich_lab19, method = "spearman")))
cat(sprintf("  sign agreement about chance (both >1 or both <1): %.1f%%\n",
            100*mean((W$enrich_lab16 > 1) == (W$enrich_lab19 > 1))))
sc("enrichment_concordance.csv", W)

## The specific claim being re-tested: lymphoid/myeloid segregation in IB1, chance in RB1.
LYM <- c("B","T","NK"); MYE <- c("myeloid")
seg <- res[a %in% LYM & b %in% MYE | a %in% MYE & b %in% LYM]
S <- seg[, .(pairs = .N, med_enrich = round(median(enrich), 3),
             min_enrich = round(min(enrich), 3), max_enrich = round(max(enrich), 3),
             n_sig = sum(q < 0.05)), by = .(samp, bin_region, labels)][order(samp, bin_region, labels)]
cat("\n=== lymphoid <-> myeloid segregation, both label sources\n")
print(S); sc("segregation_by_labels.csv", S)

allp <- res[a %in% SHARED & b %in% SHARED,
            .(pairs = .N, med_enrich = round(median(enrich), 3),
              frac_at_chance = round(mean(q >= 0.05), 3)), by = .(samp, labels)]
cat("\n=== overall: how much of each sample sits at chance?\n")
print(allp); sc("chance_fraction.csv", allp)

cat("\n=== step-19-only label (no step-16 counterpart; reported, never compared)\n")
print(res[labels == "lab19" & (a == "pDC" | b == "pDC") & q < 0.05,
          .(samp, bin_region, a, b, enrich = round(enrich,3), q = signif(q,3))][1:min(15,.N)])

writeLines(format(Sys.time()), file.path(OUT, "arch_relabel.done"))
cat("\nOK\n")
