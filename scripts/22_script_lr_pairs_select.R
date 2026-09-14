#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 22c -- pick the cell pairs to photograph, TWO WAYS on purpose.
##
## WHY LOOK AT ALL. Step 22 found ligand-receptor co-expression in proximity at 1.01-1.18
## fold. An enrichment that small either looks like something in tissue or it does not, and a
## number that size cannot be judged from the number alone. So: tight crops of real cell
## pairs, on real morphology, with the actual transcript positions drawn.
##
## TWO SELECTION RULES, BECAUSE THEY ANSWER DIFFERENT QUESTIONS AND ONLY ONE IS HONEST ABOUT
## TYPICALITY.
##
##   `cherry`  the clearest available examples -- ranked by min(ligand count, receptor count)
##             so BOTH sides are unambiguous, then by total counts. This shows what the
##             phenomenon looks like when it is unmistakable. It is CHERRY-PICKED and is
##             labelled as such in every output; it says nothing about how often it happens.
##
##   `repr`    a uniform random draw from ALL qualifying pairs, capped at two per field of
##             view so nine panels cannot all come from one lucky field. This is what a
##             typical "pair" behind a 1.06x enrichment actually looks like, and it is the
##             variant to trust when asking whether the enrichment is biology or noise.
##
## Both are drawn from the SAME candidate set, in the region and distance band where that
## pair was strongest in step 22 -- i.e. where the claim actually lives. Positivity is >= 2
## counts, the project convention throughout.
##
## THE PDCD1LG2-PDCD1 PANEL IS THE POINT OF THE EXERCISE. Step 22 measured it at 2.28x in the
## 0-10 um band of IB1's FL foci, decaying to 1.01 by 40-80 um -- the signature of transcript
## leakage rather than of biology. At 0-10 um centre-to-centre, two "cells" are closer than
## one cell diameter, so the crops should show directly whether they are two cells or one cell
## split by segmentation. That is a thing to see, not to infer.
##
## A cell that is positive for BOTH genes is excluded as a partner to itself.
##
## Usage: sbatch --mem=48G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/22_script_lr_pairs_select.R [tier]
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(Matrix); library(RANN); library(data.table)})
setDTthreads(1); options(width = 210)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
LR   <- "output/22_lr"; OUT <- file.path(LR, TIER, "snapshots")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

POS_MIN <- 2L; N_PICK <- 9L; MAX_PER_FOV <- 2L; SEED <- 22L
## `Hodgkin chemokine niche` is a REGION label, not a lineage (CLAUDE.md), and
## 22_script_lr_proximity.R excludes it from every cell-type role. The first version of this
## script did not, so 16 picked panels had a focal cell whose "type" was a region -- and
## CCL22-CCR4's representative draw was 7 of 9 such cells. The snapshots have to come from
## the same cell set as the result they illustrate.
NICHE <- "Hodgkin chemokine niche"
Q_LIKELY <- 0.05; E_NEAR <- 1.25

M  <- readRDS(file.path(LR, "lr_counts.Rds"))
CE <- fread(file.path(LR, "lr_cells.csv.gz"))
setkey(CE, cell_id); CE <- CE[colnames(M)]
RG <- fread(file.path("output/16_regions", TIER, "cell_regions.csv.gz"),
            select = c("cell_id","bin_region"))
CE <- merge(CE, RG, by = "cell_id", all.x = TRUE, sort = FALSE)
setkey(CE, cell_id); CE <- CE[colnames(M)]
stopifnot(identical(CE$cell_id, colnames(M)))

## ---- which pairs, and in which region/band, from step 22's own table -------------------
B <- fread(file.path(LR, TIER, "lr_proximity_tier1.csv"))
B[, pair := paste0(ligand, "-", receptor)]
best <- B[order(-enrich_fov)][, .SD[1], by = .(samp, bin_region, pair)]
sel <- best[(q < Q_LIKELY & enrich_fov > 1) | enrich_fov >= E_NEAR]
sel[, status := fifelse(q < Q_LIKELY & enrich_fov > 1, "likely", "near-likely")]
## one row per pair: the region/band where it is strongest
sel <- sel[order(-enrich_fov)][, .SD[1], by = pair]
sel[, c("band_lo","band_hi") := tstrsplit(band, "-", type.convert = TRUE)]
cat("=== pairs to photograph, each in the region/band where step 22 found it strongest\n")
print(sel[, .(pair, status, samp, bin_region, band, fovs, enrich = round(enrich_fov,3),
               q = signif(q,3))])
sc("snapshot_targets.csv", sel[, .(pair, ligand, receptor, status, samp, bin_region, band,
                                   band_lo, band_hi, enrich_fov, q)])

gi <- setNames(seq_len(nrow(M)), rownames(M))
set.seed(SEED)
picks <- list(); cand_stats <- list()

for (i in seq_len(nrow(sel))) {
  Lg <- sel$ligand[i]; Rc <- sel$receptor[i]
  SM <- sel$samp[i];   BR <- sel$bin_region[i]
  lo <- sel$band_lo[i]; hi <- sel$band_hi[i]
  z  <- CE[samp == SM & bin_region == BR & !is.na(celltype_cluster) &
           celltype_cluster != NICHE]
  if (!nrow(z)) next
  ii <- match(z$cell_id, colnames(M))
  lc <- as.numeric(M[gi[Lg], ii]); rc <- as.numeric(M[gi[Rc], ii])
  lp <- which(lc >= POS_MIN); rp <- which(rc >= POS_MIN)
  if (length(lp) < 5L || length(rp) < 5L) next

  ## every R+ cell within (lo, hi] um of each L+ cell
  xyL <- as.matrix(z[lp, .(x_slide_mm, y_slide_mm)]) * 1000
  xyR <- as.matrix(z[rp, .(x_slide_mm, y_slide_mm)]) * 1000
  k   <- min(40L, length(rp))
  nn  <- RANN::nn2(xyR, xyL, k = k)
  keep <- nn$nn.dists > lo & nn$nn.dists <= hi
  if (!any(keep)) next
  src <- lp[row(nn$nn.idx)[keep]]
  tgt <- rp[nn$nn.idx[keep]]
  d   <- nn$nn.dists[keep]
  ok  <- src != tgt                      # a both-positive cell is not its own partner
  src <- src[ok]; tgt <- tgt[ok]; d <- d[ok]
  if (!length(src)) next

  C <- data.table(
    pair = sel$pair[i], ligand = Lg, receptor = Rc, status = sel$status[i],
    samp = SM, bin_region = BR, band = sel$band[i],
    L_cell = z$cell_id[src], R_cell = z$cell_id[tgt],
    L_count = lc[src], R_count = rc[tgt],
    L_type = z$celltype_cluster[src], R_type = z$celltype_cluster[tgt],
    L_conf = z$confidence[src], R_conf = z$confidence[tgt],
    L_fov = z$fov_id[src], R_fov = z$fov_id[tgt],
    L_x_mm = z$x_slide_mm[src], L_y_mm = z$y_slide_mm[src],
    R_x_mm = z$x_slide_mm[tgt], R_y_mm = z$y_slide_mm[tgt],
    dist_um = d)
  C[, `:=`(mid_x_mm = (L_x_mm + R_x_mm)/2, mid_y_mm = (L_y_mm + R_y_mm)/2,
           both_min = pmin(L_count, R_count), both_sum = L_count + R_count,
           same_fov = L_fov == R_fov)]
  C <- C[same_fov == TRUE]               # a crop must sit inside one field of view
  if (!nrow(C)) next

  cand_stats[[length(cand_stats)+1]] <- C[, .(pair = pair[1], status = status[1],
    samp = samp[1], bin_region = bin_region[1], band = band[1],
    candidates = .N, fovs = uniqueN(L_fov),
    med_dist = round(median(dist_um),1),
    med_L = median(L_count), med_R = median(R_count),
    max_both_min = max(both_min))]

  ## cap per FOV, then apply each rule
  take_capped <- function(D, n) {
    D[, .SD[seq_len(min(.N, MAX_PER_FOV))], by = L_fov][seq_len(min(.N, n))]
  }
  cherry <- take_capped(C[order(-both_min, -both_sum, dist_um)], N_PICK)
  repr   <- take_capped(C[sample(.N)], N_PICK)
  if (nrow(cherry)) picks[[length(picks)+1]] <- cbind(variant = "cherry",
                          panel = seq_len(nrow(cherry)), cherry)
  if (nrow(repr))   picks[[length(picks)+1]] <- cbind(variant = "repr",
                          panel = seq_len(nrow(repr)), repr)
}

CS <- rbindlist(cand_stats); P <- rbindlist(picks, fill = TRUE)
cat("\n=== candidate pool per LR pair (what the nine are drawn FROM)\n"); print(CS)
sc("snapshot_candidate_pool.csv", CS)

cat("\n=== picked panels, per pair and variant\n")
print(P[, .(panels = .N, fovs = uniqueN(L_fov),
            dist_um = paste0(round(min(dist_um)), "-", round(max(dist_um))),
            med_both_min = median(both_min)), by = .(pair, status, variant)][order(pair, variant)])
cat("\n  `cherry` is ranked by min(L,R) counts: the clearest examples, NOT typical.\n")
cat("  `repr` is a uniform draw from the same pool: typical, and the one to judge on.\n")

sc("snapshot_picks.csv", P)
## FOVs the transcript extraction needs -- keeping it to these is what makes 22d affordable
fv <- unique(P[, .(samp, fov = L_fov)])
sc("snapshot_fovs.csv", fv)
cat(sprintf("\nFOVs needed for the transcript extract: IB1 %d, RB1 %d\n",
            fv[samp=="IB1", .N], fv[samp=="RB1", .N]))
cat("genes needed:", paste(sort(unique(c(P$ligand, P$receptor))), collapse=", "), "\n")
writeLines(format(Sys.time()), file.path(OUT, "pairs_select.done"))
cat("\nOK\n")
