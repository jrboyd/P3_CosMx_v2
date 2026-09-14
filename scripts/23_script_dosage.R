#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 23 -- IS THERE DOSAGE EVIDENCE FOR 9p24.1 GAIN?  The G5 gap.
##
## WHY THIS EXISTS. Step 17's strongest within-sample result is that IB1's Hodgkin hotspots
## carry elevated PD-L1 and PD-L2 with PDCD1 and the IFN-gamma response FLAT -- ligand-side
## without interferon, which is the pattern expected from 9p24.1 copy gain rather than from
## adaptive immune resistance. 9p24.1 gain is a defining lesion of classic Hodgkin lymphoma.
## But step 17's "9p24.1 module" is CD274 + PDCD1LG2 + JAK2, and it is labelled EXPLORATORY
## for a good reason: all three are immune-functional, so an inflammatory stimulus could raise
## them together and RNA alone cannot tell that from a copy-number event.
##
## THE DISCRIMINATOR THIS STEP ADDS, AND IT NEEDED AN ANNOTATION PACKAGE TO SEE.
## org.Hs.eg.db (installed for this step) shows the panel carries NINE genes on 9p24.1, not
## three:
##      immune / IFN-inducible   CD274, PDCD1LG2, JAK2
##      cytokine                 IL33
##      NOT immune-functional    CDC37L1, ERMP1, KDM4C, MLANA, PLPP6
## Copy gain is indifferent to what a gene does: it raises everything in the interval. IFN or
## NF-kB signalling raises the immune genes and leaves CDC37L1, ERMP1, KDM4C, MLANA and PLPP6
## alone. So the load-bearing test is the score on the NON-IMMUNE subset, and the immune
## subset is the circular one that must be reported separately.
##
## THE NULL IS THE REST OF THE GENOME, WHICH IS WHY A TARGETED PANEL IS SURVIVABLE. A 6K
## panel does not sample the genome uniformly and is heavily immune-biased, so an absolute
## 9p24.1 score means nothing. What means something is its RANK against every other cytoband
## measured on the same panel, in the same cells, against the same reference: 444 bands carry
## >= 5 panel genes and 41 chromosome arms carry >= 20. Panel-selection bias applies to the
## null exactly as it applies to 9p24.1.
##
## THE CELL-INTRINSIC CONSTRAINT, WHICH IS HOW TO READ THE PER-CELL-TYPE TABLE. Copy gain
## belongs to a CLONE. If IB1's hotspots carry 9p24.1 gain it is in the HRS cells, which this
## panel cannot call. So an elevated 9p24.1 score in the hotspots' T cells, fibroblasts or B
## cells is NOT evidence of gain in those cells -- they are not the clone. Two other things
## produce it: transcript leakage from a neighbouring high-expressing cell (step 18 measured
## exactly this displacement, and step 19 showed the cluster caller is the more contaminated
## of the two), and a genuine regional effect on expression that has nothing to do with
## dosage. A dosage reading requires the signal to be CONFINED, not spread across bystander
## lineages; a score that rises in several unrelated lineages at once is the signature of
## leakage or of a regional gradient, and is reported as such.
##
## WHAT THIS CANNOT DO, STATED UP FRONT. RNA cannot establish amplification. A high 9p24.1
## rank is consistent with copy gain and does not demonstrate it; a low rank is the more
## informative outcome, because it removes the one mechanism step 17 proposed. This is a
## dosage-CONSISTENCY test, not CNV inference -- no HMM, no segmentation, no allele
## information.
##
## CONVENTIONS OBSERVED. Region is the unit (FOV cluster bootstrap). Cell-type conditioned,
## because the hotspots are DC-rich and a band score pooled over cell types would partly
## measure composition. cp10k is compositional, so every score is a ratio against a
## within-sample reference region and the concordance with detection rate is reported.
## IB1's reference is `background` restricted to > 800 um from any hotspot, per CLAUDE.md
## fact 3. Defect burden is handled by a low-burden sensitivity re-run rather than a
## covariate, since a band score is already a ratio and the defect flattens numerator and
## null alike (CLAUDE.md fact 9).
##
## Usage: sbatch --mem=64G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/23_script_dosage.R [tier]
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(Matrix); library(data.table)})
setDTthreads(1); options(width = 200)

AG    <- commandArgs(TRUE)
TIER  <- if (length(AG) >= 1) AG[1] else "matched"
SMOKE <- length(AG) >= 2 && AG[2] == "smoke"
PB   <- file.path("output/16_ct_pseudobulk", TIER)
IF21 <- file.path("output/21_interface", TIER)
OUT  <- if (SMOKE) file.path("output/23_dosage", TIER, "smoke") else file.path("output/23_dosage", TIER)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

BAND_TARGET <- "9p24.1"
IMMUNE_9P   <- c("CD274","PDCD1LG2","JAK2")      # IFN/NF-kB inducible -- the circular subset
CYTOKINE_9P <- c("IL33")                          # secreted; excluded from the clean subset
MIN_BAND_N  <- 5L                                 # genes needed for a band to enter the null
MIN_DET     <- 0.02                               # per-gene pooled detection floor
MIN_CELLS   <- 5L
CT_TESTED   <- c("B","T","myeloid","fibro")
NBOOT       <- 2000L
if (exists("SMOKE") && SMOKE) NBOOT <- 25L
FAR_UM      <- 800                                # clean IB1 reference, CLAUDE.md fact 3
SEED        <- 23L

cat("========================================================================\n")
if (SMOKE) cat("*** SMOKE MODE -> ", OUT, " ***\n", sep = "")
cat("step 23 -- 9p24.1 dosage consistency, tier =", TIER, "\n")
cat("========================================================================\n")

CB <- fread("output/23_dosage/panel_cytoband.csv")
S  <- readRDS(file.path(PB, "pb_binct_sum.Rds"))
D  <- readRDS(file.path(PB, "pb_binct_det.Rds"))
GB <- fread(file.path(PB, "binct_meta.csv.gz")); setkey(GB, grp_bin); GB <- GB[colnames(S)]
stopifnot(identical(GB$grp_bin, colnames(S)))
GENES <- rownames(S)
GB[, row := .I]

## Clean IB1 reference: background bins > FAR_UM from any hotspot. bin_axes carries d_hot_um
## for every non-hotspot bin (step 21a).
AX <- fread(file.path(IF21, "bin_axes.csv.gz"))
AX[, bin := paste0("IB1_", bx, "_", by_)]
GB <- merge(GB, AX[, .(bin, d_hot_um)], by = "bin", all.x = TRUE)
GB[, zone := fifelse(samp == "RB1", bin_region,
              fifelse(bin_region == "background" & !is.na(d_hot_um) & d_hot_um > FAR_UM,
                      "far_background", bin_region))]
GB <- GB[cells >= MIN_CELLS]
cat("\ngroups by sample and zone:\n")
print(GB[, .(groups = .N, bins = uniqueN(bin), fovs = uniqueN(fov)),
         by = .(samp, zone)][order(samp, -groups)])

## Contrasts: disease zone vs the SAME SAMPLE's reference zone. No cross-sample subtraction.
CONTRASTS <- list(
  list(samp="IB1", a="cHL_hotspot", b="far_background", lab="IB1 hotspot vs far background"),
  list(samp="IB1", a="FL_focus",    b="far_background", lab="IB1 FL focus vs far background"),
  list(samp="RB1", a="FL_follicular", b="interfollicular", lab="RB1 follicular vs interfollicular"))

## Per-group library size, computed ONCE. The first version recomputed rowSums over all
## 6,175 genes inside every bootstrap replicate purely to obtain the denominator, which is
## why the job hit its wall clock: 2,000 reps x 12 strata x a 6,175-row rowSums. The
## denominator is the column total, which does not depend on the gene subset, so it is
## precomputed and each replicate then touches only the handful of rows it actually needs.
COLTOT <- Matrix::colSums(S)

cp <- function(idx) {  # pooled cp10k over a set of groups, all genes
  s <- Matrix::rowSums(S[, idx, drop = FALSE]); s / max(sum(s), 1) * 1e4
}
cp_rows <- function(rows, idx) {  # pooled cp10k for a few genes only -- the bootstrap path
  num <- Matrix::rowSums(S[rows, idx, drop = FALSE])
  num / max(sum(COLTOT[idx]), 1) * 1e4
}

set.seed(SEED)
res <- list(); gene_tab <- list()

for (K in CONTRASTS) for (CT in CT_TESTED) {
  ga <- GB[samp == K$samp & zone == K$a & coarse == CT]
  gb <- GB[samp == K$samp & zone == K$b & coarse == CT]
  cat(sprintf("  [%s | %s] groups %d/%d  fovs %d/%d\n", K$lab, CT,
              nrow(ga), nrow(gb), uniqueN(ga$fov), uniqueN(gb$fov)))
  if (uniqueN(ga$fov) < 6L || uniqueN(gb$fov) < 6L) { cat("    -> skipped: too few FOVs\n"); next }
  ia <- ga$row; ib <- gb$row
  cpa <- cp(ia); cpb <- cp(ib)
  deta <- Matrix::rowSums(D[, ia, drop=FALSE]) / sum(ga$cells)
  detb <- Matrix::rowSums(D[, ib, drop=FALSE]) / sum(gb$cells)
  keep <- (deta >= MIN_DET | detb >= MIN_DET)
  lfc  <- log2((cpa + 0.1) / (cpb + 0.1))
  lfcd <- log2((deta + 0.002) / (detb + 0.002))

  G <- data.table(gene = GENES, lfc = lfc, lfc_det = lfcd, keep = keep)
  G <- merge(G, CB[, .(gene, cytoband, band, arm)], by = "gene")[keep == TRUE]

  ## band scores, and the empirical null across all bands with enough genes
  BS <- G[, .(n = .N, score = mean(lfc), score_det = mean(lfc_det)), by = cytoband][n >= MIN_BAND_N]
  cat(sprintf("    genes kept %d, bands >=%d genes %d, 9p24.1 present %s\n",
              nrow(G), MIN_BAND_N, nrow(BS), nrow(BS[cytoband == BAND_TARGET]) > 0))
  if (!nrow(BS[cytoband == BAND_TARGET])) { cat("    -> skipped: 9p24.1 not in the null set\n"); next }
  ## Rank BEFORE extracting the target row. `tgt <- BS[...]` takes a COPY, so adding
  ## rank_desc to BS afterwards left tgt$rank_desc NULL, pct zero-length, and every
  ## data.table() below a 0-ROW table -- which surfaced as "nothing to report" rather than
  ## as an error. data.table copy semantics, and a silent one.
  BS[, rank_desc := frank(-score, ties.method = "min")]
  tgt <- BS[cytoband == BAND_TARGET]
  stopifnot(nrow(tgt) == 1L, is.finite(tgt$rank_desc))
  pct <- 100 * (1 - (tgt$rank_desc - 1) / nrow(BS))

  ## the discriminator: the same score on three nested gene sets
  g9 <- G[cytoband == BAND_TARGET]
  sub <- function(v) if (!length(v)) NA_real_ else mean(v)
  imm   <- g9[gene %in% IMMUNE_9P]
  clean <- g9[!gene %in% c(IMMUNE_9P, CYTOKINE_9P)]

  ## FOV cluster bootstrap on the CLEAN subset score -- the load-bearing number
  fa <- split(ia, ga$fov); fb <- split(ib, gb$fov)
  cg <- clean$gene
  crow <- match(cg, GENES)
  bs <- replicate(NBOOT, {
    ka <- unlist(fa[sample(names(fa), length(fa), TRUE)], use.names = FALSE)
    kb <- unlist(fb[sample(names(fb), length(fb), TRUE)], use.names = FALSE)
    mean(log2((cp_rows(crow, ka) + 0.1) / (cp_rows(crow, kb) + 0.1)))
  })
  bs <- bs[is.finite(bs)]

  res[[length(res)+1]] <- data.table(
    contrast = K$lab, celltype = CT,
    fovs_a = uniqueN(ga$fov), fovs_b = uniqueN(gb$fov),
    bands_in_null = nrow(BS),
    n_9p24_1 = nrow(g9), score_all = tgt$score, score_all_det = tgt$score_det,
    rank_of = tgt$rank_desc, pct_of_bands = round(pct, 1),
    n_immune = nrow(imm),  score_immune = sub(imm$lfc),
    n_clean  = nrow(clean), score_clean = sub(clean$lfc),
    score_clean_det = sub(clean$lfc_det),
    clean_lo = if (length(bs)) quantile(bs, .025) else NA_real_,
    clean_hi = if (length(bs)) quantile(bs, .975) else NA_real_,
    clean_p  = if (length(bs)) 2*min((1+sum(bs<=0))/(1+length(bs)),
                                     (1+sum(bs>=0))/(1+length(bs))) else NA_real_,
    concordant = is.finite(sub(clean$lfc)) & is.finite(sub(clean$lfc_det)) &
                 sign(sub(clean$lfc)) == sign(sub(clean$lfc_det)))
  gene_tab[[length(gene_tab)+1]] <- cbind(contrast = K$lab, celltype = CT,
    g9[, .(gene, cytoband, lfc = round(lfc,3), lfc_det = round(lfc_det,3),
           class = fifelse(gene %in% IMMUNE_9P, "immune/IFN-inducible",
                    fifelse(gene %in% CYTOKINE_9P, "cytokine", "not immune-functional")))])
}

R <- rbindlist(res); GT <- rbindlist(gene_tab)
if (!nrow(R)) { cat("\nno contrast had enough FOVs; nothing to report\n"); quit(status = 0) }
R[, q_clean := p.adjust(clean_p, "BH")]
## ---- the verdict, and it MUST consult RB1 ----------------------------------------------
## RB1 is the pure-FL recurrence: 9p24.1 gain is a classic-Hodgkin lesion and should be
## absent there, so RB1 is the construction-identical negative control, exactly as it is for
## the chemokine field in step 17. A statistic that fires in RB1 is not measuring the lesion.
## The first version of this verdict did not look at RB1 and therefore labelled RB1's own
## fibroblasts "DOSAGE-CONSISTENT" -- shipping that label would have been misleading, so the
## control is now a precondition rather than a commentary.
##
## The verdict is a pure function of columns already saved, so it can be recomputed from
## dosage_9p24.csv without repeating the bootstrap.
RB1_FIRES <- R[grepl("^RB1", contrast) & clean_p < 0.05 & score_clean > 0, .N] > 0
R[, verdict := fifelse(!is.finite(clean_p), "unestimable",
              fifelse(q_clean < 0.05 & concordant & score_clean > 0 & !RB1_FIRES,
                      "DOSAGE-CONSISTENT: non-immune 9p24.1 genes also elevated",
              fifelse(q_clean < 0.05 & concordant & score_clean > 0 & RB1_FIRES,
                      "elevated, but the no-Hodgkin control fires too: NOT dosage evidence",
              fifelse(score_immune > 0 & score_clean <= 0 & q_clean < 0.05,
                      "IMMUNE-ONLY: consistent with induction, NOT with copy gain",
                      "no evidence either way"))))]
if (RB1_FIRES) {
  cat("\n*** NEGATIVE CONTROL FIRED ***\n")
  cat("The 9p24.1 clean-subset score is significantly elevated in RB1, which carries no\n")
  cat("Hodgkin component. The statistic is therefore not specific to the lesion, and no\n")
  cat("row in this table may be read as dosage evidence. See the confinement table below:\n")
  cat("a clonal copy gain would also have to be confined to ONE lineage, and it is not.\n")
}

cat("\n\n=== 9p24.1 band score against the whole-genome band null\n")
print(R[, .(contrast, celltype, fovs_a, bands = bands_in_null,
            score_all = round(score_all,3), rank = paste0(rank_of,"/",bands_in_null),
            top_pct = pct_of_bands)])
cat("\n=== THE DISCRIMINATOR: immune vs non-immune genes on the same band\n")
print(R[, .(contrast, celltype,
            n_imm = n_immune,  imm  = round(score_immune,3),
            n_cln = n_clean,   cln  = round(score_clean,3),
            ci = sprintf("[%.2f, %.2f]", clean_lo, clean_hi),
            p = signif(clean_p,3), q = signif(q_clean,3),
            conc = concordant, verdict)])
cat("\n=== how many lineages carry the signal? (confinement is what dosage predicts)\n")
conf <- R[, .(lineages_tested = .N,
              lineages_clean_up = sum(score_clean > 0, na.rm = TRUE),
              lineages_q05 = sum(q_clean < 0.05 & score_clean > 0, na.rm = TRUE)),
          by = contrast]
print(conf)
cat("  A contrast with the score up in SEVERAL unrelated lineages is leakage- or\n")
cat("  gradient-shaped, not dosage-shaped: copy gain belongs to one clone.\n")
sc("dosage_confinement.csv", conf)

cat("\n=== per-gene detail on 9p24.1\n")
print(GT[order(contrast, celltype, class, -lfc)])

sc("dosage_9p24.csv", R)
sc("dosage_9p24_genes.csv", GT)
writeLines(format(Sys.time()), file.path(OUT, "dosage.done"))
cat("\nOK\n")
