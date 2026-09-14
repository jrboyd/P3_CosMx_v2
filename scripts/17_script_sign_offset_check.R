#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## 17 -- does the sign-consistency result survive the sample offset?
##
## THE PROBLEM THIS TESTS. `sign_consistency.csv` asks, per gene, whether IB1 exceeds RB1 in
## the same direction across all five T compartments, against a null of random direction
## (p = 2 x 0.5^5 = 0.0625 two-sided). 14 of 22 genes were consistent against 1.38 expected,
## p = 2.7e-12. That test is sound as stated -- but what it detects is a REPRODUCIBLE
## IB1-vs-RB1 direction, and a uniform sample-level offset (depth, staining, section)
## produces exactly that, in all five compartments at once, for free. So the test cannot
## distinguish "FL biology differs between the nodes" from "the two slides differ".
##
## The binomial model gives the discriminator the sign test lacks: the same genes measured
## between the two samples' REFERENCE zones (IB1 background vs RB1 interfollicular), which
## share no disease-region membership. If a gene points the same way there as it does in the
## compartment-wise comparison, its sign consistency is the offset, not the FL.
##
## Reads CSVs only.
## ---------------------------------------------------------------------------------------
suppressMessages(library(data.table))
options(width = 200)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
OUT  <- file.path("output/17_pd1", TIER)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

S <- fread(file.path(OUT, "sign_consistency.csv"))[!is.na(gene)]
B <- fread(file.path(OUT, "binomial_model.csv"))

## =======================================================================================
## CHECK 1, and the decisive one: are the 5 T compartments independent observations?
##
## The sign test's null is p = 2 x 0.5^5 = 0.0625 per gene, which assumes each compartment
## is an independent coin flip. If instead the 5 compartments mostly report the SAME
## per-gene quantity (that gene's sample-level offset) plus small noise, then "consistent
## across all 5" is ONE observation, not five, and the correct per-gene null is ~0.5.
##
## This is the project's own pseudoreplication rule (region is the unit, 154 neighbourhoods
## are not 400k observations -- step 11 lost a Tfh result to it) applied to compartments
## instead of cells. The rule was written down and then broken in a new place.
## =======================================================================================
P <- fread(file.path(OUT, "sign_consistency_per_compartment.csv"))
G <- P[, .(mu = mean(log2fc), sd = sd(log2fc), n_comp = .N,
           all_same = abs(sum(sign(log2fc))) == .N), by = gene]
G[, ratio := abs(mu) / sd]

v_between <- var(G$mu)          # variance of the per-gene means: the between-gene signal
v_within  <- mean(G$sd^2)       # mean variance ACROSS compartments within a gene
icc  <- v_between / (v_between + v_within)
neff <- G[, .N/.N * 5][1] / (1 + (5 - 1) * icc)   # Kish-style effective n per gene

cat("\n============================================================\n")
cat("CHECK 1: are the 5 T compartments 5 observations, or 1?\n")
cat("============================================================\n")
print(G[order(-ratio)])
cat(sprintf("\nvariance of per-gene means (between genes) : %.4f\n", v_between))
cat(sprintf("mean within-gene variance across compartments: %.4f\n", v_within))
cat(sprintf("intraclass correlation                       : %.3f\n", icc))
cat(sprintf("effective independent observations per gene   : %.2f  (of 5)\n", neff))
cat(sprintf("\ngenes where all 5 agree: %d; of those, |mean|/sd > 1 in %d\n",
            G[all_same == TRUE, .N], G[all_same == TRUE & ratio > 1, .N]))
cat("  -> agreement across compartments is a deterministic consequence of the gene mean\n")
cat("     being large relative to the compartment spread. It is not 5 independent draws.\n")

k <- G[all_same == TRUE, .N]; nn <- nrow(G)
cat("\n=== the global test under each null\n")
for (p0 in c(0.0625, 0.5)) {
  bt <- binom.test(k, nn, p0, alternative = "greater")
  cat(sprintf("  per-gene null p = %.4f : %d/%d, p = %.3g   %s\n", p0, k, nn, bt$p.value,
              if (bt$p.value < .05) "SIGNIFICANT" else "NOT SIGNIFICANT"))
}
cat("\n  The first null (5 independent compartments) is the one originally used and is wrong.\n")
cat("  The second (compartments are pseudo-replicates of one per-gene number, so the only\n")
cat("  random quantity is the SIGN of that number) is the defensible one. The result does\n")
cat("  not survive it. This retraction needs no appeal to the slide offset at all.\n")
fwrite(data.table(icc = round(icc, 3), n_eff_per_gene = round(neff, 2),
                  consistent = k, tested = nn,
                  p_null_5_independent = binom.test(k, nn, 0.0625, alternative = "greater")$p.value,
                  p_null_1_effective  = binom.test(k, nn, 0.5, alternative = "greater")$p.value),
       file.path(OUT, "sign_pseudoreplication.csv"))
cat("wrote", file.path(OUT, "sign_pseudoreplication.csv"), "\n")

## The reference-zone contrast, T stratum: the cleanest available estimate of the per-sample
## offset for these genes among T cells, with depth stratified out.
REF <- B[contrast == "background vs interfollicular" & stratum == "T",
         .(gene, or_ref = or, lo_ref = lo, hi_ref = hi, q_ref = q, n_ref = n_pos)]

D <- merge(S[, .(gene, direction, n_IB1_higher, med_log2, q_sign)], REF, by = "gene")
D[, sign_consistent := direction %in% c("IB1 higher in ALL", "RB1 higher in ALL")]
D[, dir_sign := fifelse(direction == "IB1 higher in ALL", 1L,
                 fifelse(direction == "RB1 higher in ALL", -1L, NA_integer_))]
D[, dir_ref := fifelse(or_ref > 1, 1L, -1L)]
D[, ref_explains := !is.na(dir_sign) & dir_sign == dir_ref]
D[, ref_is_large := q_ref < .05]

D <- D[order(-sign_consistent, -abs(log(or_ref)))]
sc("sign_offset_check.csv",
   D[, .(gene, direction, med_log2 = round(med_log2, 3), q_sign = signif(q_sign, 2),
         or_reference_zones = round(or_ref, 3),
         ci_ref = sprintf("[%.2f, %.2f]", lo_ref, hi_ref), q_ref = signif(q_ref, 2),
         reference_same_direction = ref_explains)])

cat("\n============================================================\n")
cat("CHECK 2: is sign consistency FL biology, or the sample offset?\n")
cat("============================================================\n\n")
print(D[sign_consistent == TRUE,
        .(gene, direction, med_log2 = round(med_log2, 3),
          or_reference_zones = round(or_ref, 3),
          ci_ref = sprintf("[%.2f, %.2f]", lo_ref, hi_ref), q_ref = signif(q_ref, 2),
          reference_same_direction = ref_explains)])

n_sc  <- D[sign_consistent == TRUE, .N]
n_exp <- D[sign_consistent == TRUE & ref_explains, .N]
n_big <- D[sign_consistent == TRUE & ref_explains & ref_is_large, .N]
n_esc <- D[sign_consistent == TRUE & !ref_explains, .N]

cat(sprintf("\n%d sign-consistent genes have a reference-zone estimate.\n", n_sc))
cat(sprintf("  %d of %d point the SAME WAY between the two reference zones", n_exp, n_sc))
cat(sprintf(" (%d of those at q<0.05).\n", n_big))
cat(sprintf("  %d point the OTHER way -- these are the only ones the offset cannot explain: %s\n",
            n_esc, if (n_esc) paste(D[sign_consistent & !ref_explains, gene], collapse = ", ") else "none"))

## A binomial test of the concordance itself: under "sign consistency is FL biology", the
## reference zones carry no information and should agree with the compartment direction half
## the time. Under "sign consistency is the offset", they agree almost always.
if (n_sc > 0) {
  bt <- binom.test(n_exp, n_sc, 0.5, alternative = "greater")
  cat(sprintf("\nreference/compartment direction concordance %d/%d, p = %.3g\n",
              n_exp, n_sc, bt$p.value))
  cat("  (null = sign consistency carries FL-specific information, so the reference zones\n")
  cat("   should agree at chance. Rejecting it means the offset explains the sign test.)\n")
}

cat("\n=== verdict\n")
if (n_esc == 0) {
  cat("  EVERY sign-consistent gene points the same way between the two samples' reference\n")
  cat("  zones. The sign-consistency result is a SAMPLE-LEVEL OFFSET and must not be read as\n")
  cat("  FL biology. The test was correctly executed against the wrong null.\n")
} else {
  cat(sprintf("  %d gene(s) survive: sign-consistent AND opposite in the reference zones.\n", n_esc))
  cat("  Only these are candidates for an FL-specific difference; the rest are the offset.\n")
}

writeLines(format(Sys.time()), file.path(OUT, "sign_offset_check.done"))
cat("\nOK\n")
