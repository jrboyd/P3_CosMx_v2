#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## 17 -- difference-in-differences on the log-odds scale, WITHIN cell type.
##
## Why this exists, and why it is not step 16's withdrawn C4.
##
## Step 16's C4 differenced two bin-level pseudobulk contrasts and was withdrawn because the
## two reference groups are different tissue: IB1's B-poor tissue is myeloid-rich and
## T-DEPLETED (0.53x T / 1.86x myeloid) while RB1's is a real T-zone (2.07x T / 5.09x
## myeloid). Subtracting a T-zone-referenced contrast from a myeloid-referenced one leaves
## the T-zone composition difference behind, dressed as a result.
##
## Conditioning on cell type removes exactly that objection: a contrast computed among T
## cells only cannot be moved by how many T cells there are. What remains is a per-cell
## positivity comparison with sequencing depth stratified out by decile.
##
## The assumption this still rests on, which is UNTESTABLE here (n = 1 patient):
##   the sample-level offset in log-odds is the SAME in the reference region as in the
##   focus region -- i.e. no sample x region interaction other than the biology of interest.
## A DiD is only as good as that. So this is reported as a DEFLATION analysis -- how much of
## a raw cross-sample odds ratio survives once the reference offset is removed -- and never
## as a primary result. Where the reference contrast is itself large (LAG3, CXCL13) the raw
## cross-sample number is uninterpretable and the deflated one is the only honest version.
##
## The one place the assumption is checkable is the checkpoint axis itself: the reference
## contrast is null for PDCD1/CD274/PDCD1LG2 in the T stratum, which bounds the offset.
##
## Reads only binomial_model.csv. No refit; runs in seconds off CSVs.
## ---------------------------------------------------------------------------------------
suppressMessages(library(data.table))
options(width = 210)

TIER <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "matched"
OUT  <- file.path("output/17_pd1", TIER)
IN   <- file.path(OUT, "binomial_model.csv")
stopifnot(file.exists(IN))
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }

R <- fread(IN)

## binomial_model.csv now carries TWO reference definitions (all_background /
## far_background). Merging without accounting for that would cross-join them and silently
## pair an IB1-background reference with a far-background focus. The FOCUS contrast
## (FL_focus vs FL_follicular) touches no IB1 background so it is identical in both; only the
## REFERENCE contrast differs. Assert that, then take the reference from the CLEAN definition,
## because field_extent_verdict.csv shows 40% of IB1 background carries Hodgkin field signal.
if ("refdef" %in% names(R)) {
  chk <- dcast(R[contrast == "FL_focus vs FL_follicular" & !is.na(or)],
               stratum + gene ~ refdef, value.var = "or")
  if (all(c("all_background","far_background") %in% names(chk))) {
    dmax <- chk[, max(abs(all_background - far_background), na.rm = TRUE)]
    cat(sprintf("focus contrast identical across reference definitions? max |diff| = %.2e\n", dmax))
    stopifnot(dmax < 1e-8)
  }
  REFDEF_USE <- "far_background"
  cat("using reference definition:", REFDEF_USE, "(field-contaminated background excluded)\n")
} else REFDEF_USE <- NA_character_

## Percentile bootstrap CIs are asymmetric and a few rows hit the boundary (lo == 0), where a
## log-scale SE is undefined. Those rows are dropped rather than fudged -- a DiD built on a
## boundary CI is not a number.
R[, `:=`(b = log(or), se = (log(hi) - log(lo)) / (2 * qnorm(.975)))]
R[!is.finite(se) | se <= 0 | !is.finite(b), `:=`(b = NA_real_, se = NA_real_)]

FOCUS <- "FL_focus vs FL_follicular"        # IB1 FL         vs RB1 FL
REF   <- "background vs interfollicular"    # IB1 reference  vs RB1 reference

RS <- if (!is.na(REFDEF_USE)) R[refdef == REFDEF_USE] else R
f <- RS[contrast == FOCUS, .(stratum, gene, b_f = b, se_f = se, or_f = or, q_f = q, n_f = n_pos)]
r <- RS[contrast == REF,   .(stratum, gene, b_r = b, se_r = se, or_r = or, q_r = q, n_r = n_pos)]
stopifnot(!anyDuplicated(f[, .(stratum, gene)]), !anyDuplicated(r[, .(stratum, gene)]))
D <- merge(f, r, by = c("stratum", "gene"))

## The two contrasts use DISJOINT region sets (IB1 FL + RB1 follicular vs IB1 background +
## RB1 interfollicular) and were bootstrapped over FOVs independently, so treating the two
## log-odds estimates as independent is reasonable and the DiD variance is additive.
D[, `:=`(b_d = b_f - b_r, se_d = sqrt(se_f^2 + se_r^2))]
D[, `:=`(or_did = exp(b_d),
         lo_did = exp(b_d - qnorm(.975) * se_d),
         hi_did = exp(b_d + qnorm(.975) * se_d),
         p_did  = 2 * pnorm(-abs(b_d / se_d)))]
D[, q_did := p.adjust(p_did, "BH"), by = stratum]

## WITHOUT THIS THE NULL IS UNINTERPRETABLE. The DiD variance is the SUM of two bootstrap
## variances, so the DiD is necessarily less powerful than either contrast that built it.
## "nothing survived" is only a negative result if something could have been seen; report the
## minimum detectable effect at 80% power beside every null.
D[, mde_fold := exp((qnorm(.975) + qnorm(.80)) * se_d)]

## And separate the two reasons a DiD can come back null, because they mean opposite things:
##   reference contrast LARGE -> the raw cross-sample effect was the slide. Retracted.
##   reference contrast NULL  -> no offset to remove, so the DiD ~ the raw effect and the
##                               null is a POWER failure. Unresolved, not refuted.
D[, ref_class := fifelse(q_r < .05, "reference LARGE: offset",
                                    "reference NULL: power-limited")]

## How much of the raw cross-sample effect was the reference offset. Reported as the share of
## the log-odds removed, so 100% means the whole apparent difference was the sample offset.
D[, deflated_pct := round(100 * (1 - b_d / b_f), 1)]
D[!is.finite(deflated_pct) | abs(b_f) < 1e-8, deflated_pct := NA_real_]

D <- D[order(stratum, q_did)]
sc("did_within_celltype.csv",
   D[, .(stratum, gene, n_focus = n_f, n_ref = n_r,
         or_focus_raw = round(or_f, 3), q_focus_raw = signif(q_f, 2),
         or_reference = round(or_r, 3), q_reference = signif(q_r, 2),
         or_did = round(or_did, 3),
         ci_did = sprintf("[%.2f, %.2f]", lo_did, hi_did),
         p_did = signif(p_did, 3), q_did = signif(q_did, 2),
         mde_fold = round(mde_fold, 2), ref_class, deflated_pct)])

cat("\n============================================================\n")
cat("DiD within cell type:", FOCUS, "MINUS", REF, "\n")
cat("============================================================\n")

cat("\n=== is the DiD licensed? the reference contrast on the PD-1 axis, T stratum\n")
cat("    (a null reference contrast bounds the sample offset; a large one does not)\n")
print(RS[contrast == REF & stratum == "T" & gene %in% c("PDCD1", "CD274", "PDCD1LG2"),
        .(gene, n_pos, or_reference = round(or, 3),
          ci = sprintf("[%.2f, %.2f]", lo, hi), q = signif(q, 2))])

cat("\n=== reference contrast, everything NOT null -- these genes' raw cross-sample ORs are\n")
cat("    uninterpretable without the deflation\n")
print(RS[contrast == REF & !is.na(q) & q < .05 & stratum != "all",
        .(stratum, gene, or_reference = round(or, 3), q = signif(q, 2))][order(stratum, -abs(log(or_reference)))])

for (st in c("T", "B", "myeloid")) {
  cat(sprintf("\n=== stratum '%s'\n", st))
  z <- D[stratum == st & is.finite(or_did)]
  if (!nrow(z)) { cat("  nothing estimable\n"); next }
  print(z[, .(gene, or_focus_raw = round(or_f, 3), or_reference = round(or_r, 3),
              or_did = round(or_did, 3), ci = sprintf("[%.2f, %.2f]", lo_did, hi_did),
              q_did = signif(q_did, 2), deflated_pct)])
  s <- z[q_did < .05]
  cat(sprintf("  -> %d of %d survive the DiD at q<0.05: %s\n", nrow(s), nrow(z),
              if (nrow(s)) paste(s$gene, collapse = ", ") else "none"))
}

cat("\n=== what the null can and cannot exclude\n")
print(D[stratum != "all", .(genes = .N, median_mde_fold = round(median(mde_fold), 2),
        largest_or_did = round(max(or_did), 2), n_p_under_05 = sum(p_did < .05)),
        by = .(stratum, ref_class)][order(stratum, ref_class)])
cat("\n  Point estimates >= 1.5 that the DiD could NOT exclude -- these are unresolved,\n")
cat("  not refuted, and the ones with a NULL reference contrast have no offset to blame:\n")
print(D[stratum != "all" & or_did >= 1.5,
        .(stratum, gene, or_raw = round(or_f, 2), or_reference = round(or_r, 2),
          or_did = round(or_did, 2), ci = sprintf("[%.2f, %.2f]", lo_did, hi_did),
          p_did = signif(p_did, 3), mde_fold = round(mde_fold, 2), ref_class)][order(p_did)])

cat("\n=== the PD-1 axis specifically, all strata\n")
print(D[gene %in% c("PDCD1", "CD274", "PDCD1LG2") & is.finite(or_did),
        .(stratum, gene, or_focus_raw = round(or_f, 3), or_reference = round(or_r, 3),
          or_did = round(or_did, 3), ci = sprintf("[%.2f, %.2f]", lo_did, hi_did),
          q_did = signif(q_did, 2))][order(gene, stratum)])

writeLines(format(Sys.time()), file.path(OUT, "did_logor.done"))
cat("\nOK\n")
