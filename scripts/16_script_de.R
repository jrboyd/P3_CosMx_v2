## ============================================================================
## 16_script_de.R -- the four contrasts, region as the unit.
##
## Reuses step 14's engine wholesale: pseudobulk -> cp10k AND detection rate -> vectorised
## Mann-Whitney with the region as the unit -> BH per contrast -> a concordance gate that
## only trusts a result when cp10k and detection move the same way. The two things that
## engine exists to defeat both apply here:
##
##   cp10k IS COMPOSITIONAL. It sums to 1e4, so a real rise in a few abundant genes fakes a
##   fall in everything else. Step 14 found 4 housekeeping genes "down 2x" in all 9 clusters
##   whose DETECTION rate had barely moved. Direction is therefore only believed when cp10k
##   and detection agree (14:555).
##
##   REGION IS THE UNIT. Cells inside a region are not independent observations of it.
##
## AND ONE MORE, SPECIFIC TO THIS STEP: bins inside one FOV are not independent either.
## IB1's FL_focus is 44 bins in 9 FOVs (strict) / 89 bins in 23 FOVs (matched). So every
## contrast is run TWICE -- bin as unit, and bins collapsed to FOV means -- and a gene that
## only clears at bin level is reported as unresolved, not as a finding.
##
## THE FOUR CONTRASTS
##   C1  IB1 FL_focus vs IB1 background          -- within IB1; a check against step 14
##   C2  RB1 FL_follicular vs RB1 interfollicular -- within RB1
##   C3  IB1 FL_focus vs RB1 FL_follicular        -- the direct cross-sample contrast
##   C4  the DIFFERENCE of C1 and C2 per gene     -- PRIMARY
##
## WHY C4 IS PRIMARY. C3 compares two samples whose QC, depth and damage differ; no amount
## of matching fully removes that. C4 asks instead whether FL differs from ITS OWN
## surroundings differently in the two nodes, so any per-sample technical shift that acts on
## both of a sample's regions cancels. C3 is reported as the effect size, C4 as the evidence.
##
## C4's TWO REFERENCE GROUPS ARE NOT THE SAME TISSUE, and this is stated rather than
## finessed: 16a established that IB1's B-poor tissue is myeloid-rich and T-DEPLETED (T share
## among assigned 0.53x, myeloid 1.86x) while RB1's is a genuine T-zone (T 2.07x, myeloid
## 5.09x). So C4 is computed against BOTH available references -- B_sparse and B_mid -- and a
## gene is only called if the two agree.
##
## Usage: Rscript scripts/16_script_de.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "strict"

PB   <- file.path("output/16_pseudobulk", TIER)
OUT  <- file.path("output/16_de", TIER)
for (f in c(file.path(PB,"pb_sum.Rds"), file.path(PB,"pb_det.Rds"),
            file.path(PB,"bin_meta.csv.gz")))
  if (!file.exists(f)) stop("missing ", f, " -- run 16_script_pseudobulk.R ", TIER, call.=FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

MIN_REG   <- 12L    # step 14's floor: fewer regions than this on either side is untestable
MIN_CELLS <- 40L    # a bin already cleared this in 16a; re-asserted here
DET_MIN   <- 0.30   # |log2 detection ratio| a concordant call must clear (step 14:555)
Q_CUT     <- 0.05
N_BOOT    <- 2000L
SEED      <- 1L

## Genes that DEFINE a region and are therefore circular within it. Reported, never dropped
## silently -- a reader must be able to see that IGKC is top of the FL list and know why.
SEL <- list(
  `FL foci (light chain)`      = c("IGKC","IGLC1/2"),
  `cHL hotspots (chemokine)`   = c("CCL17","CCL22"),
  `B_dense architecture`       = c("MS4A1","CD79A","CD79B","BANK1"))
SEL_ALL <- unlist(SEL, use.names = FALSE)

cat("============================================================\n")
cat("step 16d -- region-unit DE, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

S <- readRDS(file.path(PB,"pb_sum.Rds")); D <- readRDS(file.path(PB,"pb_det.Rds"))
BM <- fread(file.path(PB,"bin_meta.csv.gz"))
setkey(BM, bin)
BM <- BM[colnames(S)]
if (any(BM$bin != colnames(S))) stop("bin_meta does not align with the pseudobulk", call.=FALSE)
BM <- BM[cells >= MIN_CELLS]
S <- S[, BM$bin, drop = FALSE]; D <- D[, BM$bin, drop = FALSE]
GENES <- rownames(S)
cat(sprintf("genes %d | bins %d\n", length(GENES), ncol(S)))
print(BM[, .(bins = .N, fovs = uniqueN(fov)), by = .(samp, bin_region)][order(samp,-bins)])

## cp10k per bin (share of the bin's transcriptome) and detection rate per bin.
CP  <- as.matrix(S) %*% Diagonal(x = 1e4 / pmax(Matrix::colSums(S), 1))
DET <- as.matrix(D) %*% Diagonal(x = 1 / pmax(BM$cells, 1))
dimnames(CP) <- dimnames(DET) <- dimnames(S)

## Vectorised Mann-Whitney with the normal approximation, per gene across regions --
## 14_script_lightchain_in_out.R:481-489.
mw <- function(A, B) {
  na <- ncol(A); nb <- ncol(B)
  R  <- cbind(A, B)
  rk <- t(apply(R, 1, rank))
  Ra <- rowSums(rk[, seq_len(na), drop = FALSE])
  U  <- Ra - na*(na+1)/2
  mu <- na*nb/2; sg <- sqrt(na*nb*(na+nb+1)/12)
  z  <- (U - mu)/sg
  list(z = z, p = 2*pnorm(-abs(z)))
}

## One contrast, both units. `collapse` averages bins within an FOV first, which is the
## conservative reading -- it treats each FOV as one observation.
run1 <- function(name, ia, ib, lab_a, lab_b) {
  if (length(ia) < MIN_REG || length(ib) < MIN_REG) {
    cat(sprintf("\n%s: SKIPPED -- %d vs %d regions, below MIN_REG = %d\n",
                name, length(ia), length(ib), MIN_REG))
    return(NULL)
  }
  cp_a <- CP[, ia, drop=FALSE]; cp_b <- CP[, ib, drop=FALSE]
  dt_a <- DET[, ia, drop=FALSE]; dt_b <- DET[, ib, drop=FALSE]
  st <- mw(cp_a, cp_b)
  out <- data.table(contrast = name, group_a = lab_a, group_b = lab_b,
                    gene = GENES, n_a = length(ia), n_b = length(ib),
                    cp10k_a = rowMeans(cp_a), cp10k_b = rowMeans(cp_b),
                    det_a = rowMeans(dt_a), det_b = rowMeans(dt_b),
                    z = st$z, p = st$p)
  out[, q := p.adjust(p, "BH")]
  out[, log2fc     := log2((cp10k_a + .1)/(cp10k_b + .1))]
  out[, log2fc_det := log2((det_a + 1e-4)/(det_b + 1e-4))]
  ## FOV-collapsed rerun
  fa <- BM$fov[ia]; fb <- BM$fov[ib]
  ca <- do.call(cbind, lapply(split(seq_along(ia), fa), function(k) rowMeans(cp_a[,k,drop=FALSE])))
  cb <- do.call(cbind, lapply(split(seq_along(ib), fb), function(k) rowMeans(cp_b[,k,drop=FALSE])))
  if (ncol(ca) >= MIN_REG && ncol(cb) >= MIN_REG) {
    s2 <- mw(ca, cb)
    out[, `:=`(p_fov = s2$p, n_fov_a = ncol(ca), n_fov_b = ncol(cb))]
    out[, q_fov := p.adjust(p_fov, "BH")]
    out[, log2fc_fov := log2((rowMeans(ca) + .1)/(rowMeans(cb) + .1))]
  } else {
    out[, `:=`(p_fov = NA_real_, q_fov = NA_real_, log2fc_fov = NA_real_,
               n_fov_a = ncol(ca), n_fov_b = ncol(cb))]
    cat(sprintf("  %s: FOV-collapsed test not run (%d vs %d FOVs, need %d each)\n",
                name, ncol(ca), ncol(cb), MIN_REG))
  }
  out[, concordant := sign(log2fc) == sign(log2fc_det) & abs(log2fc_det) > DET_MIN]
  out[, is_selection := gene %in% SEL_ALL]
  out[, verdict := fifelse(!concordant, "cp10k/detection discordant",
                   fifelse(q >= Q_CUT, "not significant",
                   fifelse(is.na(q_fov), "bin-level only (no FOV test)",
                   fifelse(q_fov < Q_CUT, "significant, survives FOV collapse",
                           "bin-level only -- UNRESOLVED"))))]
  cat(sprintf("\n=== %s: %s (n=%d) vs %s (n=%d)%s\n", name, lab_a, length(ia), lab_b,
              length(ib), if (!is.na(out$n_fov_a[1]))
                sprintf(" | FOVs %d vs %d", out$n_fov_a[1], out$n_fov_b[1]) else ""))
  print(out[, .N, by = verdict][order(-N)])
  print(head(out[concordant & q < Q_CUT & !is_selection][order(-abs(log2fc)),
        .(gene, cp10k_a = round(cp10k_a,2), cp10k_b = round(cp10k_b,2),
          log2fc = round(log2fc,2), det_log2 = round(log2fc_det,2),
          q = signif(q,2), q_fov = signif(q_fov,2))], 12))
  if (out[is_selection & q < Q_CUT, .N])
    cat("  selection genes (circular, shown for audit): ",
        paste(out[is_selection & q < Q_CUT, sprintf("%s %+.2f", gene, log2fc)], collapse="; "),
        "\n", sep="")
  out[]
}

ix <- function(...) { a <- list(...); which(Reduce(`&`, a)) }
IB <- BM$samp == "IB1"; RB <- BM$samp == "RB1"
R  <- BM$bin_region

C1 <- run1("C1", ix(IB, R=="FL_focus"),      ix(IB, R=="background"),
           "IB1 FL_focus", "IB1 background")
C2 <- run1("C2", ix(RB, R=="FL_follicular"), ix(RB, R=="interfollicular"),
           "RB1 FL_follicular", "RB1 interfollicular")
C2b<- run1("C2b", ix(RB, R=="FL_follicular"), ix(RB, R=="follicular_mid"),
           "RB1 FL_follicular", "RB1 follicular_mid")
C3 <- run1("C3", ix(IB, R=="FL_focus"),      ix(RB, R=="FL_follicular"),
           "IB1 FL_focus", "RB1 FL_follicular")
C5 <- run1("C5", ix(IB, R=="cHL_hotspot"),   ix(IB, R=="background"),
           "IB1 cHL_hotspot", "IB1 background")
ALL <- rbindlist(list(C1,C2,C2b,C3,C5), fill = TRUE)
fwrite(ALL, file.path(OUT, "de_contrasts.csv.gz"), compress = "gzip")
cat("\nwrote ", file.path(OUT, "de_contrasts.csv.gz"), " (", nrow(ALL), " rows)\n", sep="")

## ---- C4: difference in differences ------------------------------------------
## Per gene, (FL vs its own surroundings in IB1) minus (the same in RB1). A bootstrap over
## REGIONS gives the interval; the point estimate alone would hide how few FL bins there are.
did <- function(name, a1, b1, a2, b2) {
  if (min(length(a1), length(b1), length(a2), length(b2)) < MIN_REG) {
    cat(sprintf("\n%s: SKIPPED -- region counts %d/%d/%d/%d below MIN_REG\n",
                name, length(a1), length(b1), length(a2), length(b2))); return(NULL) }
  lf <- function(A, B) log2((rowMeans(CP[,A,drop=FALSE]) + .1)/(rowMeans(CP[,B,drop=FALSE]) + .1))
  ld <- function(A, B) log2((rowMeans(DET[,A,drop=FALSE]) + 1e-4)/(rowMeans(DET[,B,drop=FALSE]) + 1e-4))
  obs   <- lf(a1,b1) - lf(a2,b2)
  obs_d <- ld(a1,b1) - ld(a2,b2)
  set.seed(SEED)
  bs <- vapply(seq_len(N_BOOT), function(i)
    lf(sample(a1, length(a1), TRUE), sample(b1, length(b1), TRUE)) -
    lf(sample(a2, length(a2), TRUE), sample(b2, length(b2), TRUE)), numeric(length(GENES)))
  lo <- apply(bs, 1, quantile, .025); hi <- apply(bs, 1, quantile, .975)
  ## Bootstrap p: the share of resamples on the opposite side of zero, doubled.
  pb <- 2*pmin(rowMeans(bs <= 0), rowMeans(bs >= 0))
  d <- data.table(contrast = name, gene = GENES,
                  did = obs, did_det = obs_d, ci_lo = lo, ci_hi = hi,
                  p_boot = pmax(pb, 1/N_BOOT))
  d[, q_boot := p.adjust(p_boot, "BH")]
  d[, excludes_zero := (ci_lo > 0) | (ci_hi < 0)]
  d[, concordant := sign(did) == sign(did_det)]
  d[, is_selection := gene %in% SEL_ALL]
  cat(sprintf("\n=== %s: DiD, %d genes | excluding zero %d | concordant & q<%.2f: %d\n",
              name, nrow(d), d[, sum(excludes_zero)], Q_CUT,
              d[excludes_zero & concordant & q_boot < Q_CUT, .N]))
  print(head(d[excludes_zero & concordant & q_boot < Q_CUT & !is_selection][order(-abs(did)),
        .(gene, did = round(did,2), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi),
          did_det = round(did_det,2), q_boot = signif(q_boot,2))], 15))
  d[]
}
D1 <- did("C4_vs_sparse", ix(IB,R=="FL_focus"), ix(IB,R=="background"),
                          ix(RB,R=="FL_follicular"), ix(RB,R=="interfollicular"))
D2 <- did("C4_vs_mid",    ix(IB,R=="FL_focus"), ix(IB,R=="background"),
                          ix(RB,R=="FL_follicular"), ix(RB,R=="follicular_mid"))
DD <- rbindlist(list(D1,D2), fill = TRUE)
if (nrow(DD)) {
  fwrite(DD, file.path(OUT, "de_did.csv.gz"), compress = "gzip")
  cat("\nwrote ", file.path(OUT, "de_did.csv.gz"), "\n", sep="")
  ## A gene is only called if BOTH reference groups agree -- see the header on why the two
  ## references are not the same tissue.
  w <- dcast(DD[, .(contrast, gene, did, q_boot, excludes_zero, concordant, is_selection)],
             gene + is_selection ~ contrast, value.var = c("did","q_boot","excludes_zero","concordant"))
  if (all(c("did_C4_vs_sparse","did_C4_vs_mid") %in% names(w))) {
    w[, both := excludes_zero_C4_vs_sparse & excludes_zero_C4_vs_mid &
                concordant_C4_vs_sparse & concordant_C4_vs_mid &
                q_boot_C4_vs_sparse < Q_CUT & q_boot_C4_vs_mid < Q_CUT &
                sign(did_C4_vs_sparse) == sign(did_C4_vs_mid)]
    cat(sprintf("\n=== C4 genes surviving BOTH reference groups: %d\n", w[, sum(both)]))
    print(head(w[both == TRUE & is_selection == FALSE][order(-abs(did_C4_vs_sparse)),
          .(gene, did_sparse = round(did_C4_vs_sparse,2), did_mid = round(did_C4_vs_mid,2))], 20))
    sc("de_did_consensus.csv", w[both == TRUE][order(-abs(did_C4_vs_sparse))])
  }
}
writeLines(format(Sys.time()), file.path(OUT, "de.done"))
cat("============================================================\nOK\n")
