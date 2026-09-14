## ============================================================================
## 16_script_celltype_de.R -- DE conditioned on cell type. The composition-free contrast.
##
## WHY THIS EXISTS. 16c established that the cross-sample tissue comparison is not
## identifiable: IB1's FL foci are 17.6% B cells against RB1's follicles at 50.9%, and
## matching bins on composition matched 2 of 44. The unmatched contrast (16d C3) then
## returned RB1's interfollicular markers as "up in IB1 FL", and the difference-in-differences
## returned the same list because the two reference groups are also different tissue.
##
## Conditioning on cell type removes the confound by construction rather than by matching:
## a B cell is a B cell in either node, so IB1 FL-focus B cells against RB1 follicular B cells
## has no composition term left. This is the move step 14 made when it tested expression
## IN vs OUT *within each cluster* instead of pooled.
##
## WHAT IS STILL CONFOUNDED, AND WHAT IS DONE ABOUT IT: sequencing depth. IB1's FL bins run
## ~566 median counts against RB1's follicular bins at ~722 (matched tier), and depth predicts
## sparse-gene detection at AUC ~0.70 on this platform for EVERY gene, positive controls
## included (qc_16_negprobe_floor.R). cp10k removes a uniform depth shift but detection rate
## does not, and the concordance gate needs both. So every cross-sample contrast here is run
## on depth- and size-MATCHED groups, with the achieved balance reported. An unmatched run is
## reported alongside so the effect of matching is visible rather than asserted.
##
## COARSE TYPES ONLY GET P-VALUES. B/T/myeloid clear the 5-cell floor in enough bins to make
## the bin a usable unit (IB1 FL_focus: 80 B / 69 T / 53 myeloid groups). The fine subtypes
## -- Tfh, Treg, CD8_T, Macrophage, cDC -- are 0.03-1.5% of cells, i.e. under one cell per
## bin, and are reported as pooled effect sizes with no test. That limit is printed.
##
## Usage: Rscript scripts/16_script_celltype_de.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "matched"

PB  <- file.path("output/16_ct_pseudobulk", TIER)
OUT <- file.path("output/16_ct_de", TIER)
for (f in c("pb_binct_sum.Rds","pb_binct_det.Rds","binct_meta.csv.gz",
            "pb_regct_sum.Rds","pb_regct_det.Rds","regct_meta.csv.gz"))
  if (!file.exists(file.path(PB,f))) stop("missing ", file.path(PB,f), call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

MIN_REG <- 12L; DET_MIN <- 0.30; Q_CUT <- 0.05; CALIPER <- 0.25; SEED <- 1L
SEL <- c("IGKC","IGLC1/2","CCL17","CCL22","MS4A1","CD79A","CD79B","BANK1")

cat("============================================================\n")
cat("step 16m -- cell-type-conditioned DE, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

S <- readRDS(file.path(PB,"pb_binct_sum.Rds")); D <- readRDS(file.path(PB,"pb_binct_det.Rds"))
GB <- fread(file.path(PB,"binct_meta.csv.gz")); setkey(GB, grp_bin); GB <- GB[colnames(S)]
stopifnot(identical(GB$grp_bin, colnames(S)))
GENES <- rownames(S)
GB[, `:=`(log_depth = log(nCount/cells), log_cells = log(cells))]
CP  <- as.matrix(S) %*% Diagonal(x = 1e4/pmax(Matrix::colSums(S),1))
DET <- as.matrix(D) %*% Diagonal(x = 1/pmax(GB$cells,1))
dimnames(CP) <- dimnames(DET) <- dimnames(S)
cat(sprintf("genes %d | (bin x celltype) groups %d\n", length(GENES), ncol(S)))

mw <- function(A,B) {
  na <- ncol(A); nb <- ncol(B); rk <- t(apply(cbind(A,B), 1, rank))
  U <- rowSums(rk[, seq_len(na), drop=FALSE]) - na*(na+1)/2
  z <- (U - na*nb/2)/sqrt(na*nb*(na+nb+1)/12)
  list(z=z, p=2*pnorm(-abs(z)))
}

## Greedy 1:1 match on standardised (log depth, log cells), same construction as 16c.
match_idx <- function(ia, ib) {
  V <- c("log_depth","log_cells")
  M <- rbind(as.matrix(GB[ia, ..V]), as.matrix(GB[ib, ..V]))
  mu <- colMeans(M); sd <- apply(M, 2, stats::sd); sd[sd==0] <- 1
  za <- scale(as.matrix(GB[ia, ..V]), mu, sd); zb <- scale(as.matrix(GB[ib, ..V]), mu, sd)
  set.seed(SEED); used <- rep(FALSE, length(ib)); pk <- rep(NA_integer_, length(ia))
  for (i in sample.int(length(ia))) {
    d <- sqrt(colSums((t(zb) - za[i,])^2)); d[used] <- Inf
    j <- which.min(d); if (is.finite(d[j]) && d[j] <= CALIPER) { pk[i] <- j; used[j] <- TRUE }
  }
  ok <- !is.na(pk)
  list(a = ia[ok], b = ib[pk[ok]], n_a = length(ia), matched = sum(ok),
       smd_depth_before = (mean(GB$log_depth[ia]) - mean(GB$log_depth[ib]))/sd[["log_depth"]],
       smd_depth_after  = if (!sum(ok)) NA_real_ else
         (mean(GB$log_depth[ia[ok]]) - mean(GB$log_depth[ib[pk[ok]]]))/sd[["log_depth"]])
}

run1 <- function(name, ct, ia, ib, la, lb, do_match) {
  if (length(ia) < MIN_REG || length(ib) < MIN_REG) {
    cat(sprintf("  %-28s %-8s SKIPPED (%d vs %d groups, need %d)\n", name, ct,
                length(ia), length(ib), MIN_REG)); return(NULL) }
  mi <- NULL
  if (do_match) {
    mi <- match_idx(ia, ib)
    if (mi$matched < MIN_REG) {
      cat(sprintf("  %-28s %-8s SKIPPED after matching (%d of %d within caliper)\n",
                  name, ct, mi$matched, mi$n_a)); return(NULL) }
    ia <- mi$a; ib <- mi$b
  }
  st <- mw(CP[,ia,drop=FALSE], CP[,ib,drop=FALSE])
  o <- data.table(contrast=name, celltype=ct, gene=GENES, group_a=la, group_b=lb,
                  n_a=length(ia), n_b=length(ib), matched=do_match,
                  cp10k_a=rowMeans(CP[,ia,drop=FALSE]), cp10k_b=rowMeans(CP[,ib,drop=FALSE]),
                  det_a=rowMeans(DET[,ia,drop=FALSE]), det_b=rowMeans(DET[,ib,drop=FALSE]),
                  z=st$z, p=st$p)
  o[, q := p.adjust(p,"BH")]
  o[, log2fc := log2((cp10k_a+.1)/(cp10k_b+.1))]
  o[, log2fc_det := log2((det_a+1e-4)/(det_b+1e-4))]
  fa <- split(seq_along(ia), GB$fov[ia]); fb <- split(seq_along(ib), GB$fov[ib])
  if (length(fa) >= MIN_REG && length(fb) >= MIN_REG) {
    ca <- do.call(cbind, lapply(fa, function(k) rowMeans(CP[,ia[k],drop=FALSE])))
    cb <- do.call(cbind, lapply(fb, function(k) rowMeans(CP[,ib[k],drop=FALSE])))
    s2 <- mw(ca, cb)
    o[, `:=`(p_fov=s2$p, n_fov_a=ncol(ca), n_fov_b=ncol(cb))][, q_fov := p.adjust(p_fov,"BH")]
  } else o[, `:=`(p_fov=NA_real_, q_fov=NA_real_, n_fov_a=length(fa), n_fov_b=length(fb))]
  o[, concordant := sign(log2fc)==sign(log2fc_det) & abs(log2fc_det) > DET_MIN]
  o[, is_selection := gene %in% SEL]
  o[, `:=`(smd_depth_before = if (is.null(mi)) NA_real_ else mi$smd_depth_before,
           smd_depth_after  = if (is.null(mi)) NA_real_ else mi$smd_depth_after)]
  nsig <- o[concordant & q < Q_CUT & !is_selection, .N]
  nfov <- o[concordant & q < Q_CUT & !is_selection & !is.na(q_fov) & q_fov < Q_CUT, .N]
  cat(sprintf("  %-28s %-8s n %4d vs %4d | FOV %3d vs %3d | sig %4d, FOV-robust %4d%s\n",
              name, ct, length(ia), length(ib), o$n_fov_a[1], o$n_fov_b[1], nsig, nfov,
              if (!is.null(mi)) sprintf(" | depth SMD %.2f -> %.2f", mi$smd_depth_before,
                                        mi$smd_depth_after) else ""))
  o[]
}

CTS <- c("B","T","myeloid")
IB <- GB$samp=="IB1"; RB <- GB$samp=="RB1"; R <- GB$bin_region; CC <- GB$coarse
ixx <- function(...) which(Reduce(`&`, list(...)))

ALL <- rbindlist(lapply(CTS, function(ct) rbindlist(list(
  run1("CT1 IB1 FL vs IB1 bg",     ct, ixx(IB,CC==ct,R=="FL_focus"),   ixx(IB,CC==ct,R=="background"),
       "IB1 FL_focus","IB1 background", TRUE),
  run1("CT2 RB1 foll vs RB1 inter",ct, ixx(RB,CC==ct,R=="FL_follicular"), ixx(RB,CC==ct,R=="interfollicular"),
       "RB1 FL_follicular","RB1 interfollicular", TRUE),
  run1("CT3 IB1 FL vs RB1 foll",   ct, ixx(IB,CC==ct,R=="FL_focus"),   ixx(RB,CC==ct,R=="FL_follicular"),
       "IB1 FL_focus","RB1 FL_follicular", TRUE),
  run1("CT3u UNMATCHED",           ct, ixx(IB,CC==ct,R=="FL_focus"),   ixx(RB,CC==ct,R=="FL_follicular"),
       "IB1 FL_focus","RB1 FL_follicular", FALSE),
  run1("CT4 IB1 hot vs IB1 bg",    ct, ixx(IB,CC==ct,R=="cHL_hotspot"),ixx(IB,CC==ct,R=="background"),
       "IB1 cHL_hotspot","IB1 background", TRUE)), fill=TRUE)), fill=TRUE)

if (!nrow(ALL)) stop("no contrast was testable", call. = FALSE)
fwrite(ALL, file.path(OUT,"ct_de.csv.gz"), compress="gzip")
cat("\nwrote ", file.path(OUT,"ct_de.csv.gz"), " (", nrow(ALL), " rows)\n", sep="")

for (cn in unique(ALL$contrast)) for (ct in CTS) {
  o <- ALL[contrast==cn & celltype==ct]; if (!nrow(o)) next
  h <- o[concordant & q<Q_CUT & !is_selection][order(-abs(log2fc))]
  if (!nrow(h)) next
  cat(sprintf("\n=== %s | %s cells -- top by |log2fc|\n", cn, ct))
  print(head(h[, .(gene, cp10k_a=round(cp10k_a,2), cp10k_b=round(cp10k_b,2),
                   log2fc=round(log2fc,2), det_log2=round(log2fc_det,2),
                   q=signif(q,2), q_fov=signif(q_fov,2))], 10))
}

## Does matching change the answer? If CT3 and CT3u disagree, depth was doing the work.
m <- merge(ALL[contrast=="CT3 IB1 FL vs RB1 foll", .(celltype,gene,lf_m=log2fc,sig_m=concordant&q<Q_CUT)],
           ALL[contrast=="CT3u UNMATCHED",         .(celltype,gene,lf_u=log2fc,sig_u=concordant&q<Q_CUT)],
           by=c("celltype","gene"))
if (nrow(m)) {
  cat("\n=== matched vs unmatched CT3: did depth matching change the call?\n")
  print(m[, .(genes=.N, sig_matched=sum(sig_m), sig_unmatched=sum(sig_u),
              both=sum(sig_m&sig_u), lost_by_matching=sum(sig_u&!sig_m),
              gained=sum(sig_m&!sig_u), cor_log2fc=round(cor(lf_m,lf_u),3)), by=celltype])
  sc("ct3_matching_effect.csv", m)
}

## ---- fine subtypes: pooled effect sizes, no test ----------------------------
RS <- readRDS(file.path(PB,"pb_regct_sum.Rds")); RD <- readRDS(file.path(PB,"pb_regct_det.Rds"))
GR <- fread(file.path(PB,"regct_meta.csv.gz")); setkey(GR, grp_reg); GR <- GR[colnames(RS)]
rcp <- as.matrix(RS) %*% Diagonal(x=1e4/pmax(Matrix::colSums(RS),1))
rdt <- as.matrix(RD) %*% Diagonal(x=1/pmax(GR$cells,1))
dimnames(rcp) <- dimnames(rdt) <- dimnames(RS)
FINE <- c("Tfh","Treg","CD8_T","CD4_conv","Macrophage","cDC","Monocyte")
FOC <- list(IB1 = c("FL_focus","cHL_hotspot","background"),
            RB1 = c("FL_follicular","interfollicular"))
keep <- GR[celltype %in% FINE & mapply(function(s,r) r %in% FOC[[s]], samp, bin_region) &
           cells >= 30, which = TRUE]
if (length(keep)) {
  fine <- rbindlist(lapply(keep, function(j) data.table(
    samp=GR$samp[j], bin_region=GR$bin_region[j], celltype=GR$celltype[j],
    cells=GR$cells[j], bins=GR$bins[j], fovs=GR$fovs[j], gene=GENES,
    cp10k=round(rcp[,j],3), det=round(rdt[,j],4))))
  fwrite(fine, file.path(OUT,"fine_celltype_pooled.csv.gz"), compress="gzip")
  cat("\nwrote ", file.path(OUT,"fine_celltype_pooled.csv.gz"),
      " -- POOLED effect sizes only, no region-unit test is possible for these\n", sep="")
  cat("\n=== fine-subtype groups available (cells >= 30)\n")
  print(GR[keep][order(samp, bin_region, -cells), .(samp, bin_region, celltype, cells, bins, fovs)])
}
writeLines(format(Sys.time()), file.path(OUT,"ct_de.done"))
cat("============================================================\nOK\n")
