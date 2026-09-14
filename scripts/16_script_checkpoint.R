## ============================================================================
## 16_script_checkpoint.R -- PD-1 / PD-L1 / PD-L2, in detail.
##
## THE CALIBRATION THIS RESTS ON (scripts/sbatch/qc_16_negprobe_floor.R):
##   * At >= 1 count these genes are near the noise floor: CD274 is 2.7x the mean blank
##     probe, PDCD1 2.9x, PDCD1LG2 2.6x. A >=1 call is not a call.
##   * At >= 2 counts the blank rate collapses to 0.05% and the same genes are 25-31x it --
##     specific, but sparse (about 1-2% of cells). THR = 2 is therefore not a tuning choice,
##     it is the smallest threshold at which the measurement means anything.
##   * Depth predicts positivity at AUC 0.64-0.76 for EVERY gene, IGKC (0.693) and CD3E
##     (0.695) included. So AUC ~0.70 is a property of thresholding sparse counts here, not
##     evidence a gene is noise -- but it does mean any comparison across regions of unequal
##     depth must be depth-matched. Both are done below.
##
## FOUR CROSS-CUTS, and one of them cannot carry a p-value:
##   1. region x cell type, pooled cp10k + detection            effect sizes, both samples
##   2. within B / T / myeloid, bin as the unit, depth-matched   REGION-UNIT TESTS
##   3. per-cell positivity at >= 2 counts, with the negprobe    calibrated prevalence
##      false-positive rate printed beside every number
##   4. spatial PD-1 -> PD-L1 pairing against a density- AND     is the axis engaged?
##      celltype-matched permutation null
##
## Cross-cut 1 is pooled over cells and therefore pseudo-replicated; it is reported as an
## effect size only. Cross-cut 2 is where the p-values live. The fine subtypes (Tfh, Treg,
## Macrophage, cDC) are under one cell per bin and can never support cross-cut 2 -- said here
## so no table implies otherwise.
##
## Tfh is the validated PD-1 source: 16b called it WITHOUT PDCD1 (CXCL13/ICOS/CXCR5/BTLA/
## CD40LG) and it still came top of its siblings on held-out PDCD1. Treg was NOT validated
## (rank 2 behind Tfh on CTLA4) and is flagged accordingly wherever it appears.
##
## Usage: Rscript scripts/16_script_checkpoint.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(RANN); library(data.table); library(ggplot2); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "matched"

REG_F <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT_F  <- "output/16_celltypes/cell_celltypes.csv.gz"
CNT_F <- "output/16_regions/panel_counts_all.Rds"
NEG_F <- "output/16_regions/negprobe_counts_all.Rds"
CTPB  <- file.path("output/16_ct_pseudobulk", TIER)
OUT   <- file.path("output/16_checkpoint", TIER)
FIG   <- file.path(OUT, "figures")
for (f in c(REG_F, CT_F, CNT_F, NEG_F)) if (!file.exists(f)) stop("missing ", f, call.=FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

THR      <- 2L        # calibrated, not chosen -- see the header
THR_SWEEP<- c(1L,2L,3L)
PRIMARY  <- c("PDCD1","CD274","PDCD1LG2")
CONTEXT  <- c("CTLA4","LAG3","TIGIT","HAVCR2","TOX","ENTPD1","TNFRSF9","BTLA","VSIR",
              "LGALS9","IDO1","ICOS","CXCL13","FOXP3","CD8A","GZMB","IFNG","CCR4",
              "CXCL9","CXCL10","CXCL11","CD80","CD86")
N_PERM   <- 100L; MIN_POS <- 30L; SEED <- 1L
FOC <- list(IB1 = c("FL_focus","cHL_hotspot","background"),
            RB1 = c("FL_follicular","interfollicular"))
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
DIV_HI<-"#2a78d6"; DIV_LO<-"#eb6834"

cat("============================================================\n")
cat("step 16e -- PD-1 / PD-L1 / PD-L2, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

rg <- fread(REG_F, select=c("cell_id","samp","bin_region","x_slide_mm","y_slide_mm",
                            "is_defect","unbinned","fov_id","nCount_RNA"))
ct <- fread(CT_F, select=c("cell_id","celltype","lineage16"))
d  <- merge(rg, ct, by="cell_id", sort=FALSE)[is_defect==FALSE & unbinned==FALSE]
d  <- d[mapply(function(s,r) r %in% FOC[[s]], samp, bin_region)]
M  <- readRDS(CNT_F); NP <- readRDS(NEG_F)
jj <- match(d$cell_id, colnames(M)); stopifnot(!anyNA(jj))
M <- M[, jj, drop=FALSE]; NP <- NP[, jj, drop=FALSE]
GOI <- intersect(c(PRIMARY, CONTEXT), rownames(M))
cat(sprintf("cells %d | genes %d | negprobes %d\n", nrow(d), length(GOI), nrow(NP)))
print(d[, .N, by=.(samp, bin_region)][order(samp,-N)])

## ---- 1. pooled cp10k + detection, region x cell type -------------------------
pool <- rbindlist(lapply(unique(paste(d$samp, d$bin_region, d$celltype, sep="|")), function(k) {
  p <- strsplit(k, "|", fixed=TRUE)[[1]]
  ii <- which(d$samp==p[1] & d$bin_region==p[2] & d$celltype==p[3])
  if (length(ii) < 30) return(NULL)
  tot <- sum(d$nCount_RNA[ii])
  data.table(samp=p[1], bin_region=p[2], celltype=p[3], cells=length(ii), gene=GOI,
             cp10k = round(1e4*as.numeric(Matrix::rowSums(M[GOI,ii,drop=FALSE]))/tot, 3),
             pct_det = round(100*as.numeric(Matrix::rowMeans(M[GOI,ii,drop=FALSE] > 0)), 3),
             pct_pos = round(100*as.numeric(Matrix::rowMeans(M[GOI,ii,drop=FALSE] >= THR)), 3),
             med_depth = as.numeric(median(d$nCount_RNA[ii])))
}))
sc("checkpoint_pooled.csv", pool)
cat(sprintf("\n=== PD-1 axis by cell type and region, cp10k (pooled -- EFFECT SIZES, no test)\n"))
for (g in PRIMARY) {
  cat(sprintf("\n-- %s\n", g))
  print(dcast(pool[gene==g & celltype %in% c("B","Tfh","Treg","CD8_T","CD4_conv",
              "Macrophage","cDC","Monocyte","T_unresolved")],
              celltype ~ samp + bin_region, value.var="cp10k"))
}

## ---- 2. region-unit tests within coarse cell type ---------------------------
if (file.exists(file.path(CTPB,"pb_binct_sum.Rds"))) {
  S <- readRDS(file.path(CTPB,"pb_binct_sum.Rds")); D <- readRDS(file.path(CTPB,"pb_binct_det.Rds"))
  GB <- fread(file.path(CTPB,"binct_meta.csv.gz")); setkey(GB,grp_bin); GB <- GB[colnames(S)]
  GB[, `:=`(log_depth=log(nCount/cells), log_cells=log(cells))]
  CP <- as.matrix(S[GOI,,drop=FALSE]) %*% Diagonal(x=1e4/pmax(Matrix::colSums(S),1))
  DT <- as.matrix(D[GOI,,drop=FALSE]) %*% Diagonal(x=1/pmax(GB$cells,1))
  dimnames(CP) <- dimnames(DT) <- list(GOI, GB$grp_bin)
  mtch <- function(ia,ib,cal=0.25) {
    V <- c("log_depth","log_cells")
    Mx <- rbind(as.matrix(GB[ia,..V]), as.matrix(GB[ib,..V]))
    mu <- colMeans(Mx); sdv <- apply(Mx,2,stats::sd); sdv[sdv==0] <- 1
    za <- scale(as.matrix(GB[ia,..V]),mu,sdv); zb <- scale(as.matrix(GB[ib,..V]),mu,sdv)
    set.seed(SEED); used <- rep(FALSE,length(ib)); pk <- rep(NA_integer_,length(ia))
    for (i in sample.int(length(ia))) { dd <- sqrt(colSums((t(zb)-za[i,])^2)); dd[used] <- Inf
      j <- which.min(dd); if (is.finite(dd[j]) && dd[j] <= cal) { pk[i] <- j; used[j] <- TRUE } }
    ok <- !is.na(pk); list(a=ia[ok], b=ib[pk[ok]])
  }
  CON <- list(
    list("IB1 FL vs IB1 bg", "IB1","FL_focus","IB1","background"),
    list("IB1 hot vs IB1 bg","IB1","cHL_hotspot","IB1","background"),
    list("RB1 foll vs inter", "RB1","FL_follicular","RB1","interfollicular"),
    list("IB1 FL vs RB1 foll","IB1","FL_focus","RB1","FL_follicular"))
  RU <- rbindlist(lapply(c("B","T","myeloid"), function(cc) {
   rbindlist(lapply(CON, function(z) {
    ia <- which(GB$samp==z[[2]] & GB$bin_region==z[[3]] & GB$coarse==cc)
    ib <- which(GB$samp==z[[4]] & GB$bin_region==z[[5]] & GB$coarse==cc)
    if (min(length(ia),length(ib)) < 12L) return(NULL)
    m <- mtch(ia,ib); if (length(m$a) < 12L) return(NULL)
    rbindlist(lapply(GOI, function(g) {
      w <- suppressWarnings(wilcox.test(CP[g,m$a], CP[g,m$b]))
      data.table(contrast=z[[1]], celltype=cc, gene=g, n_a=length(m$a), n_b=length(m$b),
                 cp10k_a=mean(CP[g,m$a]), cp10k_b=mean(CP[g,m$b]),
                 det_a=mean(DT[g,m$a]), det_b=mean(DT[g,m$b]), p=w$p.value)
    }), fill=TRUE)
   }), fill=TRUE)
  }), fill=TRUE)
  if (nrow(RU)) {
    RU[, q := p.adjust(p,"BH"), by=.(contrast,celltype)]
    RU[, log2fc := log2((cp10k_a+.1)/(cp10k_b+.1))]
    RU[, log2fc_det := log2((det_a+1e-4)/(det_b+1e-4))]
    RU[, concordant := sign(log2fc)==sign(log2fc_det) & abs(log2fc_det) > .3]
    sc("checkpoint_region_unit.csv", RU)
    cat("\n=== REGION-UNIT tests, depth-matched, within cell type (q<0.05 & concordant)\n")
    h <- RU[q < 0.05 & concordant][order(contrast, celltype, -abs(log2fc))]
    if (nrow(h)) print(h[, .(contrast, celltype, gene, n_a, n_b, cp10k_a=round(cp10k_a,2),
                             cp10k_b=round(cp10k_b,2), log2fc=round(log2fc,2), q=signif(q,2))])
    else cat("  none. No checkpoint gene differs between regions within a cell type\n",
             "  once depth is matched and cp10k/detection are required to agree.\n", sep="")
    cat("\n  PD-1 axis specifically (all results, significant or not):\n")
    print(RU[gene %in% PRIMARY][order(contrast,celltype,gene),
             .(contrast, celltype, gene, cp10k_a=round(cp10k_a,2), cp10k_b=round(cp10k_b,2),
               log2fc=round(log2fc,2), q=signif(q,2), concordant)])
  }
}

## ---- 3. per-cell positivity, with the blank-probe floor beside it ------------
pos <- rbindlist(lapply(THR_SWEEP, function(t) rbindlist(lapply(
  unique(paste(d$samp,d$bin_region,sep="|")), function(k) {
    p <- strsplit(k,"|",fixed=TRUE)[[1]]
    ii <- which(d$samp==p[1] & d$bin_region==p[2])
    nd <- as.numeric(Matrix::rowMeans(NP[,ii,drop=FALSE] >= t))
    data.table(thr=t, samp=p[1], bin_region=p[2], cells=length(ii), gene=GOI,
      pct_pos = round(100*as.numeric(Matrix::rowMeans(M[GOI,ii,drop=FALSE] >= t)),4),
      neg_floor_pct = round(100*mean(nd),4), neg_worst_pct = round(100*max(nd),4))
  }))))
pos[, ratio_vs_floor := round(pct_pos/pmax(neg_floor_pct,1e-9),1)]
sc("positivity_by_region.csv", pos)
cat(sprintf("\n=== per-cell positivity at >= %d counts, against the blank-probe floor\n", THR))
print(dcast(pos[thr==THR & gene %in% PRIMARY], samp+bin_region+neg_floor_pct ~ gene,
            value.var=c("pct_pos","ratio_vs_floor")))

## ---- 4. is the axis engaged? spatial PD-1 -> PD-L1 pairing -------------------
## Observed median distance from each PD-1+ cell to the nearest PD-L1+ cell, against a null
## that permutes the PD-L1 label WITHIN strata of (local density decile x cell type). Both
## strata matter: density because a dense region shortens every distance, and cell type
## because PD-L1 sits on macrophages and DCs, so a label shuffled across cell types would
## move it out of the myeloid compartment and manufacture an effect. This is step 10's
## adjacency null (10_script_hrs_call.R:234-266) with the cell-type stratum added.
pair1 <- function(S, R) {
  ii <- which(d$samp==S & d$bin_region==R)
  if (length(ii) < 500) return(NULL)
  xy <- as.matrix(d[ii, .(x_slide_mm*1000, y_slide_mm*1000)])
  p1 <- as.vector(M["PDCD1", ii] >= THR)
  for (lig in intersect(c("CD274","PDCD1LG2"), rownames(M))) {
    p2 <- as.vector(M[lig, ii] >= THR)
    if (sum(p1) < MIN_POS || sum(p2) < MIN_POS) next
    nn <- RANN::nn2(xy, k = min(21L, length(ii)))
    dens <- nn$nn.dists[, ncol(nn$nn.dists)]          # distance to the 20th neighbour
    strat <- paste(cut(dens, quantile(dens, seq(0,1,.1)), include.lowest=TRUE),
                   d$celltype[ii])
    obs <- median(RANN::nn2(xy[p2,,drop=FALSE], xy[p1,,drop=FALSE], k=1)$nn.dists[,1])
    set.seed(SEED)
    nul <- vapply(seq_len(N_PERM), function(b) {
      pp <- unlist(lapply(split(seq_along(ii), strat), function(k)
        if (!length(k)) integer(0) else sample(k, size=sum(p2[k]))), use.names=FALSE)
      if (length(pp) < MIN_POS) return(NA_real_)
      median(RANN::nn2(xy[pp,,drop=FALSE], xy[p1,,drop=FALSE], k=1)$nn.dists[,1])
    }, 0)
    nul <- nul[is.finite(nul)]; if (!length(nul)) next
    assign("PAIR", rbind(get0("PAIR", ifnotfound=data.table()), data.table(
      samp=S, bin_region=R, ligand=lig, cells=length(ii), n_pd1=sum(p1), n_lig=sum(p2),
      obs_um=round(obs,1), null_um=round(median(nul),1), null_sd=round(sd(nul),2),
      ratio_null_obs=round(median(nul)/obs,3),
      p_emp=(1+sum(nul <= obs))/(length(nul)+1))), envir=.GlobalEnv)
  }
  invisible(NULL)
}
PAIR <- data.table()
for (S in names(FOC)) for (R in FOC[[S]]) pair1(S, R)
if (nrow(PAIR)) {
  PAIR[, verdict := fifelse(p_emp < 0.05 & ratio_null_obs > 1, "closer than chance",
                    fifelse(p_emp > 0.95 & ratio_null_obs < 1, "further than chance", "no signal"))]
  sc("spatial_pairing.csv", PAIR)
  cat(sprintf("\n=== spatial PD-1 -> ligand pairing at >= %d counts\n", THR))
  cat("    ratio > 1 and p < 0.05 = PD-1+ cells sit CLOSER to ligand+ cells than a\n")
  cat("    density- and celltype-matched shuffle allows\n")
  print(PAIR[order(samp, bin_region, ligand)])
} else cat("\n=== spatial pairing: too few positive cells in every region at this threshold\n")
writeLines(format(Sys.time()), file.path(OUT,"checkpoint.done"))
cat("============================================================\nOK\n")
