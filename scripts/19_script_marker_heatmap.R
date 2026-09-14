## ============================================================================
## 19_script_marker_heatmap.R -- marker evidence for the step-19 labels, TWICE.
##
## The labels were BUILT from marker evidence, so a heatmap of the defining genes is
## circular by construction -- it can only confirm what the caller already assumed. This
## project's rule is that such a result is reported twice, with and without the defining
## genes, and only the second is evidence. So:
##
##   PANEL A  the DEFINING genes -- LINEAGE_SETS_COARSE, plus the pDC and Hodgkin-chemokine
##            modules and CR2/FCER2. Circular. It shows the caller did what it claims, and
##            nothing more.
##   PANEL B  HELD-OUT genes -- canonical markers of the same populations that no rule in
##            step 19 ever looked at. This is the evidence.
##
## The sharpest single cell in panel B is TNFRSF8 (CD30) against the Hodgkin chemokine
## niche. That niche was named on CCL17/CCL22/BATF3 alone; CD30 is the defining HRS antigen
## and was never used. Note the ceiling: per-cell CD30 calling is a closed avenue here
## (TNFRSF8 maxes at 6 counts, 78% singletons), so this is a POOLED cluster-level readout,
## which is exactly the level at which step 11 showed the chemokine field does work.
##
## cp10k is the POOLED share of the label's transcriptome (the step-11 convention), and
## detection rate is reported beside it because cp10k is compositional and fakes depression
## wherever something else is abundant.
##
## Usage: sbatch --mem=200G --cpus-per-task=2 --time=3:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_marker_heatmap.R
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); library(patchwork)
  setDTthreads(2)
})
source("scripts/lib/lineage.R")
OUT <- "output/19_celltypes"; FIG <- file.path(OUT,"figures")
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=170,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

DEFINING <- c(LINEAGE_SETS_COARSE,
  list(`pDC (module)`      = c("CLEC4C","LILRA4","IL3RA","IRF7","PLD4","GZMB","IRF8"),
       `Hodgkin (module)`  = c("CCL17","CCL22","BATF3","IL22RA2"),
       `FDC/naive (module)`= c("CR2","FCER2","IGHD","TCL1A","SELL")))
HELDOUT <- list(
  B            = c("CD22","TNFRSF13C","PAX5","BLK","VPREB3"),
  T            = c("CD2","CD5","LCK","ZAP70","IL7R"),
  NK           = c("KLRB1","PRF1","KLRF1","FGFBP2","SPON2"),
  myeloid      = c("CSF1R","MRC1","MSR1","TYROBP","AIF1","MNDA","CTSS"),
  endo         = c("RAMP2","CDH5","ESAM","PLVAP"),
  `stromal/FDC`= c("LUM","COL3A1","COL6A2","POSTN","FBLN1","SPARC","THY1"),
  plasma       = c("FKBP11","SLAMF7"),
  pDC          = c("TCF4","SPIB","BCL11A"),
  `Hodgkin/HRS`= c("TNFRSF8","CD274","PDCD1LG2","IL13"))

F_IN <- file.path(OUT,"hc_cells_final.csv.gz")
if (!file.exists(F_IN)) F_IN <- file.path(OUT,"hc_cells_labelled_high.csv.gz")
cat("reading ", F_IN, "\n", sep=""); d <- fread(F_IN, na.strings=c("NA",""))
## depth is needed for the matched by-sample panel
d <- fread("output/16_regions/cells_meta_all.csv.gz",
           select=c("cell_id","nCount_RNA"))[d, on=.(cell_id)]
seu <- readRDS(file.path(OUT,"hc_02_clustered.Rds"))
ca  <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
setkey(d,cell_id); d <- d[J(colnames(ca))]
tot <- seu$nCount_RNA
GEN <- unique(c(unlist(DEFINING), unlist(HELDOUT)))
miss <- setdiff(GEN, rownames(ca))
if (length(miss)) cat("not on panel, dropped: ", paste(miss, collapse=", "), "\n", sep="")
GEN <- intersect(GEN, rownames(ca))
sub <- as(as.matrix(ca[GEN,,drop=FALSE]),"dgCMatrix"); rownames(sub) <- GEN
rm(seu, ca); invisible(gc())

LAB <- d[!is.na(celltype_cluster), unique(celltype_cluster)]
prof <- rbindlist(lapply(LAB, function(L) {
  j <- which(d$celltype_cluster == L)
  data.table(label=L, gene=GEN,
             cp10k = as.numeric(Matrix::rowSums(sub[,j,drop=FALSE]))/sum(tot[j])*1e4,
             det   = as.numeric(Matrix::rowSums(sub[,j,drop=FALSE]>0))/length(j))
}))
sc("marker_profile_by_label.csv", prof)

## the same profile split by sample, so a label's marker evidence can be checked for
## CONSISTENCY between the two sections rather than assumed
MIN_CELLS <- 500L
key <- d[!is.na(celltype_cluster), .(cells=.N), by=.(label=celltype_cluster, samp=merge_sample)]
prof2 <- rbindlist(lapply(seq_len(nrow(key)), function(i) {
  L <- key$label[i]; S <- key$samp[i]
  j <- which(d$celltype_cluster == L & d$merge_sample == S)
  if (length(j) < MIN_CELLS) return(NULL)
  data.table(label=L, samp=S, gene=GEN,
             cp10k = as.numeric(Matrix::rowSums(sub[,j,drop=FALSE]))/sum(tot[j])*1e4,
             det   = as.numeric(Matrix::rowSums(sub[,j,drop=FALSE]>0))/length(j))
}))
sc("marker_profile_by_label_sample.csv", prof2)
cat("\n=== (label, sample) cells; those under ", MIN_CELLS, " are not drawn ===\n", sep="")
print(dcast(key, label ~ samp, value.var="cells", fill=0L))

## composition percentages, for the column labels
COMP  <- d[!is.na(celltype_cluster), .(n=.N), by=.(label=celltype_cluster)][, pct := 100*n/sum(n)]
COMPS <- d[!is.na(celltype_cluster), .(n=.N), by=.(label=celltype_cluster, samp=merge_sample)][
            , pct := 100*n/sum(n), by=samp]
lab_pct  <- setNames(sprintf("%s\n%.1f%%", COMP$label, COMP$pct), COMP$label)
samp_pct <- setNames(sprintf("%s  %.1f%%", COMPS$samp, COMPS$pct), paste(COMPS$label, COMPS$samp))

ORD <- c("B cell","Plasmablast","Plasma cell","T cell","NK cell","Myeloid","pDC",
         "Endothelial","Stromal/FDC","Hodgkin chemokine niche")
ORD <- intersect(ORD, LAB)
mk <- function(sets, title, sub_) {
  g <- rbindlist(lapply(names(sets), function(k)
        data.table(set=k, gene=intersect(sets[[k]], GEN))))
  g <- g[!duplicated(gene)]
  p <- prof[g, on=.(gene)][label %in% ORD]
  p[, z := (cp10k-mean(cp10k))/pmax(sd(cp10k),1e-9), by=gene]
  p[, label := factor(label, levels=ORD)]
  p[, gene := factor(gene, levels=rev(g$gene))]
  p[, set := factor(set, levels=names(sets))]
  ggplot(p, aes(label, gene, fill=z)) + geom_tile(colour="#fcfcfb", linewidth=.3) +
    geom_point(aes(size=det), shape=21, stroke=0, colour="#00000000", show.legend=TRUE) +
    scale_fill_gradient2(low="#184f95", mid="#f0efec", high="#d03b3b", midpoint=0,
                         limits=c(-1.5,1.5), oob=scales::squish, name="z across labels") +
    scale_size_area(max_size=2.2, name="detection rate") +
    scale_x_discrete(labels = function(v) lab_pct[v]) +
    facet_grid(set ~ ., scales="free_y", space="free_y", switch="y") +
    labs(title=title, subtitle=sub_, x=NULL, y=NULL) +
    theme_minimal(base_size=10) +
    theme(plot.background=element_rect(fill="#fcfcfb",colour=NA),
          panel.grid=element_blank(),
          axis.text.x=element_text(angle=40,hjust=1,colour="#52514e",size=7,lineheight=.9),
          axis.text.y=element_text(colour="#52514e",size=7),
          strip.text.y.left=element_text(angle=0,colour="#52514e",size=7,face="bold"),
          strip.placement="outside", plot.title=element_text(colour="#0b0b0b",face="bold"),
          plot.subtitle=element_text(colour="#52514e",size=8), legend.position="bottom")
}
pA <- mk(DEFINING, "A. DEFINING genes -- circular",
         "These genes named the labels. This panel can only show the caller did what it claims.\nColumn labels carry each label's share of the high-confidence set.")
pB <- mk(HELDOUT,  "B. HELD-OUT genes -- the evidence",
         "Canonical markers no step-19 rule ever looked at. TNFRSF8 (CD30) vs the Hodgkin niche is the sharpest test:\nthat niche was named on CCL17/CCL22/BATF3 alone.")
sf("marker_heatmap.png", pA + pB + plot_layout(widths=c(1,1)), 15, 11)

## ---- variant: IB1 / RB1 as minor columns ------------------------------------
## z is taken across ALL label x sample columns jointly, not within sample, so that a
## label whose marker profile differs between the two sections shows it. Per-sample
## composition is on the x ticks; a label's share of IB1 is not its share of RB1, and the
## object is 63.6% RB1, so the unsplit percentages are dominated by RB1 throughout.
mk2 <- function(sets, title, sub_) {
  g <- rbindlist(lapply(names(sets), function(k)
        data.table(set=k, gene=intersect(sets[[k]], GEN))))
  g <- g[!duplicated(gene)]
  p <- prof2[g, on=.(gene)][label %in% ORD]
  p[, z := (cp10k-mean(cp10k))/pmax(sd(cp10k),1e-9), by=gene]
  p[, label := factor(label, levels=ORD)]
  p[, gene := factor(gene, levels=rev(g$gene))]
  p[, set := factor(set, levels=names(sets))]
  ## x must be keyed on label|samp, not on samp. With x = samp every facet carries the same
  ## two discrete values, so scale_x_discrete can only give them ONE label and the per-sample
  ## percentages silently collapse to a bare "IB1"/"RB1". Free x scales then show each facet
  ## only its own two levels.
  p[, xkey := paste(label, samp, sep="|")]
  XL <- setNames(samp_pct[paste(COMPS$label, COMPS$samp)],
                 paste(COMPS$label, COMPS$samp, sep="|"))
  ggplot(p, aes(xkey, gene, fill=z)) + geom_tile(colour="#fcfcfb", linewidth=.3) +
    scale_fill_gradient2(low="#184f95", mid="#f0efec", high="#d03b3b", midpoint=0,
                         limits=c(-1.5,1.5), oob=scales::squish, name="z across label x sample") +
    scale_x_discrete(labels=function(v) XL[v]) +
    facet_grid(set ~ label, scales="free", space="free", switch="y") +
    labs(title=title, subtitle=sub_, x=NULL, y=NULL) +
    theme_minimal(base_size=10) +
    theme(plot.background=element_rect(fill="#fcfcfb",colour=NA), panel.grid=element_blank(),
          panel.spacing.x=unit(2,"pt"),
          axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour="#52514e",size=6),
          axis.text.y=element_text(colour="#52514e",size=7),
          strip.text.x=element_text(colour="#52514e",size=6.5,face="bold",angle=90),
          strip.text.y.left=element_text(angle=0,colour="#52514e",size=7,face="bold"),
          strip.placement="outside", plot.title=element_text(colour="#0b0b0b",face="bold"),
          plot.subtitle=element_text(colour="#52514e",size=8), legend.position="bottom")
}
cat("\n=== per-sample composition (% of that sample's high-confidence cells) ===\n")
print(dcast(COMPS, label ~ samp, value.var="pct")[, lapply(.SD, function(x)
      if (is.numeric(x)) round(x,2) else x)])
qA <- mk2(DEFINING, "A. DEFINING genes -- circular", "IB1 and RB1 as minor columns; x ticks carry each label's share of THAT sample.")
qB <- mk2(HELDOUT,  "B. HELD-OUT genes -- the evidence",
          "A label whose two columns disagree is not the same population in the two sections.")
sf("marker_heatmap_by_sample_pooled.png", qA + qB + plot_layout(widths=c(1,1)), 20, 11)

## ---- the CORRECTED by-sample panel ------------------------------------------
## The pooled version above is not a valid CROSS-SAMPLE comparison and the preserved
## *_v1_POOLED_SUPERSEDED.png files exist to show why. IB1's B label is 85% argmax-
## unconfirmed at a median 221 counts against RB1's 51.3% at 658, so a pooled cp10k compares
## two different populations and the paler IB1 column is that difference, not biology.
##
## Fixed two ways at once:
##   CONFIRMED  only cells argmax INDEPENDENTLY assigns to the label's family, which removes
##              label purity. pDC and the Hodgkin niche have no argmax family, so they cannot
##              be confirmed and are flagged rather than silently pooled.
##   DEPTH-STANDARDISED  cp10k computed within nCount deciles (deciles set on the pooled
##              confirmed cells of that label) and averaged with equal weight. Direct
##              standardisation: neither sample's depth distribution can drive the result.
##              Share falls steeply with depth -- 119.5 cp10k in decile 1 against 46.0 in
##              decile 10 -- so without this the shallower sample looks HIGHER, which is the
##              Simpson's-paradox reversal seen in b_discrepancy.
CONFIRM <- c("B cell"="B", "Plasmablast"="PC", "Plasma cell"="PC", "T cell"="T",
             "NK cell"="NK", "Myeloid"="myeloid", "Endothelial"="endo", "Stromal/FDC"="fibro")
MIN_DEC <- 30L
dm <- d[!is.na(celltype_cluster) & celltype_cluster %in% ORD]
dm[, fam := CONFIRM[celltype_cluster]]
dm[, confirmed := !is.na(fam) & lineage_argmax == fam]
cat("\n=== confirmable cells per label (argmax independently agrees) ===\n")
print(dcast(dm[, .(n=.N, conf=sum(confirmed)), by=.(celltype_cluster, merge_sample)],
            celltype_cluster ~ merge_sample, value.var=c("n","conf")))
prof3 <- rbindlist(lapply(ORD, function(L) {
  x <- dm[celltype_cluster == L]
  unmatched <- all(is.na(x$fam))
  if (!unmatched) x <- x[confirmed == TRUE]
  if (nrow(x) < 200) return(NULL)
  x[, dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                 include.lowest=TRUE, labels=FALSE)]
  per <- rbindlist(lapply(split(seq_len(nrow(x)), list(x$merge_sample, x$dec), drop=TRUE), function(ii) {
    if (length(ii) < MIN_DEC) return(NULL)
    j <- match(x$cell_id[ii], colnames(sub))
    data.table(samp=x$merge_sample[ii][1], dec=x$dec[ii][1], gene=GEN,
               cp10k=as.numeric(Matrix::rowSums(sub[,j,drop=FALSE]))/sum(tot[j])*1e4)
  }))
  ## a decile counts only if BOTH samples have enough cells in it
  keep <- per[, .(ns=uniqueN(samp)), by=.(dec,gene)][ns==2L, unique(dec)]
  per <- per[dec %in% keep]
  if (!nrow(per)) return(NULL)
  per[, .(label=L, unmatched=unmatched, cp10k=mean(cp10k)), by=.(samp, gene)]
}))
sc("marker_profile_matched.csv", prof3)
mk3 <- function(sets, title, sub_) {
  g <- rbindlist(lapply(names(sets), function(k)
        data.table(set=k, gene=intersect(sets[[k]], GEN))))[!duplicated(gene)]
  p <- prof3[g, on=.(gene)]
  p[, z := (cp10k-mean(cp10k))/pmax(sd(cp10k),1e-9), by=gene]
  p[, lab2 := fifelse(unmatched, paste0(label, " *"), label)]
  p[, lab2 := factor(lab2, levels=unique(p[order(match(label,ORD)), lab2]))]
  p[, gene := factor(gene, levels=rev(g$gene))]
  p[, set := factor(set, levels=names(sets))]
  p[, xkey := paste(lab2, samp, sep="|")]
  XL <- setNames(p$samp, p$xkey)
  ggplot(p, aes(xkey, gene, fill=z)) + geom_tile(colour="#fcfcfb", linewidth=.3) +
    scale_fill_gradient2(low="#184f95", mid="#f0efec", high="#d03b3b", midpoint=0,
                         limits=c(-1.5,1.5), oob=scales::squish, name="z across label x sample") +
    scale_x_discrete(labels=function(v) XL[v]) +
    facet_grid(set ~ lab2, scales="free", space="free", switch="y") +
    labs(title=title, subtitle=sub_, x=NULL, y=NULL) +
    theme_minimal(base_size=10) +
    theme(plot.background=element_rect(fill="#fcfcfb",colour=NA), panel.grid=element_blank(),
          panel.spacing.x=unit(2,"pt"),
          axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour="#52514e",size=6),
          axis.text.y=element_text(colour="#52514e",size=7),
          strip.text.x=element_text(colour="#52514e",size=6.5,face="bold",angle=90),
          strip.text.y.left=element_text(angle=0,colour="#52514e",size=7,face="bold"),
          strip.placement="outside", plot.title=element_text(colour="#0b0b0b",face="bold"),
          plot.subtitle=element_text(colour="#52514e",size=8), legend.position="bottom")
}
rA <- mk3(DEFINING, "A. DEFINING genes -- circular",
          "argmax-confirmed cells only, depth-standardised across nCount deciles.")
rB <- mk3(HELDOUT,  "B. HELD-OUT genes -- the evidence",
          paste0("Label purity and depth are both removed, so a difference between a label's two columns is about the CELLS, not the label.\n",
                 "* = pDC and Hodgkin chemokine niche cannot be argmax-confirmed (no family in LINEAGE_SETS_COARSE); those columns stay pooled."))
sf("marker_heatmap_by_sample.png", rA + rB + plot_layout(widths=c(1,1)), 20, 11)
cat("\n(the pooled version is preserved as *_v1_POOLED_SUPERSEDED.png -- see\n",
    " figures/README_marker_heatmap_versions.md)\n", sep="")

cat("\n=== HELD-OUT check: does each label top its own held-out set? ===\n")
for (k in names(HELDOUT)) {
  gs <- intersect(HELDOUT[[k]], GEN); if (!length(gs)) next
  x <- prof[gene %in% gs & label %in% ORD, .(mean_cp10k=mean(cp10k)), by=label][order(-mean_cp10k)]
  cat(sprintf("  %-13s top label: %-25s (%.2f)   2nd: %-25s (%.2f)\n", k,
              x$label[1], x$mean_cp10k[1], x$label[2], x$mean_cp10k[2]))
}
cat("\n=== TNFRSF8 (CD30) by label -- never used to name anything ===\n")
print(prof[gene=="TNFRSF8" & label %in% ORD,
           .(label, cp10k=round(cp10k,3), det=round(det,4))][order(-cp10k)])
cat("\ndone\n")
