## ============================================================================
## 19_script_b_discrepancy.R -- why are B markers weaker in IB1's B cells than RB1's?
##
## The by-sample heatmap shows MS4A1/CD79A/BANK1 and the held-out CD22/PAX5/BLK all paler in
## IB1's B column than RB1's. Three things could cause that and they demand different
## responses, so they have to be separated rather than argued about:
##
##   (1) LABEL PURITY   IB1's "B cell" label contains more non-B cells, so the pooled share
##                      falls without any single cell changing. Tested by restricting to
##                      cells argmax INDEPENDENTLY calls B.
##   (2) DILUTION       IB1's B cells sit in tissue that is 49.4% myeloid+stromal against
##                      RB1's 11.1%, so they carry more foreign transcripts and every share
##                      falls mechanically. cp10k is compositional -- this is CLAUDE.md
##                      fact 6, and the project's rule is to gate on cp10k and detection
##                      rate agreeing in direction. Tested by comparing ABSOLUTE per-cell
##                      marker counts, which no amount of foreign transcript can lower.
##   (3) REAL           IB1's B cells genuinely express less. Only believable once (1) and
##                      (2) are excluded.
##
## THE CROSS-SAMPLE NEIGHBOURHOOD-MATCHED VERSION IS NOT AVAILABLE and must not be attempted:
## CLAUDE.md records that IB1's most B-dense quartile is less B-rich than RB1's least B-dense
## quartile, so only 2 of 89 bins find a composition-matched partner. Matching IB1 and RB1 B
## cells on neighbourhood is not identifiable on this data.
##
## SO THE DECISIVE TEST IS WITHIN IB1, which is the design this project prefers anyway (the
## CCR4 field gradient is the model). Within IB1 alone, does a B cell's marker signal fall as
## its neighbourhood becomes less B? If the SHARE falls while the ABSOLUTE COUNT stays flat,
## dilution is demonstrated outright, with no cross-sample comparison needed.
##
## Region is the unit: FOV-level bootstrap on every interval.
##
## Usage: sbatch --mem=250G --cpus-per-task=2 --time=4:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_b_discrepancy.R
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); library(RANN)
  setDTthreads(2)
})
source("scripts/lib/lineage.R")
OUT <- "output/19_celltypes"; FIG <- file.path(OUT,"figures")
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}
UMPX <- 120.280945e-6*1000; K_NB <- 30L; BOOT <- 400L; set.seed(1L)

B_DEF  <- c("MS4A1","CD79A","CD79B","BANK1")
B_HELD <- c("CD22","TNFRSF13C","PAX5","BLK","VPREB3")
FOREIGN<- unlist(LINEAGE_SETS_COARSE[c("T","myeloid","endo","fibro")], use.names=FALSE)

d <- fread(file.path(OUT,"hc_cells_labelled_high.csv.gz"), na.strings=c("NA",""))
md <- fread("output/16_regions/cells_meta_all.csv.gz", select=c("cell_id","fov_id","nCount_RNA"))
ce <- rbindlist(lapply(c("IB1","RB1"), function(S)
  fread(sprintf("output/18_leakage/%s_centroids.csv.gz",S), select=c("cell_id","cx_px","cy_px"))))
d <- ce[md[d, on=.(cell_id)], on=.(cell_id)]
d[, `:=`(X=cx_px*UMPX, Y=-cy_px*UMPX)]

seu <- readRDS(file.path(OUT,"hc_02_clustered.Rds"))
ca <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
setkey(d, cell_id); d <- d[J(colnames(ca))]
G <- intersect(c(B_DEF,B_HELD,FOREIGN), rownames(ca))
sub <- as(as.matrix(ca[G,,drop=FALSE]),"dgCMatrix"); rownames(sub) <- G
rm(seu, ca); invisible(gc())
gsum <- function(gs) as.numeric(Matrix::colSums(sub[intersect(gs,G),,drop=FALSE]))
d[, `:=`(b_def = gsum(B_DEF), b_held = gsum(B_HELD), foreign = gsum(FOREIGN))]
d[, `:=`(b_def_share = b_def/nCount_RNA, b_held_share = b_held/nCount_RNA,
         foreign_share = foreign/nCount_RNA,
         b_def_det = as.integer(b_def > 0), b_held_det = as.integer(b_held > 0))]

BC <- d[celltype_cluster == "B cell"]
cat(sprintf("cells labelled B cell: IB1 %d, RB1 %d\n",
            BC[merge_sample=="IB1",.N], BC[merge_sample=="RB1",.N]))

cat("\n=== (1) LABEL PURITY: is IB1's B label less pure? ===\n")
print(BC[, .(cells=.N, argmax_B = round(100*mean(lineage_argmax=="B"),1),
             argmax_resolves = round(100*mean(lineage_argmax %in% names(LINEAGE_SETS_COARSE)),1)),
         by=merge_sample])

cat("\n=== (2) SHARE vs ABSOLUTE COUNT vs DETECTION, all B-labelled cells ===\n")
f <- function(x) x[, .(cells=.N, nCount=round(median(nCount_RNA)),
        b_def_count = round(mean(b_def),3), b_def_share = round(1e4*mean(b_def_share),1),
        b_def_det   = round(100*mean(b_def_det),1),
        b_held_count= round(mean(b_held),3), b_held_share= round(1e4*mean(b_held_share),1),
        b_held_det  = round(100*mean(b_held_det),1),
        foreign_share = round(1e4*mean(foreign_share),1)), by=merge_sample][order(merge_sample)]
print(f(BC))
cat("\n-- restricted to cells argmax INDEPENDENTLY calls B (removes cause 1) --\n")
print(f(BC[lineage_argmax=="B"]))
cat("\n-- and depth-matched: within nCount deciles, argmax-B only --\n")
x <- BC[lineage_argmax=="B"]
x[, dec := cut(nCount_RNA, quantile(nCount_RNA,0:10/10), include.lowest=TRUE, labels=FALSE)]
print(dcast(x[, .(count=round(mean(b_def),2), share=round(1e4*mean(b_def_share),1),
                  det=round(100*mean(b_def_det),1)), by=.(merge_sample, dec)],
            dec ~ merge_sample, value.var=c("count","share","det")))

cat("\n=== (3) WITHIN-IB1 GRADIENT: the identifiable test ===\n")
cat("For each B-labelled cell, the fraction of its 30 nearest neighbours that are NOT B.\n")
cat("If SHARE falls with neighbourhood while COUNT stays flat, dilution is demonstrated.\n\n")
d[, is_B := as.integer(celltype_cluster == "B cell")]
nb <- rbindlist(lapply(c("IB1","RB1"), function(S) {
  z <- d[merge_sample==S]
  r <- RANN::nn2(z[, .(X,Y)], z[, .(X,Y)], k=K_NB+1L)
  nonB <- 1 - (rowSums(matrix(z$is_B[r$nn.idx[,-1]], nrow=nrow(z))) / K_NB)
  data.table(cell_id=z$cell_id, nonB_frac=nonB)
}))
d <- nb[d, on=.(cell_id)]
BC <- d[celltype_cluster == "B cell" & lineage_argmax == "B"]
BC[, nbin := cut(nonB_frac, c(-.01,.1,.25,.5,.75,.9,1.01),
                 labels=c("<10%","10-25%","25-50%","50-75%","75-90%",">90%"))]
g <- BC[, .(cells=.N, nCount=round(median(nCount_RNA)),
            b_count = round(mean(b_def),2), b_share = round(1e4*mean(b_def_share),1),
            b_det   = round(100*mean(b_def_det),1),
            foreign_share = round(1e4*mean(foreign_share),1)), by=.(merge_sample, nbin)][
            order(merge_sample, nbin)]
print(g[cells >= 200])
sc("b_discrepancy_gradient.csv", g)

## FOV-bootstrapped slope of share and of count on the non-B fraction, within each sample
fovs <- unique(BC[, .(merge_sample, fov_id)])
sl <- function(z) c(coef(lm(1e4*b_def_share ~ nonB_frac, z))[2], coef(lm(b_def ~ nonB_frac, z))[2])
res <- rbindlist(lapply(c("IB1","RB1"), function(S) {
  z <- BC[merge_sample==S]; fv <- unique(z$fov_id)
  o <- sl(z)
  bs <- replicate(BOOT, { f <- data.table(fov_id=sample(fv, length(fv), replace=TRUE))
                          sl(f[z, on=.(fov_id), nomatch=NULL, allow.cartesian=TRUE]) })
  data.table(samp=S, metric=c("share (cp10k)","absolute count"), slope=round(o,3),
             lo=round(apply(bs,1,quantile,.025),3), hi=round(apply(bs,1,quantile,.975),3))
}))
cat("\n=== slope on the non-B neighbour fraction, FOV-bootstrapped ===\n")
print(res)
sc("b_discrepancy_slopes.csv", res)
cat("\nA negative SHARE slope with a flat/positive COUNT slope is dilution, not lower expression.\n")

p <- ggplot(g[cells>=200], aes(nbin, group=merge_sample, colour=merge_sample)) +
  geom_line(aes(y=b_share/b_share[1]*100), linewidth=.7) +
  geom_line(aes(y=b_count/b_count[1]*100), linewidth=.7, linetype="22") +
  facet_wrap(~merge_sample) +
  scale_colour_manual(values=c(IB1="#2a78d6", RB1="#eb6834"), guide="none") +
  labs(title="B-marker signal against neighbourhood composition, within sample",
       subtitle="solid = share of the cell's transcriptome (cp10k), dashed = absolute count. Both indexed to 100 at the leftmost bin.\nShare falling while count holds is dilution by foreign transcripts, not lower B expression.",
       x="fraction of the 30 nearest neighbours that are NOT B", y="relative to the most B-dense bin (%)") +
  theme_minimal(11) + theme(plot.background=element_rect(fill="#fcfcfb",colour=NA),
    panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
    plot.subtitle=element_text(size=8), axis.text.x=element_text(angle=30,hjust=1))
sf("b_discrepancy.png", p, 11, 5)
cat("\ndone\n")
