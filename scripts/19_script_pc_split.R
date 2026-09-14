## ============================================================================
## 19_script_pc_split.R -- the Plasma cell cluster sits in TWO places. Are they the same thing?
##
## hc_clusters_high 24 carries the Plasma cell label and appears as two separate components in
## the embedding: an isolated island, and an offshoot hanging off the B mass. Step 19e found
## this label to be the single most neighbourhood-driven one in the object -- halo fold
## **12.69**, worse than the 8.4x step 18 measured for the InSituType PC call. The question is
## whether the halo lives in ONE of the two components.
##
## The prediction if it does:
##   island    real plasma cells -- high OWN plasma markers, light-chain saturated, argmax
##             agrees, and the label does not depend on who the neighbours are
##   offshoot  halo victims -- low own plasma markers, B markers intact, sitting close to real
##             plasma cells, argmax calls them B, and membership tracks the NEIGHBOUR's load
##
## If instead both components look alike, the halo is spread through the cluster and the whole
## label has to be treated as suspect rather than split.
##
## The components are found by k-means on the two embedding coordinates, which is legitimate
## here only because the question is "are these two visually separate blobs different?" -- the
## split is defined by position and then interrogated with evidence that is not position.
##
## Usage: sbatch --mem=250G --cpus-per-task=2 --time=3:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_pc_split.R
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); library(RANN)
  setDTthreads(2)
})
source("scripts/lib/lineage.R")
OUT <- "output/19_celltypes"; FIG <- file.path(OUT,"figures")
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}
UMPX <- 120.280945e-6*1000; MIN_DONOR <- 2L; D_MAX <- 20; set.seed(1L)
PC_G <- LINEAGE_SETS_COARSE$PC; B_G <- LINEAGE_SETS_COARSE$B
LC_G <- c("IGKC","IGLC1/2"); EXTRA <- c("JCHAIN","SEC11C","POU2AF1","TNFRSF17","CLEC4C")

d  <- fread(file.path(OUT,"hc_cells_labelled_high.csv.gz"), na.strings=c("NA",""))
md <- fread("output/16_regions/cells_meta_all.csv.gz", select=c("cell_id","fov_id","nCount_RNA"))
ce <- rbindlist(lapply(c("IB1","RB1"), function(S)
  fread(sprintf("output/18_leakage/%s_centroids.csv.gz",S), select=c("cell_id","cx_px","cy_px"))))
d <- ce[md[d, on=.(cell_id)], on=.(cell_id)]
d[, `:=`(X=cx_px*UMPX, Y=-cy_px*UMPX)]

seu <- readRDS(file.path(OUT,"hc_02_clustered.Rds"))
ca <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
setkey(d, cell_id); d <- d[J(colnames(ca))]
G <- intersect(c(PC_G,B_G,LC_G,EXTRA), rownames(ca))
sub <- as(as.matrix(ca[G,,drop=FALSE]),"dgCMatrix"); rownames(sub) <- G
rm(seu, ca); invisible(gc())
gs <- function(x) as.numeric(Matrix::colSums(sub[intersect(x,G),,drop=FALSE]))
d[, `:=`(k_PC=gs(PC_G), k_B=gs(B_G), lc=gs(LC_G), jchain=gs("JCHAIN"), clec4c=gs("CLEC4C"))]

## ---- 1. split the cluster on position ---------------------------------------
PCL <- d[cl == 24L]
if (!nrow(PCL)) stop("cluster 24 is not the Plasma cell cluster in this build", call.=FALSE)
## k = 3, NOT 2. With k = 2 the two DENSE blobs -- one at the top of the map, one on the
## edge of the B mass -- were merged into a single component while the diffuse scatter that
## covers the whole embedding became the other. That inverted the question being asked. The
## naming below is by COMPACTNESS rather than by distance to the B centroid, because the two
## dense blobs sit almost equidistant from it (3.25 vs 3.09) and a distance rule picked
## between them arbitrarily.
km <- kmeans(as.matrix(PCL[, .(U1,U2)]), centers=3, nstart=50)
PCL[, comp := km$cluster]
cen <- PCL[, .(n=.N, U1=round(mean(U1),2), U2=round(mean(U2),2),
               spread=round(mean(sqrt((U1-mean(U1))^2+(U2-mean(U2))^2)),2)), by=comp]
bcen <- d[celltype_cluster=="B cell", .(U1=median(U1), U2=median(U2))]
cen[, dist_to_B := round(sqrt((U1-bcen$U1)^2 + (U2-bcen$U2)^2), 2)]
setorder(cen, -spread)
cen[1, name := "diffuse scatter"]
rest <- cen[-1][order(-U2)]
cen[comp == rest$comp[1], name := "top island"]
cen[comp == rest$comp[2], name := "B-edge blob"]
cat("=== the two components ===\n"); print(cen)
PCL <- cen[, .(comp, name)][PCL, on=.(comp)]
cat(sprintf("\nk-means separation: between/total SS = %.3f\n", km$betweenss/km$totss))

## ---- 2. OWN evidence vs NEIGHBOUR evidence ---------------------------------
pool <- d[, .(cell_id, samp=merge_sample, X, Y, k_PC)]
don  <- pool[k_PC >= MIN_DONOR]
nn <- rbindlist(lapply(c("IB1","RB1"), function(S) {
  q <- pool[samp==S]; dd <- don[samp==S]
  r <- RANN::nn2(dd[, .(X,Y)], q[, .(X,Y)], k=2L)
  self <- dd$cell_id[r$nn.idx[,1]] == q$cell_id
  data.table(cell_id=q$cell_id,
             donPC_dist = ifelse(self, r$nn.dists[,2], r$nn.dists[,1]),
             donPC_k    = dd$k_PC[ifelse(self, r$nn.idx[,2], r$nn.idx[,1])])
}))
d <- nn[d, on=.(cell_id)]
PCL <- nn[PCL, on=.(cell_id)]

cat("\n=== OWN vs NEIGHBOUR plasma evidence, by component ===\n")
print(PCL[, .(cells=.N, pct_IB1=round(100*mean(merge_sample=="IB1"),1),
              nCount=round(median(nCount_RNA)),
              own_k_PC=round(mean(k_PC),2), own_k_B=round(mean(k_B),2),
              lightchain=round(mean(lc),1), JCHAIN=round(mean(jchain),2),
              argmax_PC=round(100*mean(lineage_argmax=="PC"),1),
              argmax_B =round(100*mean(lineage_argmax=="B"),1),
              no_call  =round(100*mean(!lineage_argmax %in% names(LINEAGE_SETS_COARSE)),1),
              donor_k_PC=round(mean(donPC_k, na.rm=TRUE),2),
              donor_dist=round(median(donPC_dist, na.rm=TRUE),1)), by=name][order(name)])

## step-07 plasma tier / distance to a real plasma cell, where available
pm <- fread("output/07_lightchain/IB1/plasma_markers.csv.gz",
            select=c("cell_id","tier","d_plasma_um","lc_total","pct_lc"))
x <- pm[PCL, on=.(cell_id), nomatch=NULL]
cat(sprintf("\n=== step-07 light-chain readout (available for %d of %d) ===\n", nrow(x), nrow(PCL)))
print(x[, .(cells=.N, med_d_plasma=round(median(d_plasma_um),1),
            med_lc=as.numeric(median(lc_total)), med_pct_lc=round(median(pct_lc),1),
            pct_plasma_tier=round(100*mean(tier=="plasma"),1),
            pct_halo_tier  =round(100*mean(tier=="halo"),1)), by=name][order(name)])

## ---- 3. the halo diagnostic, per component ---------------------------------
cat("\n=== HALO DIAGNOSTIC per component ===\n")
cat("Among cells with NO plasma marker of their own, P(cell is in this component) as the\n")
cat("nearest plasma-rich neighbour's load rises. Flat = a real population.\n\n")
z <- d[!is.na(donPC_dist) & donPC_dist <= D_MAX & k_PC == 0]
z <- PCL[, .(cell_id, name)][z, on=.(cell_id)]
z[, dq := cut(donPC_k, c(-Inf,2,3,5,9,Inf), labels=c("2","3","4-5","6-9","10+"))]
for (nm in c("top island","B-edge blob","diffuse scatter")) {
  t <- z[!is.na(dq), .(n=.N, p=round(100*mean(!is.na(name) & name==nm), 4)), by=dq][order(dq)]
  t[, fold := round(p/p[1], 2)]
  cat("--", nm, "\n"); print(t)
}
sc("pc_split_components.csv", PCL[, .(cell_id, merge_sample, name, comp, U1, U2, k_PC, k_B, lc,
                                      lineage_argmax, donPC_k, donPC_dist, nCount_RNA)])

SURF<-"#fcfcfb"
p <- ggplot(d[celltype_cluster %in% c("B cell","Plasma cell")][sample(.N, min(.N,250000L))],
            aes(U1,U2)) +
  geom_point(colour="#dcdcd6", size=.05, alpha=.4) +
  geom_point(data=PCL, aes(colour=name), size=.25, alpha=.7) +
  scale_colour_manual(values=c(`top island`="#c51b7d", `B-edge blob`="#00857a",
                               `diffuse scatter`="#e08214"), name=NULL) +
  guides(colour=guide_legend(override.aes=list(size=3,alpha=1))) +
  coord_fixed() +
  labs(title="The Plasma cell cluster in two positions",
       subtitle="grey = B cell and Plasma cell labels; coloured = the three components of cluster 24",
       x="UMAP 1", y="UMAP 2") +
  theme_minimal(11) + theme(plot.background=element_rect(fill=SURF,colour=NA),
    panel.grid.minor=element_blank(), plot.title=element_text(face="bold"), legend.position="top")
sf("pc_split.png", p, 7.5, 7)
cat("\ndone\n")
