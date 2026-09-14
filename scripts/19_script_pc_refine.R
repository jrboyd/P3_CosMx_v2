## ============================================================================
## 19_script_pc_refine.R -- is the Plasma cell halo flag real, and split the label.
##
## 19e flagged Plasma cell at halo fold 12.69, worse than the 8.4x step 18 measured for the
## InSituType PC call. But 19_script_pc_split.R showed halo fold tracks POSITIVELY with a
## component's own plasma evidence (top island 15.68 at light chain 26.1, B-edge blob 8.08 at
## 24.1, diffuse scatter 6.49 at 5.3). Under leakage the most-affected component would carry
## the LEAST independent evidence; here it carries the most. That is the signature of tissue
## architecture -- plasma cells form foci, which step 07 already formalised with its
## plasma/halo/B_clean tiers -- and it is exactly the confound 19e cannot resolve on its own.
##
## THE TEST. Take cells with ZERO plasma markers of their own sitting near a plasma-rich
## neighbour. They are either (a) genuine plasma cells whose markers dropped out, or (b) other
## cells wearing a PC label they did not earn. Those two have opposite fingerprints:
##
##                                    genuine, dropped out      contaminated non-plasma
##   own B markers                    LOW  (plasma downregulates MS4A1/CD79A)   HIGH
##   cell area                        LARGE (plasma cells are big)              typical
##   own light chain                  HIGH                                     LOW
##
## LIGHT CHAIN IS NOT A CLEAN ARBITER HERE and is reported with that caveat: step 18 showed
## the leaked mass IS light chain, so a cell's "own" IGKC may itself have bled in. B-marker
## load and cell area do not have that problem -- nothing about plasma-cell leakage raises a
## bystander's MS4A1 or makes it bigger -- so those two carry the verdict.
##
## THE SPLIT. Independently of the halo question, cluster 24 is three populations that
## resolution 0.8 failed to separate (k-means between/total SS = 0.833):
##   top island       19,520  76.3% IB1  -> Plasma cell
##   B-edge blob       5,056   7.7% IB1  -> Plasmablast: highest own PC markers (2.06), highest
##                                          JCHAIN (3.52), 52.1% step-07 plasma tier,
##                                          d_plasma 0.0um, AND highest B markers (2.80) with
##                                          38.7% argmax-B. Plasmacytoid differentiation of
##                                          the FL clone; it sits on the B-mass edge because
##                                          it is genuinely intermediate.
##   diffuse scatter   2,381  60.9% IB1  -> dropped: 272 median counts, no evidence of
##                                          anything, 41.8% step-07 halo tier, scattered over
##                                          the whole embedding.
##
## Usage: sbatch --mem=250G --cpus-per-task=2 --time=3:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_pc_refine.R
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(RANN); setDTthreads(2)
})
source("scripts/lib/lineage.R")
OUT <- "output/19_celltypes"; FIG <- file.path(OUT,"figures")
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}
D_MAX <- 20

## pc_split wrote per-cell components and the marker loads it computed
CP <- fread(file.path(OUT,"pc_split_components.csv"))
d  <- fread(file.path(OUT,"hc_cells_labelled_high.csv.gz"), na.strings=c("NA",""))
md <- fread("output/16_regions/cells_meta_all.csv.gz", select=c("cell_id","Area.um2","nCount_RNA"))
d  <- md[d, on=.(cell_id)]

## The donor geometry and marker loads for ALL cells: recomputed here the same way pc_split
## did, because the test needs every cell, not just cluster 24's.
UMPX <- 120.280945e-6*1000; MIN_DONOR <- 2L
ce <- rbindlist(lapply(c("IB1","RB1"), function(S)
  fread(sprintf("output/18_leakage/%s_centroids.csv.gz",S), select=c("cell_id","cx_px","cy_px"))))
d <- ce[d, on=.(cell_id)][, `:=`(X=cx_px*UMPX, Y=-cy_px*UMPX)]
suppressPackageStartupMessages(library(Seurat))
seu <- readRDS(file.path(OUT,"hc_02_clustered.Rds"))
ca <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
setkey(d, cell_id); d <- d[J(colnames(ca))]
G <- intersect(c(LINEAGE_SETS_COARSE$PC, LINEAGE_SETS_COARSE$B, "IGKC","IGLC1/2"), rownames(ca))
sub <- as(as.matrix(ca[G,,drop=FALSE]),"dgCMatrix"); rownames(sub) <- G
rm(seu, ca); invisible(gc())
gsx <- function(x) as.numeric(Matrix::colSums(sub[intersect(x,G),,drop=FALSE]))
d[, `:=`(k_PC=gsx(LINEAGE_SETS_COARSE$PC), k_B=gsx(LINEAGE_SETS_COARSE$B),
         lc=gsx(c("IGKC","IGLC1/2")))]
don <- d[k_PC >= MIN_DONOR]
nn <- rbindlist(lapply(c("IB1","RB1"), function(S) {
  q <- d[merge_sample==S]; dd <- don[merge_sample==S]
  r <- RANN::nn2(dd[, .(X,Y)], q[, .(X,Y)], k=2L)
  self <- dd$cell_id[r$nn.idx[,1]] == q$cell_id
  data.table(cell_id=q$cell_id, donPC_dist=ifelse(self,r$nn.dists[,2],r$nn.dists[,1]),
             donPC_k=dd$k_PC[ifelse(self,r$nn.idx[,2],r$nn.idx[,1])]) }))
d <- nn[d, on=.(cell_id)]

## ---- THE TEST ---------------------------------------------------------------
z <- d[k_PC == 0 & !is.na(donPC_dist) & donPC_dist <= D_MAX]
z[, dq := cut(donPC_k, c(-Inf,2,3,5,9,Inf), labels=c("2","3","4-5","6-9","10+"))]
z[, is_PC := celltype_cluster == "Plasma cell"]
cat("=== cells with ZERO plasma markers of their own, near a plasma-rich neighbour ===\n")
cat("If the PC-labelled ones are genuine plasma cells with dropout, their own B markers stay\n")
cat("LOW and their area stays LARGE as donor load rises. If they are contaminated bystanders,\n")
cat("B markers RISE.\n\n")
t1 <- z[!is.na(dq), .(n=.N, pct_PC_label=round(100*mean(is_PC),3)), by=dq][order(dq)]
t1[, fold := round(pct_PC_label/pct_PC_label[1], 2)]
print(t1)
cat("\n-- and, among the PC-LABELLED ones in each bin, what are they? --\n")
t2 <- z[!is.na(dq) & is_PC, .(n=.N, own_k_B=round(mean(k_B),2),
        area_um2=round(mean(Area.um2),1), own_lc=round(mean(lc),1),
        nCount=round(median(nCount_RNA))), by=dq][order(dq)]
print(t2)
cat("\n-- reference: NON-PC-labelled cells in the same bins --\n")
print(z[!is.na(dq) & !is_PC, .(n=.N, own_k_B=round(mean(k_B),2),
        area_um2=round(mean(Area.um2),1), own_lc=round(mean(lc),1)), by=dq][order(dq)])
sc("pc_halo_test.csv", merge(t1, t2, by="dq", suffixes=c("_all","_PClab")))

cat("\n=== the same gradient, stratified by the cell's OWN light chain ===\n")
cat("(caveat: light chain is the transcript step 18 showed actually leaks, so a high value\n")
cat(" here may itself be bled in. Read the B-marker and area columns above first.)\n\n")
z[, lcb := cut(lc, c(-Inf,0,2,10,Inf), labels=c("0","1-2","3-10",">10"))]
print(dcast(z[!is.na(dq) & !is.na(lcb), .(pct=round(100*mean(is_PC),3), n=.N), by=.(lcb,dq)][n>=200],
            lcb ~ dq, value.var="pct"))

## ---- THE SPLIT --------------------------------------------------------------
NEW <- c(`top island`="Plasma cell", `B-edge blob`="Plasmablast", `diffuse scatter`=NA)
d <- CP[, .(cell_id, pc_comp=name)][d, on=.(cell_id)]
d[, celltype_refined := celltype_cluster]
d[celltype_cluster == "Plasma cell", celltype_refined := NEW[pc_comp]]
d[celltype_cluster == "Plasma cell" & pc_comp == "diffuse scatter", cluster_status := "mixed"]
cat("\n=== refined labels ===\n")
print(d[, .N, by=celltype_refined][order(-N)])

## halo fold recomputed on the refined labels
for (L in c("Plasma cell","Plasmablast")) {
  y <- d[k_PC == 0 & !is.na(donPC_dist) & donPC_dist <= D_MAX]
  y[, dq := cut(donPC_k, c(-Inf,2,3,5,9,Inf), labels=c("2","3","4-5","6-9","10+"))]
  t <- y[!is.na(dq), .(n=.N, p=100*mean(celltype_refined == L, na.rm=TRUE)), by=dq][order(dq)]
  cat(sprintf("\nhalo fold, refined '%s': %.2f  (was 12.69 for the pooled label)\n",
              L, t$p[nrow(t)]/max(t$p[1],1e-9)))
  print(t[, .(dq, n, p=round(p,3))])
}
## halo_suspect is CLEARED for Plasma cell and not set for Plasmablast. The 12.69x flag was
## measuring plasma-cell spatial aggregation, not leakage: among cells with zero plasma
## markers of their own, PC-labelled cells carry HALF the B markers and 5-10x the light chain
## of non-PC-labelled cells IN THE SAME donor-load bins -- identical neighbourhood exposure,
## so the light chain is theirs -- and conditioning on their own light chain collapses the
## gradient from 12.62x to 1.62x. A raw halo fold is not a verdict for any spatially
## aggregated population; splitting the label RAISED the island's fold to 15.76 while making
## it more clearly plasma, which is the same point.
d[celltype_refined %in% c("Plasma cell","Plasmablast"), halo_suspect := FALSE]
d[, halo_note := NA_character_]
d[celltype_refined %in% c("Plasma cell","Plasmablast"),
  halo_note := "raw halo fold high (15.8/8.1) but shown to be spatial aggregation, not leakage: see pc_halo_test.csv"]

## hc_cells_final.csv.gz is the canonical per-cell table from here on. `celltype_cluster` is
## the refined label; the pre-split value is kept as celltype_cluster_pooled so the split is
## auditable and reversible.
setnames(d, "celltype_cluster", "celltype_cluster_pooled")
d[, celltype_cluster := celltype_refined]
sc("hc_cells_final.csv.gz",
   d[, .(cell_id, merge_sample, cl, pc_comp, celltype_cluster, celltype_cluster_pooled,
         cluster_status, is_niche, halo_suspect, halo_note,
         celltype_argmax, lineage_argmax, U1, U2)])
cat("\ndone\n")
