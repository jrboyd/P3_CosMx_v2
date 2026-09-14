## Two checks on step 14.
##
## A. IS THE "DOWN IN NEIGHBORHOODS" SET REAL, OR AMBIENT BACKGROUND? PHKG2, A1BG, DUSP2
##    and GPX1 come out depressed inside neighborhoods in nearly every cluster, and PHKG2
##    sits at ~100 cp10k in out-tissue -- about 1% of all transcripts for a phosphorylase
##    kinase subunit, which is not credible. If out-tissue simply carries more ambient
##    signal, the NEGATIVE PROBES will show the same in/out shift. They are not biology,
##    so any shift in them is pure background and calibrates the gene-level result.
##
## B. HOW BIASED IS THE 154 -> 68 NEIGHBORHOOD ATTRITION? Step 07 built neighborhoods on
##    the PRE-QC object (2.54M cells); pearson clusters exist only post-QC (400k), so 86
##    neighborhoods have too little surviving tissue. Report which classes were lost.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);setDTthreads(1)})
OUT <- "output/14_lightchain_in_out/IB1"
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote ",file.path(OUT,n),"\n",sep="")}

## ---- B: attrition (cheap, do it first) --------------------------------------
nbh <- fread("output/07_lightchain/IB1/lightchain_neighborhoods_r40um.csv")
comp <- fread(file.path(OUT,"neighborhood_composition.csv"))
kept <- as.integer(sub("^nb_","",comp$region))
nbh[, retained := nb %in% kept]
cat("=== neighborhood attrition, 154 -> ", sum(nbh$retained), "\n", sep="")
att <- nbh[, .(total=.N, retained=sum(retained), pct=round(100*mean(retained),1),
               med_cells_kept   = as.numeric(median(n_cells[retained])),
               med_cells_lost   = as.numeric(median(n_cells[!retained])),
               med_pct_i_kept   = round(median(pct_cluster_i[retained]),1),
               med_pct_i_lost   = round(median(pct_cluster_i[!retained]),1)), by=class]
print(att[order(-total)]); sc("neighborhood_attrition.csv", att)
cat("\noverall: kept median n_cells ", median(nbh[retained==TRUE]$n_cells),
    " vs lost ", median(nbh[retained==FALSE]$n_cells), "\n", sep="")
cat("kappa_frac kept ", round(median(nbh[retained==TRUE]$kappa_frac),3),
    " vs lost ", round(median(nbh[retained==FALSE]$kappa_frac),3), "\n", sep="")
w <- wilcox.test(nbh[retained==TRUE]$n_cells, nbh[retained==FALSE]$n_cells)
cat("size difference kept vs lost: p = ", signif(w$p.value,3), "\n", sep="")

## ---- A: negative probes in vs out --------------------------------------------
o <- readRDS("output/03_filter_and_process/IB1_03_HCA_cell_typed.Rds")
rn <- rownames(o@meta.data)
md <- as.data.table(o@meta.data[, c("x_slide_mm","y_slide_mm","nCount_RNA",
                                    "nCount_negprobes","inSitu_clusters"), drop=FALSE])
md[, cell_id := rn]
NEG <- as(LayerData(o[["negprobes"]], layer="counts"), "dgCMatrix")
RNAM <- as(LayerData(o[["RNA"]], layer="counts"), "dgCMatrix")
rm(o); invisible(gc())
md <- md[!is.na(x_slide_mm) & nCount_RNA > 0]
lc <- fread("output/07_lightchain/IB1/lightchain_cells_r40um.csv")
mem <- lc[nb %in% nbh$nb][cell_id %in% md$cell_id]
q <- RANN::nn2(as.matrix(mem[, .(x_slide_mm*1000, y_slide_mm*1000)]),
               as.matrix(md[, .(x_slide_mm*1000, y_slide_mm*1000)]), k=1)
md[, d := as.numeric(q$nn.dists[,1])]
md[, grp := fifelse(d <= 40, "in", fifelse(d > 100, "out", "buffer"))]
md <- md[grp %in% c("in","out") & inSitu_clusters != "i" & !(cell_id %in% mem$cell_id)]
i <- match(md$cell_id, colnames(NEG))
cat("\n=== negative probes and total signal, in vs out\n")
res <- md[, .(cells=.N, med_nCount=as.numeric(median(nCount_RNA)),
              neg_per_cell = round(mean(nCount_negprobes),3),
              neg_cp10k = round(1e4*sum(nCount_negprobes)/sum(nCount_RNA),2)), by=grp]
print(res); sc("negprobe_in_out.csv", res)
r <- res[grp=="out", neg_cp10k] / res[grp=="in", neg_cp10k]
cat(sprintf("out:in negative-probe rate = %.2fx\n", r))
FOC <- c("PHKG2","A1BG","DUSP2","GPX1","CXCL13","CLU","VCAM1","MS4A1","CD79A")
gi <- match(FOC, rownames(RNAM)); names(gi) <- FOC
sub <- RNAM[gi[!is.na(gi)], i, drop=FALSE]
tot <- md$nCount_RNA
gg <- rbindlist(lapply(rownames(sub), function(g) {
  v <- as.numeric(sub[g,])
  data.table(gene=g,
             cp10k_in  = round(1e4*sum(v[md$grp=="in"]) /sum(tot[md$grp=="in"]),2),
             cp10k_out = round(1e4*sum(v[md$grp=="out"])/sum(tot[md$grp=="out"]),2),
             pct_det_in  = round(100*mean(v[md$grp=="in"]>0),2),
             pct_det_out = round(100*mean(v[md$grp=="out"]>0),2),
             max_count = max(v))
}))
gg[, log2fc := round(log2((cp10k_in+.1)/(cp10k_out+.1)),2)]
cat("\n=== focus genes: is a high cp10k carried by real counts or by breadth?\n")
print(gg[order(log2fc)]); sc("focus_gene_detection.csv", gg)
cat("\nIf PHKG2/A1BG sit at high cp10k with LOW max counts and HIGH detection breadth,\n",
    "they behave like ambient background rather than expression.\n", sep="")
