## ============================================================================
## 19_script_identify_clusters.R -- what ARE the unlabelable low-resolution clusters?
##
## SUPERSEDED, KEPT AS THE RECORD OF A DIAGNOSIS. This was a one-off run against
## hc_clusters_low, and its two findings are now part of the pipeline proper:
##   cluster 4  ->  the Hodgkin compartment (CCL17 36.7x, CCL22 30.4x, BATF3 15.3x, 85.8% IB1)
##   cluster 7  ->  pDC (CLEC4C 11.4x, LILRA4 7.6x, IL3RA 4.8x, IRF7 11.0x)
## Both are now named by the cluster-level module rules in 19_script_b_substructure.R and
## 19_script_cluster_labels.R, and both were validated afterwards on held-out markers.
##
## Do not re-run it expecting current cluster ids: it predates the per-cell low-information
## floor added to 19a, so the object was re-clustered and the ids in its output no longer
## refer to the same cells. It is retained because it is the evidence behind two labels that
## the coarse marker sets could never have produced, and because its QC profile section is
## the worked example of a cluster looking like borderline QC (dim, noisy, sample-specific)
## while actually being real biology.
##
## hc_clusters_low 4 and 7 both fail 19d's purity/margin test and carry no label. The
## question is whether that is a QC failure -- dim, noisy tissue that scraped past the
## per-cell floor -- or a real population the coarse marker sets cannot name.
##
## Their headline QC numbers already point opposite ways against the object baselines
## (frac_resolved 0.579, median nCount 541, pooled negprobe rate 0.000442):
##   cluster 4   nCount 414 (0.77x), negprobe 1.40x, 85.8% IB1, resolves 0.426
##   cluster 7   nCount 612 (1.13x), negprobe 1.15x, 53.2% IB1, resolves 0.464
## Cluster 7 is BETTER sequenced than the object average, so whatever is wrong with it is
## not depth. This script asks the rest of the questions.
##
## Usage: sbatch --mem=200G --cpus-per-task=2 --time=3:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_identify_clusters.R
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(RANN); setDTthreads(2)
})
OUT <- "output/19_celltypes"
HC  <- fread(file.path(OUT,"hc_cells.csv.gz"))
LB  <- fread(file.path(OUT,"hc_cells_labelled_high.csv.gz"), na.strings=c("NA",""))
MD  <- fread("output/16_regions/cells_meta_all.csv.gz",
             select=c("cell_id","samp","fov_id","nCount_RNA","nCount_negprobes",
                      "Area.um2","inSitu_clusters"))
NF  <- rbindlist(lapply(c("IB1","RB1"), function(S)
        fread(sprintf("output/18_leakage/%s_centroids.csv.gz",S),
              select=c("cell_id","cx_px","cy_px","nFeature_RNA","aspect","circularity"))))
TGT <- c(4L, 7L)
d <- NF[MD[HC[, .(cell_id, merge_sample, hc_clusters_low, hc_clusters_high, U1, U2)],
           on=.(cell_id)], on=.(cell_id)]
d <- LB[, .(cell_id, celltype_cluster, lineage_argmax)][d, on=.(cell_id)]
setnames(d, "hc_clusters_low", "cl")
cat(sprintf("cells %d | target clusters %s\n", nrow(d), paste(TGT, collapse=", ")))

cat("\n=== 1. QC profile vs the rest of the object ===\n")
d[, grp := fifelse(cl %in% TGT, paste0("cluster ", cl), "all other clusters")]
print(d[, .(cells=.N,
            nCount_p25 = as.numeric(quantile(nCount_RNA,.25)),
            nCount_med = as.numeric(median(nCount_RNA)),
            nFeature_med = as.numeric(median(nFeature_RNA)),
            negprobe_pooled = round(sum(nCount_negprobes)/sum(nCount_RNA), 6),
            area_med = round(median(Area.um2),1),
            aspect_med = round(median(aspect),2),
            circ_med = round(median(circularity),2),
            pct_IB1 = round(100*mean(merge_sample=="IB1"),1)), by=grp][order(grp)])
cat("\n-- how close are they to the per-cell floor (nCount 50 / nFeature 30)? --\n")
print(d[, .(pct_under_100ct = round(100*mean(nCount_RNA<100),2),
            pct_under_60ft  = round(100*mean(nFeature_RNA<60),2)), by=grp][order(grp)])

cat("\n=== 2. vendor InSituType composition (i = the over-fixation defect, excluded upstream) ===\n")
print(dcast(d[, .N, by=.(grp, inSitu_clusters)], grp ~ inSitu_clusters, value.var="N", fill=0L))

cat("\n=== 3. where do they sit in the FINER clustering? ===\n")
for (k in TGT) {
  x <- d[cl==k, .N, by=hc_clusters_high][order(-N)][1:5]
  x[, pct := round(100*N/d[cl==k,.N],1)]
  cat("\n-- low cluster", k, "-> hc_clusters_high:\n"); print(x)
  cat("   labels there: "); print(d[cl==k, .N, by=celltype_cluster][order(-N)][1:4])
}

cat("\n=== 4. spatial coherence: is this a place in the tissue, or scattered cells? ===\n")
## fraction of a cluster's cells whose nearest neighbour (same sample) is also in it,
## against the fraction expected if the cluster were sprinkled at random
d[, `:=`(X = cx_px*0.120280945, Y = -cy_px*0.120280945)]
for (k in c(TGT, 0L, 1L)) {
  fr <- sapply(c("IB1","RB1"), function(S) {
    z <- d[merge_sample==S]; if (nrow(z) < 1000) return(NA_real_)
    r <- RANN::nn2(z[, .(X,Y)], z[, .(X,Y)], k=2L)
    mean(z$cl[r$nn.idx[,2]][z$cl==k] == k)
  })
  exp <- sapply(c("IB1","RB1"), function(S) d[merge_sample==S, mean(cl==k)])
  cat(sprintf("  cluster %d: nn-same %.3f (IB1) / %.3f (RB1)  vs random %.3f / %.3f  -> %.1fx / %.1fx\n",
              k, fr[1], fr[2], exp[1], exp[2], fr[1]/exp[1], fr[2]/exp[2]))
}

cat("\n=== 5. FOV concentration (a QC artefact tends to sit in particular FOVs) ===\n")
for (k in TGT) {
  z <- d[cl==k, .N, by=.(merge_sample, fov_id)][order(-N)]
  tot <- d[cl==k, .N]
  cat(sprintf("  cluster %d: %d FOVs; top 10 hold %.1f%% of it; top FOV %.1f%%\n",
              k, nrow(z), 100*sum(z$N[1:min(10,nrow(z))])/tot, 100*z$N[1]/tot))
}
nf <- d[, .(n_fov = uniqueN(paste(merge_sample, fov_id))), by = grp]
cat("\n  FOVs occupied:\n"); print(nf)

cat("\n=== 6. what is expressed there? top genes by pooled cp10k vs the rest ===\n")
seu <- readRDS(file.path(OUT,"hc_02_clustered.Rds"))
ca  <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
setkey(d, cell_id); d <- d[J(colnames(ca))]
tot <- seu$nCount_RNA
for (k in TGT) {
  i <- which(d$cl == k); j <- which(!d$cl %in% TGT)
  a <- as.numeric(BPCells::rowSums(ca[, i])) / sum(tot[i]) * 1e4
  b <- as.numeric(BPCells::rowSums(ca[, j])) / sum(tot[j]) * 1e4
  g <- data.table(gene = rownames(ca), cp10k = round(a,2), rest = round(b,2),
                  fold = round((a+.01)/(b+.01), 2))
  cat("\n-- cluster ", k, ": top 25 enriched genes --\n", sep="")
  print(g[cp10k >= 1][order(-fold)][1:25])
}
cat("\ndone\n")
