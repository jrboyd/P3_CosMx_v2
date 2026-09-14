## ============================================================================
## 19_script_t_substructure.R -- split the T-cell label. Sub-clustering, then composition.
##
## T is ONE high-confidence cluster of 142,361 cells, and steps 16 and 17 are built on Tfh.
## B was split by a module rule and pDC and the Hodgkin niche were named by module rules,
## because each of those was already a SEPARATE cluster that the coarse vocabulary could not
## name. T is the opposite problem: one cluster containing several populations. A module rule
## cannot help -- there is nothing to score against. So the T cells are re-clustered among
## themselves and the sub-clusters are then labelled by argmax composition, which is the same
## logic 19d uses at the top level.
##
## THE ARGMAX T SUBSETS ALREADY EXIST and are reused, not recomputed: step 16 ran a second
## hierarchical stage inside T (Tfh 25,025, Treg 28,326, CD8_T 22,991, CD4_conv 54,827,
## T_unresolved 183,583 over the whole sample). 59% unresolved is exactly the gap
## sub-clustering is meant to close.
##
## THE CHECKPOINT GENES ARE THE VALIDATION, AND THEY ARE CLEAN. Step 16 deliberately held
## PDCD1 out of the Tfh set and CTLA4 out of the Treg set, so that "Tfh express PD-1" could
## not become a restatement of the definition. That makes them perfect independent tests
## here: if the sub-cluster labelled Tfh tops PDCD1 and the one labelled Treg tops CTLA4,
## neither gene having been used to define anything, the split is real.
##
## Usage: sbatch --mem=200G --cpus-per-task=8 --time=6:00:00 \
##          scripts/sbatch/19_hc_object_wrapper.sh   # (thread pinning matters for UMAP)
##   or:  sbatch --mem=200G --cpus-per-task=8 --time=6:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_t_substructure.R
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(scPearsonPCA)
  library(ggplot2); library(patchwork); setDTthreads(2)
})
source("scripts/lib/lineage.R")
OUT <- "output/19_celltypes"; FIG <- file.path(OUT,"figures")
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=160,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}
T_SETS <- list(Tfh      = c("CXCL13","ICOS","CXCR5","BTLA","CD40LG"),
               Treg     = c("FOXP3","IL2RA","IKZF2"),
               CD8_T    = c("CD8A","GZMK","GZMB"),
               CD4_conv = c("CD4","IL7R","CCR7"))
HELD_OUT <- c("PDCD1","CTLA4","TIGIT","LAG3","TOX","CD8B","FOXP3")  # PDCD1/CTLA4 never used
CONF_PUR <- 0.40; CONF_MRG <- 0.15; CONF_ENR <- 2.0; WEAK_ENR <- 1.5; WEAK_PUR <- 0.25

## PER-SAMPLE MODE. The joint sub-clustering separated SAMPLE, not T subset: against a 26.2%
## IB1 baseline the sub-clusters ran 4.6% to 82.7% IB1, while every one of them had CD4_conv
## as its modal argmax subset at enrichment 0.75-1.26 -- no enrichment at all. IB1 is a
## composite cHL node and RB1 is an FL node, so their T cells genuinely differ, and that
## between-sample variance dominates the within-T variance the split needs. Running each
## sample on its own removes it. If subsets still fail to separate here, they are not
## separable by clustering on this panel and the per-cell argmax call is the only route.
args <- commandArgs(trailingOnly = TRUE)
SAMP <- if (length(args) >= 1) args[1] else "both"
SUB_RDS <- file.path(OUT, sprintf("t_subclustered_%s.Rds", SAMP))
lab <- fread(file.path(OUT,"hc_cells_final.csv.gz"), na.strings=c("NA",""))
if (SAMP != "both") lab <- lab[merge_sample == SAMP]
cat("sample scope: ", SAMP, "\n", sep="")
ct  <- fread("output/16_celltypes/cell_celltypes.csv.gz", select=c("cell_id","celltype","lineage16"))
tc  <- lab[celltype_cluster == "T cell", cell_id]
cat(sprintf("T-labelled high-confidence cells: %d\n", length(tc)))

if (!file.exists(SUB_RDS)) {
  cat("\nreading the high-confidence object\n")
  seu <- readRDS(file.path(OUT,"hc_02_clustered.Rds"))
  seu <- subset(seu, cells = tc)
  cat(sprintf("subset: %d cells\n", ncol(seu)))
  ca <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
  tcnt <- Matrix::colSums(ca); gf <- scPearsonPCA::gene_frequency(ca)
  ## HVGs are recomputed WITHIN the T cells. The object-level 900 were chosen to separate
  ## lineages, so they carry little of the variation that distinguishes one T cell from
  ## another; reusing them would sub-cluster on the wrong axes.
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 900)
  hv <- VariableFeatures(seu)
  ch <- as(ca[hv, ], "dgCMatrix")
  cat(sprintf("  %d HVGs x %d cells, %.2f%% nonzero\n", nrow(ch), ncol(ch),
              100*length(ch@x)/(as.numeric(nrow(ch))*ncol(ch))))
  po <- sparse_quasipoisson_pca_seurat(ch, totalcounts = tcnt, grate = gf[hv],
                                       scale.max = 10, do.scale = TRUE, do.center = TRUE)
  rm(ch); invisible(gc())
  uo <- scPearsonPCA::make_umap(po)
  seu[["tpca"]]   <- po$reduction.data
  seu[["tumap"]]  <- uo$ump
  seu[["tgraph"]] <- Seurat::as.Graph(uo$grph)
  seu <- Seurat::FindClusters(seu, graph = "tgraph", cluster.name = "t_sub", resolution = .6)
  seu$seurat_clusters <- NULL
  saveRDS(seu, SUB_RDS)
} else { cat("\nreusing ", SUB_RDS, "\n", sep=""); seu <- readRDS(SUB_RDS) }

emb <- as.data.table(Embeddings(seu, "tumap"), keep.rownames="cell_id")
setnames(emb, 2:3, c("T1","T2"))
d <- data.table(cell_id = colnames(seu), t_sub = as.character(seu$t_sub),
                samp = as.character(seu$merge_sample))
d <- ct[emb[d, on=.(cell_id)], on=.(cell_id)]
cat("\n=== sub-clusters ===\n"); print(d[, .N, by=t_sub][order(-N)])

## ---- label the sub-clusters by argmax T-subset composition ------------------
SUBS <- names(T_SETS)
d[, resolved := celltype %in% SUBS]
base <- d[resolved == TRUE, .N, by=celltype][, .(celltype, base=N/sum(N))]
cat("\n=== argmax T-subset baseline among resolved T cells ===\n"); print(base[order(-base)])
comp <- d[resolved == TRUE, .N, by=.(t_sub, celltype)][, share := N/sum(N), by=t_sub]
r <- comp[order(-share), .(top1=celltype[1], purity=round(share[1],3),
                           top2=celltype[2], margin=round(share[1]-share[2],3)), by=t_sub]
tot <- d[, .(cells=.N, frac_resolved=round(mean(resolved),3),
             pct_IB1=round(100*mean(samp=="IB1"),1)), by=t_sub]
r <- base[r, on=c(celltype="top1")]; setnames(r, "celltype", "top1")
r <- tot[r, on=.(t_sub)]
r[, enrichment := round(purity/base, 2)]
r[, status := fifelse((purity>=CONF_PUR & margin>=CONF_MRG) | (enrichment>=CONF_ENR & purity>=WEAK_PUR),
              "confident", fifelse(enrichment>=WEAK_ENR & purity>=WEAK_PUR, "weak_support", "mixed"))]
r[, label := fifelse(status=="mixed", NA_character_, paste0("T, ", top1))]
setorder(r, -cells)
cat("\n=== sub-cluster labels ===\n"); print(r)
sc(sprintf("t_subcluster_labels_%s.csv", SAMP), r)

## ---- the held-out test -----------------------------------------------------
ca <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
G <- intersect(c(unlist(T_SETS), HELD_OUT), rownames(ca))
sub <- as(as.matrix(ca[G,,drop=FALSE]),"dgCMatrix"); rownames(sub) <- G
setkey(d, cell_id); d <- d[J(colnames(ca))]
ncnt <- seu$nCount_RNA
prof <- rbindlist(lapply(unique(d$t_sub), function(k) {
  j <- which(d$t_sub == k)
  data.table(t_sub=k, gene=G,
             cp10k=as.numeric(Matrix::rowSums(sub[,j,drop=FALSE]))/sum(ncnt[j])*1e4)
}))
prof <- r[, .(t_sub, label)][prof, on=.(t_sub)]
sc(sprintf("t_subcluster_markers_%s.csv", SAMP), prof)
cat("\n=== HELD-OUT validation: PDCD1 and CTLA4 were never used to define anything ===\n")
for (g in c("PDCD1","CTLA4","CD8B","FOXP3")) {
  x <- prof[gene==g & !is.na(label)][order(-cp10k)]
  if (!nrow(x)) next
  cat(sprintf("  %-6s top: %-16s (%.2f)   2nd: %-16s (%.2f)\n", g,
              x$label[1], x$cp10k[1], x$label[2], x$cp10k[2]))
}

## ---- figure -----------------------------------------------------------------
d <- r[, .(t_sub, label, status)][d, on=.(t_sub)]
d[, lab := fifelse(is.na(label), paste0("<", status, ">"), label)]
PAL <- c(`T, Tfh`="#c51b7d", `T, Treg`="#0b6b4f", `T, CD8_T`="#d03b3b",
         `T, CD4_conv`="#184f95", `<mixed>`="#b9b8b2")
set.seed(1L); s <- d[sample(.N, min(.N, 120000L))]
lb <- s[!is.na(label), .(T1=median(T1), T2=median(T2)), by=lab]
p1 <- ggplot(s, aes(T1,T2,colour=lab)) + geom_point(size=.12, alpha=.5) +
  scale_colour_manual(values=PAL, name=NULL) + coord_fixed() +
  guides(colour=guide_legend(override.aes=list(size=3,alpha=1))) +
  ggrepel::geom_text_repel(data=lb, aes(label=lab), colour="#0b0b0b", size=3.2,
                           fontface="bold", show.legend=FALSE, seed=1, bg.color="#fcfcfb") +
  labs(title=sprintf("T-cell sub-clustering -- %d cells", nrow(d)),
       subtitle="labels from argmax T-subset composition; step 16 left 59% of T cells unresolved",
       x="T-UMAP 1", y="T-UMAP 2")
hm <- prof[!is.na(label)]
hm[, z := (cp10k-mean(cp10k))/pmax(sd(cp10k),1e-9), by=gene]
hm[, isheld := gene %in% HELD_OUT]
p2 <- ggplot(hm, aes(label, gene, fill=z)) + geom_tile(colour="#fcfcfb", linewidth=.3) +
  scale_fill_gradient2(low="#184f95", mid="#f0efec", high="#d03b3b", midpoint=0,
                       limits=c(-1.5,1.5), oob=scales::squish, name="z across sub-clusters") +
  facet_grid(isheld ~ ., scales="free_y", space="free_y", switch="y",
             labeller=as_labeller(c(`TRUE`="HELD OUT", `FALSE`="defining"))) +
  labs(title="marker evidence", subtitle="PDCD1 and CTLA4 were held out of the definitions on purpose",
       x=NULL, y=NULL)
th <- theme_minimal(11) + theme(plot.background=element_rect(fill="#fcfcfb",colour=NA),
  panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
  plot.subtitle=element_text(size=8), axis.text.x=element_text(angle=35,hjust=1),
  strip.text=element_text(face="bold",size=8), legend.position="bottom")
sf(sprintf("t_substructure_%s.png", SAMP), (p1 + th | p2 + th) + plot_layout(widths=c(1.2,1)), 15, 7)

cat("\n=== per sample ===\n")
print(dcast(d[, .N, by=.(lab, samp)][, pct := round(100*N/sum(N),1), by=samp],
            lab ~ samp, value.var="pct"))
sc(sprintf("t_cells_labelled_%s.csv.gz", SAMP), d[, .(cell_id, samp, t_sub, celltype_t = label,
                                    t_status = status, celltype_argmax_t = celltype, T1, T2)])
cat("\ndone\n")
