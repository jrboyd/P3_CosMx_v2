## ============================================================================
## qc_19_t_umap_subsets.R -- EXPLORATORY, not a numbered step.
##
## Where do step 16's argmax T subsets land on step 19's T-cell UMAP?
##
## NOTHING IS RECOMPUTED. Step 19 already did the subset-and-reprocess:
## 19_script_t_substructure.R takes the T-labelled cells out of the
## high-confidence object, recomputes 900 HVGs WITHIN T, runs a quasipoisson
## PCA, a UMAP and Louvain, and dumps the embedding to
## output/19_celltypes/t_cells_labelled_{IB1,RB1}.csv.gz as T1/T2, with the
## step-16 argmax T-subset call already joined on as celltype_argmax_t.
## This script only plots that table.
##
## PER-SAMPLE objects, not the joint one. The joint run separated SAMPLE and
## not subset (CLAUDE.md, closed avenues); its UMAP correlates 0.15/0.08 with
## the per-sample ones, so the per-sample embeddings are the honest view.
##
## Usage: Rscript --vanilla scripts/sbatch/qc_19_t_umap_subsets.R    (light, no sbatch)
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); setDTthreads(1)
})
OUT <- "output/19_celltypes"; FIG <- file.path(OUT, "figures")

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#d03b3b"; DIV_MID<-"#f0efec"; DIV_LO<-"#184f95"
theme_viz <- function(base=11, xgrid=FALSE) {
  th <- theme_minimal(base_size=base) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.major.y=element_line(colour=GRID,linewidth=.25),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.line.x=element_line(colour=AXIS,linewidth=.3), axis.ticks=element_blank(),
          axis.text=element_text(colour=MUTED,size=base-2),
          axis.title=element_text(colour=INK2,size=base-1),
          plot.title=element_text(colour=INK,face="bold",size=base+2),
          plot.subtitle=element_text(colour=INK2,size=base-1),
          plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
          strip.text=element_text(colour=INK,face="bold",size=base-1),
          legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left", legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF);cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

SUBS   <- c("Tfh","Treg","CD4_conv","CD8_T")
NOCALL <- c("T_unresolved","ambiguous","all_markers_0","below_background")
SUB_COLS <- c(Tfh="#c2185b", Treg="#1b7f5c", CD4_conv="#e08a1e", CD8_T="#1f4fa8")
## sample hues are fixed project-wide
SAMP_COLS <- c(IB1="#2a78d6", RB1="#eb6834")

read_one <- function(s) {
  d <- fread(file.path(OUT, sprintf("t_cells_labelled_%s.csv.gz", s)))
  d[, samp := s][, t_sub := factor(t_sub)]
  ## three buckets: a resolved step-16 T subset, a step-16 non-call, or a
  ## step-16 call to a DIFFERENT lineage (caller disagreement, 70.2% agreement fact)
  d[, grp := fifelse(celltype_argmax_t %in% SUBS, celltype_argmax_t,
             fifelse(celltype_argmax_t %in% NOCALL, "no call", "other lineage"))]
  d[]
}
dl <- lapply(c("IB1","RB1"), read_one); names(dl) <- c("IB1","RB1")

for (s in names(dl)) {
  d <- dl[[s]]
  cat(sprintf("\n=== %s: %d T cells, %d sub-clusters ===\n", s, nrow(d), uniqueN(d$t_sub)))
  print(d[, .N, by=grp][order(-N)])
}

## ---- panel A: the UMAP itself, by sub-cluster ------------------------------
p_sub <- function(s) {
  d <- dl[[s]][sample(.N)]
  ggplot(d, aes(T1, T2, colour=t_sub)) +
    geom_point(size=.18, alpha=.5, stroke=0) +
    scale_colour_brewer(palette="Dark2", name="sub-cluster") +
    coord_fixed() + theme_viz(10) +
    guides(colour=guide_legend(override.aes=list(size=2.5, alpha=1), nrow=1)) +
    labs(title=sprintf("%s -- step 19 T-cell UMAP", s),
         subtitle=sprintf("%s T cells, %d Louvain sub-clusters (res 0.6)", format(nrow(d), big.mark=","), uniqueN(d$t_sub)),
         x="T-UMAP 1", y="T-UMAP 2")
}

## ---- panel B: where each step-16 argmax subset lands -----------------------
## grey = every other T cell in the same sample, so each facet is a highlight
## against the full embedding rather than a re-scaled subset.
p_hi <- function(s) {
  d  <- dl[[s]]
  lv <- c(SUBS, "no call", "other lineage")
  hi <- rbindlist(lapply(lv, function(g) d[grp == g][, facet := g]))
  hi[, facet := factor(facet, levels=lv)]
  bg <- rbindlist(lapply(lv, function(g) d[, .(T1, T2)][, facet := g]))
  bg[, facet := factor(facet, levels=lv)]
  n  <- hi[, .N, by=facet]
  lab <- n[, .(facet, txt=sprintf("n = %s", format(N, big.mark=",")))]
  cols <- c(SUB_COLS, "no call"="#8d8b85", "other lineage"="#5d3fa8")
  ggplot(hi, aes(T1, T2)) +
    geom_point(data=bg, colour="#e6e5df", size=.10, alpha=.5, stroke=0) +
    geom_point(aes(colour=facet), size=.45, alpha=.9, stroke=0) +
    geom_text(data=lab, aes(x=Inf, y=-Inf, label=txt), hjust=1.1, vjust=-.7,
              size=3, colour=INK2, inherit.aes=FALSE) +
    scale_colour_manual(values=cols, guide="none") +
    facet_wrap(~facet, nrow=2) + coord_fixed() + theme_viz(10) +
    labs(subtitle=sprintf("%s -- step-16 argmax T-subset call, highlighted on the same embedding", s),
         x="T-UMAP 1", y="T-UMAP 2")
}

## ---- enrichment: is any sub-cluster actually a subset? --------------------
## share of a sub-cluster's RESOLVED cells that are subset k, over that subset's
## share of all resolved cells in the sample. log2, so 0 = at baseline.
enr <- rbindlist(lapply(names(dl), function(s) {
  d <- dl[[s]][celltype_argmax_t %in% SUBS]
  base <- d[, .N, by=celltype_argmax_t][, .(celltype_argmax_t, base=N/sum(N))]
  e <- d[, .N, by=.(t_sub, celltype_argmax_t)][, share := N/sum(N), by=t_sub]
  e <- base[e, on=.(celltype_argmax_t)]
  e[, `:=`(samp=s, enr=share/base, l2=log2(share/base))]
  ## per-sub-cluster resolved n, for reading the noise
  e[, n_resolved := sum(N), by=t_sub]
  e[]
}))
setorder(enr, samp, t_sub, -enr)
enr[, celltype_argmax_t := factor(celltype_argmax_t, levels=SUBS)]
cat("\n=== argmax-subset enrichment per sub-cluster ===\n")
print(enr[, .(samp, t_sub, subset=celltype_argmax_t, n=N, n_resolved,
              share=round(share,3), base=round(base,3), enr=round(enr,2))])
sc("qc_19_t_umap_subset_enrichment.csv", enr)

LIM <- max(abs(enr$l2))
p_enr <- ggplot(enr, aes(celltype_argmax_t, t_sub, fill=l2)) +
  geom_tile(colour=SURF, linewidth=.8) +
  geom_text(aes(label=sprintf("%.2f", enr)), size=3, colour=INK) +
  scale_fill_gradient2(low=DIV_LO, mid=DIV_MID, high=DIV_HI, midpoint=0,
                       limits=c(-LIM, LIM), name="log2 enrichment vs sample baseline") +
  facet_wrap(~samp, scales="free_y") + theme_viz() +
  theme(panel.grid.major.y=element_blank(), axis.line.x=element_blank()) +
  labs(subtitle="fold enrichment of each argmax subset within each sub-cluster (1.00 = at baseline)",
       x="step-16 argmax T subset", y="step 19 sub-cluster")

## ---- what IS the "other lineage" mass? ------------------------------------
## In RB1 these concentrate in one corner of the embedding, so name them.
oth <- rbindlist(lapply(names(dl), function(s) {
  d <- dl[[s]][grp == "other lineage"]
  d[, .N, by=.(samp, t_sub, celltype_argmax_t)][, share := N/sum(N), by=.(samp, t_sub)]
}))
setorder(oth, samp, t_sub, -N)
cat("\n=== step-16 calls behind the \"other lineage\" bucket, top 3 per sub-cluster ===\n")
print(oth[, head(.SD, 3), by=.(samp, t_sub)])
sc("qc_19_t_umap_other_lineage.csv", oth)

## per sub-cluster: fraction of ALL its cells that step 16 calls a non-T lineage
frac <- rbindlist(lapply(names(dl), function(s) dl[[s]][, .(
  cells=.N, pct_other=round(100*mean(grp=="other lineage"),1),
  pct_nocall=round(100*mean(grp=="no call"),1),
  pct_Tsubset=round(100*mean(grp %in% SUBS),1)), by=.(samp, t_sub)][order(t_sub)]))
cat("\n=== per sub-cluster: how step 16 sees it ===\n"); print(frac)
sc("qc_19_t_umap_subcluster_calls.csv", frac)

p_oth <- ggplot(frac, aes(t_sub, pct_other, fill=samp)) +
  geom_col(width=.7) + geom_text(aes(label=sprintf("%.0f%%", pct_other)), vjust=-.4, size=3, colour=INK2) +
  scale_fill_manual(values=SAMP_COLS, guide="none") +
  facet_wrap(~samp, scales="free_x") + theme_viz(10, xgrid=FALSE) +
  labs(subtitle="% of each sub-cluster that step 16 calls a NON-T lineage",
       x="sub-cluster", y="% non-T argmax call")

## ---- assemble --------------------------------------------------------------
for (s in names(dl)) {
  p <- (p_sub(s) | p_hi(s)) + plot_layout(widths=c(1, 1.75))
  sf(sprintf("qc_19_t_umap_subsets_%s.png", s), p, 16, 7)
}
sf("qc_19_t_umap_subset_enrichment.png", (p_enr / p_oth) + plot_layout(heights=c(1.3,1)), 11, 8)

## ---- and the frank read of it ---------------------------------------------
cat("\n=== resolved fraction: how much of the UMAP the step-16 subsets can even colour ===\n")
print(rbindlist(lapply(names(dl), function(s) dl[[s]][, .(samp=s, cells=.N,
      pct_resolved_subset=round(100*mean(celltype_argmax_t %in% SUBS),1),
      pct_no_call=round(100*mean(celltype_argmax_t %in% NOCALL),1),
      pct_other_lineage=round(100*mean(!celltype_argmax_t %in% c(SUBS,NOCALL)),1))])))
cat("\ndone\n")
