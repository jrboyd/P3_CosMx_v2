## Is the IB1 CCL17/CCL22 elevation Hodgkin biology, or a CELLULAR COMPOSITION artifact?
##
## THE PROBLEM. Step 11 found CCL17+CCL22 elevated 28.9x over the RB1 null and called
## 29% of IB1 bins "niche" -- which is a quarter of the tissue, not a focal niche. The
## composition table shows why that needs checking:
##       IB1 background   f_B  6.8%   f_myeloid 17.1%
##       RB1 all          f_B 38.9%   f_myeloid  3.3%
## RB1's node is effaced by follicular lymphoma and is overwhelmingly B cells. IB1
## retains far more non-B tissue. CCL17/TARC and CCL22/MDC are made by dendritic cells
## and are hallmarks of the interfollicular T/DC zone, so "IB1 has more CCL17 than RB1"
## may simply mean "IB1 has more T-zone tissue left, because RB1's is replaced by FL".
## That is a composition confound, not Hodgkin biology.
##
## THE TEST. If the elevation is Hodgkin-specific it should survive within lineage --
## an IB1 myeloid cell should carry more CCL17 than an RB1 myeloid cell, and so on. If
## instead each lineage carries the same amount in both samples and only their
## PROPORTIONS differ, the whole-sample difference is composition and the "niche" is
## just where the T/DC zone is.
##
##   sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 \
##     scripts/sbatch/run_rscript.sh scripts/sbatch/qc_11_niche_composition_confound.R
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); setDTthreads(1)
})

IN_RDS <- "output/04_merge/merged_processed.Rds"
OUT    <- "output/11_hodgkin_niche"; FIG <- file.path(OUT, "figures")
SAMPLES <- c("IB1","RB1"); DEFECT_CL <- "i"
CHEMO  <- c("CCL17","CCL22"); EXTRA <- c("BATF3","TNFRSF8","CD3E","CD68","MS4A1")

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
SAMP_COLS <- c(IB1 = "#2a78d6", RB1 = "#eb6834")
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
          legend.text=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left")
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF)
                       cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

cat("============================================================\n")
cat("step 11 QC -- is the CCL17/CCL22 elevation composition, or Hodgkin?\n")
cat("============================================================\n")

obj <- readRDS(IN_RDS); rn <- rownames(obj@meta.data)
md <- as.data.table(obj@meta.data[, c("merge_sample","nCount_RNA","inSitu_clusters",
                                      "lineage"), drop = FALSE])
md[, cell_id := rn]
cnt <- LayerData(obj[["RNA"]], layer = "counts")
G <- intersect(c(CHEMO, EXTRA), rownames(cnt))
M <- as(cnt[G, , drop = FALSE], "dgCMatrix"); colnames(M) <- rn
rm(obj, cnt); invisible(gc())

md <- md[merge_sample %in% SAMPLES & nCount_RNA > 0 & inSitu_clusters != DEFECT_CL]
gi <- match(md$cell_id, colnames(M)); M <- M[, gi, drop = FALSE]
md[, samp := factor(merge_sample, levels = SAMPLES)]
md[, lin  := as.character(lineage)]
md[, chemo := as.vector(Matrix::colSums(M[intersect(CHEMO, G), , drop = FALSE]))]
md[, chemo_cp10k := 1e4 * chemo / nCount_RNA]
md[, dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                include.lowest = TRUE, labels = FALSE)]

## ---- 1. composition ----------------------------------------------------------
comp <- dcast(md[, .N, by = .(samp, lin)], lin ~ samp, value.var = "N", fill = 0)
comp[, `:=`(pct_IB1 = round(100*IB1/sum(IB1), 1), pct_RB1 = round(100*RB1/sum(RB1), 1))]
cat("\n=== cellular composition differs enormously between the samples\n")
print(comp[order(-pct_IB1), .(lin, pct_IB1, pct_RB1)])
sc("niche_confound_composition.csv", comp)

## ---- 2. within-lineage, depth-standardised -----------------------------------
## The decisive table. Depth deciles are pooled so both samples sit on one scale.
w <- md[, .(w = as.numeric(.N)), by = dec][, w := w/sum(w)]
per <- merge(md[, .(v = mean(chemo_cp10k), n = .N), by = .(samp, lin, dec)], w, by = "dec")
per <- per[, .(chemo_cp10k = sum(v*w)/sum(w), cells = sum(n)), by = .(samp, lin)]
pw  <- dcast(per, lin ~ samp, value.var = c("chemo_cp10k","cells"))
pw[, ratio := round(chemo_cp10k_IB1 / pmax(chemo_cp10k_RB1, 1e-9), 2)]
pw[, `:=`(chemo_cp10k_IB1 = round(chemo_cp10k_IB1, 2),
          chemo_cp10k_RB1 = round(chemo_cp10k_RB1, 2))]
setorder(pw, -chemo_cp10k_IB1)
cat("\n=== CCL17+CCL22 per cell, WITHIN lineage, depth-standardised\n")
cat("    Hodgkin biology -> ratio well above 1 in most lineages.\n")
cat("    Composition confound -> ratios near 1, and only the proportions differ.\n")
print(pw)
sc("niche_confound_within_lineage.csv", pw)

## ---- 3. reweight IB1 to RB1's composition, and vice versa --------------------
## Direct standardisation on lineage. If the whole-sample gap collapses once the two
## samples are given the SAME cell mix, the gap was composition.
tot_raw <- merge(md[, .(v = mean(chemo_cp10k), n = .N), by = .(samp, dec)], w, by = "dec")[
             , .(raw = sum(v*w)/sum(w)), by = samp]
mix <- md[, .N, by = .(samp, lin)][, p := N/sum(N), by = samp]
std_to <- function(target) {
  tp <- mix[samp == target, .(lin, p)]
  merge(per, tp, by = "lin")[, .(std = sum(chemo_cp10k * p)/sum(p)), by = samp]
}
s_ib <- std_to("IB1")[, .(samp, std_to_IB1_mix = round(std, 2))]
s_rb <- std_to("RB1")[, .(samp, std_to_RB1_mix = round(std, 2))]
res <- merge(merge(tot_raw[, .(samp, raw = round(raw, 2))], s_ib, by = "samp"), s_rb, by = "samp")
cat("\n=== whole-sample CCL17+CCL22, raw vs composition-standardised\n")
print(res)
r_raw <- res[samp=="IB1", raw] / res[samp=="RB1", raw]
r_ib  <- res[samp=="IB1", std_to_IB1_mix] / res[samp=="RB1", std_to_IB1_mix]
r_rb  <- res[samp=="IB1", std_to_RB1_mix] / res[samp=="RB1", std_to_RB1_mix]
cat(sprintf("\nIB1:RB1 ratio  raw %.2fx | standardised to IB1's mix %.2fx | to RB1's mix %.2fx\n",
            r_raw, r_ib, r_rb))
sc("niche_confound_standardised.csv",
   data.table(ratio_raw = round(r_raw,2), ratio_std_IB1mix = round(r_ib,2),
              ratio_std_RB1mix = round(r_rb,2)))

## ---- verdict -------------------------------------------------------------------
surv <- mean(c(r_ib, r_rb)) / r_raw
cat("\n============================================================\n")
cat(sprintf("Composition-standardised ratio retains %.0f%% of the raw ratio.\n", 100*surv))
## braces on every arm: a bare `else` on its own line after a closing paren is a
## parse error at top level in R
r_std <- mean(c(r_ib, r_rb))
if (r_std >= 2) {
  cat("VERDICT: the elevation SURVIVES composition standardisation -- it is not merely\n",
      "         a consequence of RB1's node being effaced by FL. Hodgkin-associated\n",
      "         chemokine signalling in IB1 is supported.\n", sep = "")
} else if (r_std >= 1.3) {
  cat("VERDICT: PARTIAL. Some elevation survives, but a large share of the raw gap was\n",
      "         cellular composition. Report the standardised number, not the raw one.\n", sep="")
} else {
  cat("VERDICT: the elevation is largely COMPOSITION. IB1 retains T/DC-zone tissue that\n",
      "         RB1 has lost to FL, and CCL17/CCL22 track that zone. The step-11 'niche'\n",
      "         is then the interfollicular zone, not a Hodgkin niche -- relabel it.\n", sep="")
}

pl <- melt(pw[, .(lin, IB1 = chemo_cp10k_IB1, RB1 = chemo_cp10k_RB1)], id.vars = "lin",
           variable.name = "samp", value.name = "cp10k")
p1 <- ggplot(pl, aes(reorder(lin, cp10k), cp10k, fill = samp)) +
  geom_col(position = position_dodge(.8), width = .74) + coord_flip() +
  scale_fill_manual(values = SAMP_COLS, name = NULL) +
  labs(title = "CCL17+CCL22 per cell, within each lineage",
       subtitle = paste0("Depth-standardised, common scale. If the bars match within ",
                         "lineage, the whole-sample\ndifference is cellular composition ",
                         "rather than Hodgkin biology."),
       x = NULL, y = "CCL17+CCL22, counts per 10k",
       caption = sprintf("Data: %s/niche_confound_within_lineage.csv", OUT)) +
  theme_viz(xgrid = TRUE)
sf("fig_11_confound_within_lineage.png", p1, 8, 5.5)
cat("============================================================\nOK\n")
