## Diagnostic for step 10: is the cHL "POS" programme actually Hodgkin-specific, or is
## it a generic activated-lymphoid programme that happens to be well measured?
##
## WHY. scripts/10_script_hrs_call.R scores 21 POS genes and holds TNFRSF8 (CD30) out
## as the validator. CD30 turned out to be anemic (max 6 counts, 78% singletons) while
## the POS genes are robustly expressed -- so the validator cannot arbitrate, and the
## published conclusion rests on a test that may simply lack power.
##
## The score is z-scored WITHIN SAMPLE, which makes it mean-0 in each sample by
## construction and therefore structurally incapable of showing an IB1-vs-RB1
## difference in the programme. This script does the comparison the score cannot:
##   1. per-gene detection and expression on a COMMON scale, depth-standardised
##   2. IB1 vs RB1 for the programme as a whole, depth-standardised
##   3. which lineages actually express these genes -- if it is T cells and myeloid,
##      the programme is not HRS-specific no matter how well it is measured
##   4. how much each gene contributes to the score
##
##   sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 \
##     scripts/sbatch/run_rscript.sh scripts/sbatch/qc_10_chl_program.R
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); setDTthreads(1)
})

IN_RDS <- "output/04_merge/merged_processed.Rds"
SEG_F  <- "output/10_hrs/hrs_segments.csv.gz"
OUT    <- "output/10_hrs"; FIG <- file.path(OUT, "figures")
SAMPLES <- c("IB1","RB1"); DEFECT_CL <- "i"

POS <- c("CCL17","CCL22","IRF4","TRAF1","CD274","PDCD1LG2","FUT4","BIRC3","NFKB2",
         "RELB","CD40","TNFAIP3","GATA3","ID2","LGALS1","IL13","IL21","BATF3",
         "STAT6","JAK2","SOCS1")
CD30 <- "TNFRSF8"
SHOWN <- c("CCL17","IRF4","TRAF1")   # the three drawn in the step-10 views

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
cat("step 10 diagnostic -- is the cHL POS programme Hodgkin-specific?\n")
cat("============================================================\n")

obj <- readRDS(IN_RDS)
rn  <- rownames(obj@meta.data)
md  <- as.data.table(obj@meta.data[, c("merge_sample","nCount_RNA","inSitu_clusters",
                                       "lineage"), drop = FALSE])
md[, cell_id := rn]
cnt <- LayerData(obj[["RNA"]], layer = "counts")
G <- intersect(c(POS, CD30), rownames(cnt))
M <- as(cnt[G, , drop = FALSE], "dgCMatrix"); colnames(M) <- rn
rm(obj, cnt); invisible(gc())

md <- md[merge_sample %in% SAMPLES & nCount_RNA > 0 & inSitu_clusters != DEFECT_CL]
gi <- match(md$cell_id, colnames(M)); stopifnot(!anyNA(gi))
M  <- M[, gi, drop = FALSE]
md[, samp := factor(merge_sample, levels = SAMPLES)]
## depth deciles on the POOLED set, so the two samples are standardised to one scale
md[, dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                include.lowest = TRUE, labels = FALSE)]
cat("\ncells:\n"); print(md[, .(cells = .N, med_depth = as.numeric(median(nCount_RNA))), by = samp])

## ---- 1. per gene, common scale, depth-standardised --------------------------
## Counts per 10k rather than a within-sample z, which is exactly what the score uses
## and exactly what makes an IB1-vs-RB1 comparison impossible.
cp10k <- Matrix::t(Matrix::t(M) / md$nCount_RNA) * 1e4
wts <- md[, .(w = as.numeric(.N)), by = dec][, w := w/sum(w)]
per_gene <- rbindlist(lapply(G, function(g) {
  v <- as.vector(cp10k[g, ])
  d <- data.table(samp = md$samp, dec = md$dec, v = v, det = v > 0)
  s <- d[, .(cp = mean(v), pct = 100*mean(det), n = .N), by = .(samp, dec)]
  s <- merge(s, wts, by = "dec")
  o <- s[, .(cp10k_std = sum(cp*w)/sum(w), pct_det_std = sum(pct*w)/sum(w)), by = samp]
  data.table(gene = g, samp = o$samp, cp10k = round(o$cp10k_std, 3),
             pct_det = round(o$pct_det_std, 2))
}))
pw <- dcast(per_gene, gene ~ samp, value.var = c("cp10k","pct_det"))
pw[, ratio_IB1_RB1 := round(cp10k_IB1 / pmax(cp10k_RB1, 1e-9), 2)]
pw[, in_views := gene %in% SHOWN][, is_validator := gene == CD30]
setorder(pw, -cp10k_IB1)
cat("\n=== per gene, depth-standardised counts per 10k and detection %\n")
cat("    (a Hodgkin-specific gene should have ratio_IB1_RB1 well above 1)\n")
print(pw)
sc("chl_program_per_gene.csv", pw)

## ---- 2. the programme as a whole --------------------------------------------
prog <- as.vector(Matrix::colSums(cp10k[intersect(POS, G), , drop = FALSE]))
d <- data.table(samp = md$samp, dec = md$dec, prog = prog,
                cd30 = as.vector(cp10k[CD30, ]), lineage = as.character(md$lineage))
s <- merge(d[, .(v = mean(prog), c = mean(cd30), n = .N), by = .(samp, dec)], wts, by = "dec")
tot <- s[, .(prog_cp10k_std = round(sum(v*w)/sum(w), 2),
             cd30_cp10k_std = round(sum(c*w)/sum(w), 4), cells = sum(n)), by = samp]
cat("\n=== whole POS programme, depth-standardised, common scale\n"); print(tot)
cat(sprintf("IB1:RB1 programme ratio %.2f | CD30 ratio %.2f\n",
            tot[samp=="IB1", prog_cp10k_std] / tot[samp=="RB1", prog_cp10k_std],
            tot[samp=="IB1", cd30_cp10k_std] / max(tot[samp=="RB1", cd30_cp10k_std], 1e-9)))
sc("chl_program_sample_compare.csv", tot)

## ---- 3. which lineages carry it ---------------------------------------------
lin <- merge(d[, .(prog = mean(prog), n = .N), by = .(samp, dec, lineage)], wts, by = "dec")
lin <- lin[, .(prog_cp10k_std = round(sum(prog*w)/sum(w), 2), cells = sum(n)),
           by = .(samp, lineage)]
cat("\n=== who expresses the POS programme (depth-standardised)\n")
cat("    If T/myeloid carry it, the programme is not HRS-specific however well measured\n")
print(dcast(lin, lineage ~ samp, value.var = "prog_cp10k_std")[order(-IB1)])
sc("chl_program_by_lineage.csv", lin)

## ---- 4. what actually drives the score --------------------------------------
seg <- fread(SEG_F, select = c("cell_id","samp","hrs_score","z_pos"))
sj  <- merge(md[, .(cell_id, samp2 = samp)], seg, by = "cell_id")
gm  <- match(sj$cell_id, md$cell_id)
contrib <- rbindlist(lapply(intersect(POS, G), function(g) {
  v <- as.vector(cp10k[g, gm])
  data.table(gene = g,
             cor_with_score = round(suppressWarnings(cor(v, sj$hrs_score,
                                     method = "spearman", use = "complete.obs")), 3),
             pct_det = round(100*mean(v > 0), 2))
}))[order(-cor_with_score)]
cat("\n=== each POS gene's rank correlation with the final score (IB1+RB1)\n")
print(contrib)
sc("chl_program_score_contribution.csv", contrib)

## ---- figures -----------------------------------------------------------------
pl <- melt(pw[, .(gene, IB1 = cp10k_IB1, RB1 = cp10k_RB1, in_views, is_validator)],
           id.vars = c("gene","in_views","is_validator"),
           variable.name = "samp", value.name = "cp10k")
p1 <- ggplot(pl, aes(reorder(gene, cp10k), cp10k, fill = samp)) +
  geom_col(position = position_dodge(.8), width = .74) + coord_flip() +
  scale_y_log10() + scale_fill_manual(values = SAMP_COLS, name = NULL) +
  labs(title = "The scored cHL genes are well measured; the validator is not",
       subtitle = paste0("Depth-standardised counts per 10k, common scale. TNFRSF8 (CD30) ",
                         "is the held-out validator\nand sits orders of magnitude below the ",
                         "genes that actually drive the score."),
       x = NULL, y = "counts per 10k (log scale, depth-standardised)",
       caption = sprintf("Data: %s/chl_program_per_gene.csv", OUT)) + theme_viz(xgrid = TRUE)
sf("fig_10_chl_program_per_gene.png", p1, 8.5, 8)

p2 <- ggplot(lin, aes(reorder(lineage, prog_cp10k_std), prog_cp10k_std, fill = samp)) +
  geom_col(position = position_dodge(.8), width = .74) + coord_flip() +
  scale_fill_manual(values = SAMP_COLS, name = NULL) +
  labs(title = "Which cells express the scored cHL programme",
       subtitle = paste0("Depth-standardised. A programme carried by T and myeloid cells in ",
                         "BOTH samples is a generic\nactivated-lymphoid signature, not a ",
                         "Hodgkin-specific one."),
       x = NULL, y = "POS programme, counts per 10k",
       caption = sprintf("Data: %s/chl_program_by_lineage.csv", OUT)) + theme_viz(xgrid = TRUE)
sf("fig_10_chl_program_by_lineage.png", p2, 8, 5.5)

cat("\n============================================================\n")
cat("Read chl_program_per_gene.csv and chl_program_by_lineage.csv together.\n")
cat("If the POS genes are expressed at similar levels in IB1 and RB1, and largely by\n")
cat("T/myeloid cells, then the score is measuring a generic activated-lymphoid\n")
cat("programme well -- not Hodgkin biology -- and no threshold on it can find HRS.\n")
cat("============================================================\nOK\n")
