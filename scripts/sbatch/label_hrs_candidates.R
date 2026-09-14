# Label a likely HRS (Hodgkin Reed-Sternberg) set in IB1 BEFORE integration, so
# integration can be judged on whether it preserves them.
#
# HRS cells are the hardest preservation test available: genuinely rare (0.1-2%
# of cells in cHL tissue), biologically the point of the experiment, and they sit
# in a T-cell/myeloid-rich microenvironment that integration could easily smear
# them into. RB1 is the built-in SPECIFICITY control -- if the definition flags a
# comparable fraction of RB1, it is not identifying HRS.
#
# Phenotype encoded here (classic Hodgkin lymphoma):
#   POSITIVE  CD30/TNFRSF8 (defining), TARC/CCL17 (HRS-secreted), MUM1/IRF4,
#             NF-kB program (TRAF1, BIRC3, NFKB2, RELB, CD40, TNFAIP3),
#             immune evasion (CD274/PD-L1, PDCD1LG2), CD15/FUT4, GATA3, ID2, LGALS1
#   NEGATIVE  B program is LOST (MS4A1, CD19, CD79A/B, CD22, FCER2), not T
#             (CD3D/E, TRAC, CD2), and CD45/PTPRC negative
#   MORPHOLOGY  Reed-Sternberg cells are giant and often multinucleated ->
#             large Area.um2 and NucArea
#
# Two tiers, because at ~500-700 counts/cell dropout makes any single-marker gate
# lose most true positives:
#   TIER A (strict seed): CD30+ AND >=1 corroborating positive AND B/T markers all
#          zero AND area in the upper quartile. High specificity, low sensitivity.
#   TIER B (score): per-sample z-scored composite (mean positive z - mean negative
#          z), top TOP_PCT within the sample, plus above-median area.
# Tier A is the preservation label if it is large enough; both are written out.
#
# Usage: sbatch --mem=200G scripts/sbatch/run_rscript.sh scripts/sbatch/label_hrs_candidates.R

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(data.table); setDTthreads(1)
})
OUT <- "output/04_merge"; FIG <- file.path(OUT, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
TOP_PCT <- 0.01

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
S_IB1<-"#2a78d6"; S_RB1<-"#eb6834"; CRIT<-"#d03b3b"
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
          legend.position="top", legend.justification="left",
          legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF);cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

CORE   <- c("TNFRSF8")
CORROB <- c("CCL17","IRF4","TRAF1","CD274","PDCD1LG2","FUT4","BIRC3","NFKB2","RELB","CD40",
            "TNFAIP3","GATA3","ID2","LGALS1")
NEG_B  <- c("MS4A1","CD19","CD79A","CD79B","CD22","FCER2")
NEG_T  <- c("CD3D","CD3E","TRAC","CD2")
NEG_L  <- c("PTPRC")
POS <- c(CORE, CORROB); NEG <- c(NEG_B, NEG_T, NEG_L)

obj <- readRDS(file.path(OUT, "merged_processed.Rds"))
md <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
cnt <- LayerData(obj[["RNA"]], layer = "counts")
genes <- intersect(c(POS, NEG), rownames(cnt))
cat("markers found:", length(genes), "of", length(c(POS,NEG)), "\n")
m <- as(cnt[genes, ], "dgCMatrix")     # ~40 x 1.8M, materialised off BPCells
rm(obj, cnt); invisible(gc())
stopifnot(identical(colnames(m), md$cell_id))
# meta.data carries 126 columns including its own `cell` column, which would
# collide with keep.rownames and break ggplot ("duplicate columns")
md <- md[, .(cell_id, merge_sample, Area.um2, NucArea, nCount_RNA, inSitu_clusters,
             lineage, strongB, pearson_clusters_low,
             pearson_clusters_high, x_slide_mm, y_slide_mm)]

# our own log-normalisation from counts (the merged object carries counts only)
tc <- pmax(md$nCount_RNA, 1)
ln <- t(log1p(t(m) / tc * median(tc)))
md[, samp := factor(merge_sample, levels = c("IB1","RB1"))]

# z-score WITHIN sample: IB1 and RB1 differ ~4x in detection, so a shared z would
# just re-encode sample identity
zmat <- matrix(0, nrow(ln), ncol(ln), dimnames = dimnames(ln))
for (s in levels(md$samp)) {
  i <- which(md$samp == s)
  sub <- as.matrix(ln[, i, drop = FALSE])
  mu <- rowMeans(sub); sdv <- apply(sub, 1, sd); sdv[sdv == 0] <- 1
  zmat[, i] <- (sub - mu) / sdv
}
pos_g <- intersect(POS, genes); neg_g <- intersect(NEG, genes)
md[, hrs_score := colMeans(zmat[pos_g, , drop=FALSE]) - colMeans(zmat[neg_g, , drop=FALSE])]

det <- function(g) as.vector(m[g, ] > 0)
md[, `:=`(cd30_pos = det("TNFRSF8"),
          n_corrob = colSums(m[intersect(CORROB, genes), , drop=FALSE] > 0),
          negB_zero = colSums(m[intersect(NEG_B, genes), , drop=FALSE] > 0) == 0,
          negT_zero = colSums(m[intersect(NEG_T, genes), , drop=FALSE] > 0) == 0)]
md[, area_q75 := Area.um2 >= quantile(Area.um2, .75), by = samp]
md[, area_med := Area.um2 >= median(Area.um2), by = samp]
md[, nuc_med  := NucArea >= median(NucArea), by = samp]

md[, hrs_tierA := cd30_pos & n_corrob >= 1 & negB_zero & negT_zero & area_q75 & nuc_med]
md[, score_cut := quantile(hrs_score, 1 - TOP_PCT), by = samp]
md[, hrs_tierB := hrs_score >= score_cut & area_med]
md[, hrs_label := ifelse(hrs_tierA, "tierA_strict",
                  ifelse(hrs_tierB, "tierB_score", "other"))]

cat("\n=== candidate counts (RB1 is the specificity control) ===\n")
print(md[, .(cells = .N,
             tierA = sum(hrs_tierA), tierA_pct = round(100*mean(hrs_tierA), 4),
             tierB = sum(hrs_tierB), tierB_pct = round(100*mean(hrs_tierB), 3),
             overlap = sum(hrs_tierA & hrs_tierB)), by = samp])
cat("\nIB1:RB1 tierA rate ratio =",
    round(md[samp=="IB1", mean(hrs_tierA)] / max(md[samp=="RB1", mean(hrs_tierA)], 1e-9), 2), "\n")

cat("\n=== marker detection %: tierA candidates vs rest, per sample ===\n")
dt <- rbindlist(lapply(levels(md$samp), function(s) {
  i <- md$samp == s
  rbindlist(lapply(genes, function(g) {
    v <- m[g, i] > 0; a <- md$hrs_tierA[i]
    data.table(sample=s, gene=g, set=ifelse(g %in% POS, "positive", "negative"),
               pct_in_tierA = round(100*mean(v[a]), 1),
               pct_in_rest  = round(100*mean(v[!a]), 1))
  }))
}))
dt[, ratio := round(pct_in_tierA / pmax(pct_in_rest, 0.01), 2)]
print(dt[sample == "IB1"][order(set, -ratio)], nrows = 50)
sc("hrs_marker_validation.csv", dt)

cat("\n=== morphology (IB1) ===\n")
print(md[samp == "IB1", .(n=.N, med_Area=round(median(Area.um2),1), med_NucArea=round(median(NucArea),1),
        med_nCount=median(nCount_RNA), pct_cluster_i=round(100*mean(inSitu_clusters=="i"),1)),
        by = .(hrs_label)])

cat("\n=== pre-integration cluster distribution of IB1 tierA ===\n")
for (g in c("pearson_clusters_low","pearson_clusters_high")) {
  x <- md[samp=="IB1" & hrs_tierA == TRUE, .N, by=c(g)][order(-N)]
  x[, pct := round(100*N/sum(N),1)]; cat("--", g, "\n"); print(head(x, 8))
}

keep <- c("cell_id","merge_sample","hrs_score","hrs_tierA","hrs_tierB","hrs_label",
          "cd30_pos","n_corrob","negB_zero","negT_zero","Area.um2","NucArea","nCount_RNA",
          "inSitu_clusters","lineage","strongB",
          "pearson_clusters_low","pearson_clusters_high","x_slide_mm","y_slide_mm")
sc("hrs_candidates.csv", md[hrs_label != "other", ..keep])
sc("hrs_scores_all_cells.csv", md[, .(cell_id, merge_sample, hrs_score, hrs_label)])

## ---- figures ----
p1 <- ggplot(md, aes(hrs_score, colour = samp)) +
  stat_density(geom = "line", linewidth = 0.7, position = "identity") +
  geom_vline(data = unique(md[, .(samp, score_cut)]), aes(xintercept = score_cut, colour = samp),
             linetype = "22", linewidth = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = c(IB1=S_IB1, RB1=S_RB1), name = NULL) +
  labs(title = "HRS composite score: IB1 vs RB1 (specificity control)",
       subtitle = paste0("Positive HRS program minus lost-B/T/CD45 program, z-scored within sample.\n",
                         "Dashed lines = each sample's top ", 100*TOP_PCT, "% cut (tier B)."),
       x = "HRS composite score", y = "density",
       caption = "Data: hrs_scores_all_cells.csv") + theme_viz()
sf("fig_04_hrs_score_distribution.png", p1, 8, 4.6)

dd <- dt[sample == "IB1"][order(set, -ratio)]
dd[, gene := factor(gene, levels = rev(gene))]
p2 <- ggplot(melt(dd[, .(gene, set, tierA = pct_in_tierA, rest = pct_in_rest)],
                  id.vars = c("gene","set"), variable.name = "grp", value.name = "pct"),
             aes(pct, gene, fill = grp)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(tierA = S_IB1, rest = MUTED), name = NULL,
                    labels = c("tier A candidates", "all other IB1 cells")) +
  labs(title = "Do the tier-A candidates actually look like HRS cells?",
       subtitle = "Marker detection rate in IB1. Positives should be far higher in candidates; negatives near zero by construction.",
       x = "% of cells detecting the gene", y = NULL,
       caption = "Data: hrs_marker_validation.csv") +
  theme_viz(xgrid = TRUE) + theme(panel.grid.major.y = element_blank(),
                                  strip.text.y = element_text(angle = 0))
sf("fig_04_hrs_marker_validation.png", p2, 9, 9)

ib <- md[samp == "IB1"]
p3 <- ggplot() +
  geom_point(data = ib[hrs_label == "other"][sample(.N, min(.N, 120000))],
             aes(x_slide_mm, y_slide_mm), colour = GRID, size = 0.05) +
  geom_point(data = ib[hrs_tierB == TRUE], aes(x_slide_mm, y_slide_mm),
             colour = "#9ec5f4", size = 0.35) +
  geom_point(data = ib[hrs_tierA == TRUE], aes(x_slide_mm, y_slide_mm),
             colour = CRIT, size = 0.9) +
  coord_equal() +
  labs(title = "Where the HRS candidates sit in IB1",
       subtitle = paste0("Red = tier A strict (n=", ib[, sum(hrs_tierA)], "); ",
                         "light blue = tier B score (n=", ib[, sum(hrs_tierB)], "); grey = other cells.\n",
                         "Clustered rather than scattered placement is what a real HRS population should look like."),
       x = "x (mm)", y = "y (mm)", caption = "Data: hrs_candidates.csv") +
  theme_viz() + theme(panel.grid.major.y = element_blank())
sf("fig_04_hrs_spatial_IB1.png", p3, 7.5, 8)
cat("\nOK\n")
