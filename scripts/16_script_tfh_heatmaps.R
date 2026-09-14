## ============================================================================
## 16_script_tfh_heatmaps.R -- marker heatmaps to pin down what "Tfh" is here.
##
## The Tfh label in this project is a marker-argmax call, not a cluster: cells first called
## `T` on CD3D/E/G/TRAC, then within those an argmax over {Tfh, Treg, CD8_T, CD4_conv}
## (scripts/16_script_celltypes.R). PDCD1 and CTLA4 are deliberately EXCLUDED from the
## definition so the PD-1 readout is not circular. These heatmaps are how a reader checks
## that call against the evidence instead of taking it on trust.
##
## TWO PANELS, both from the same pseudobulk:
##   A  columns = cell type (major) x sample (minor)   -- is the profile the same in both nodes?
##   B  columns = region (major) x cell type (minor)   -- is it the same in every region?
##
## For B the regional categories already encode the sample (IB1 FL_focus, RB1 FL_follicular,
## ...), because the two nodes do not share a region vocabulary -- RB1 has no Hodgkin hotspot
## and its light-chain readout is degenerate. So `sample` would be a redundant minor axis
## there and cell type is the informative one.
##
## MARKERS ARE GROUPED BY ROLE and the groups are drawn as separate row blocks:
##   defining     what the Tfh call was made ON -- high here is definitional, not evidence
##   held out     PDCD1, CTLA4, TOX -- the call never saw these, so THIS is the evidence
##   siblings     the other T subsets' defining markers -- Tfh should be low here
##   other lineage anchors -- Tfh should be low here too, or the call is picking up spillover
##
## SCALE: z-score of log1p(cp10k) per gene across the columns shown, because the markers span
## two orders of magnitude in absolute level and an unscaled heatmap would show only CD3E and
## MS4A1. The raw cp10k goes to CSV alongside. Columns with < MIN_CELLS cells are dropped, and
## every column carries its n in the axis label -- a 36-cell column and a 180,000-cell column
## should not look equally solid.
##
## Usage: Rscript scripts/16_script_tfh_heatmaps.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(ggplot2); library(patchwork); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER      <- if (length(args) >= 1) args[1] else "matched"
## Minimum cells for a column to be drawn. A pseudobulk of 36 cells and one of 180,000 look
## equally solid on a heatmap, which is the whole problem -- so this is a real parameter.
##
## 100 IS CHOSEN DELIBERATELY, not as a default. It is the largest cutoff that KEEPS Tfh in
## IB1's FL foci (125 cells), which is the column the FL question turns on, while still
## dropping the ones that are genuinely unusable: Tfh in IB1's cHL hotspot (36), Treg in IB1
## FL (48), CD8_T in IB1 FL (47), and cDC everywhere outside IB1 background (5-148). At 500
## the IB1 FL_focus and cHL_hotspot columns vanish from every block except B/Monocyte/fibro,
## which turns the figure into IB1-background-vs-RB1 and loses the point. The dropped columns
## and their n are printed either way, so the cost is never silent.
MIN_CELLS <- if (length(args) >= 2) as.integer(args[2]) else 100L

PB  <- file.path("output/16_ct_pseudobulk", TIER)
OUT <- file.path("output/16_tfh", TIER)
FIG <- file.path(OUT, "figures")
for (f in c("pb_regct_sum.Rds","regct_meta.csv.gz"))
  if (!file.exists(file.path(PB,f))) stop("missing ", file.path(PB,f), call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

CTS  <- c("Tfh","Treg","CD4_conv","CD8_T","T_unresolved","B","PC","Macrophage","cDC",
          "Monocyte","NK","endo","fibro")
REGS <- c("IB1|FL_focus","IB1|cHL_hotspot","IB1|background",
          "RB1|FL_follicular","RB1|interfollicular")
## The PD-1 axis gets its OWN row block rather than being split between "held out" and
## "ligands" as it was. It is the question of step 17, so the receptor and both ligands need
## to be readable together -- and all three are held out of every cell-type definition in this
## project, so high values here are a prediction the calls made, never a restatement.
MARKERS <- list(
  `PD-1 axis (held out)` = c("PDCD1","CD274","PDCD1LG2"),
  `Tfh: defining`   = c("CXCL13","CXCR5","ICOS","BTLA","CD40LG"),
  `other held out`  = c("CTLA4","TOX","TNFRSF9"),
  `Treg`            = c("FOXP3","IL2RA","IKZF2"),
  `CD4 / CD8`       = c("CD4","IL7R","CCR7","CD8A","GZMK","GZMB"),
  `T core`          = c("CD3D","CD3E","CD3G","TRAC"),
  `exhaustion`      = c("LAG3","TIGIT","HAVCR2"),
  `other lineages`  = c("MS4A1","CD79A","MZB1","CD68","CD163","BATF3","NKG7","PECAM1","COL1A1"))
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"
DIV_LO<-"#2a78d6"; DIV_MID<-"#f4f2ec"; DIV_HI<-"#eb6834"

cat("============================================================\n")
cat("step 16n -- Tfh marker heatmaps, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

S  <- readRDS(file.path(PB,"pb_regct_sum.Rds"))
GR <- fread(file.path(PB,"regct_meta.csv.gz")); setkey(GR, grp_reg); GR <- GR[colnames(S)]
stopifnot(identical(GR$grp_reg, colnames(S)))
GENES <- rownames(S)
MG <- rbindlist(lapply(names(MARKERS), function(k)
        data.table(block = k, gene = MARKERS[[k]])))[gene %in% GENES]
absent <- setdiff(unlist(MARKERS, use.names=FALSE), GENES)
if (length(absent)) message("not on the panel, dropped: ", paste(absent, collapse=", "))
MG[, block := factor(block, levels = names(MARKERS))]
MG[, gene  := factor(gene,  levels = MG$gene)]
cat(sprintf("markers shown: %d in %d blocks | pseudobulk columns: %d\n",
            nrow(MG), uniqueN(MG$block), ncol(S)))

## cp10k per column (share of that group's transcriptome), then z per gene across the
## columns actually drawn -- computed separately for A and B so each panel is self-scaled.
## as.matrix() on the OUTSIDE too: `dense %*% Diagonal()` returns a dgeMatrix, which
## as.data.table() cannot coerce, and the failure surfaces three steps later in the melt.
CP <- as.matrix(as.matrix(S) %*% Diagonal(x = 1e4/pmax(Matrix::colSums(S),1)))
dimnames(CP) <- dimnames(S)

zmat <- function(cp) {
  m <- log1p(cp)
  t(apply(m, 1, function(v) { s <- stats::sd(v); if (!is.finite(s) || s == 0) v*0 else (v-mean(v))/s }))
}

hm <- function(d, xlab, ttl, sub, ncol_major) {
  lim <- max(abs(d$z), na.rm = TRUE)
  ggplot(d, aes(minor, gene, fill = z)) +
    geom_tile(colour = SURF, linewidth = .35) +
    facet_grid(block ~ major, scales = "free", space = "free", switch = "y") +
    scale_fill_gradient2(low = DIV_LO, mid = DIV_MID, high = DIV_HI, midpoint = 0,
                         limits = c(-lim, lim), name = "z of log1p(cp10k)\nper gene") +
    scale_y_discrete(limits = rev) +
    labs(x = xlab, y = NULL, title = ttl, subtitle = sub,
         caption = sprintf("Data: output/16_tfh/%s/ | tier '%s' | columns with < %d cells dropped; n in each label",
                           TIER, TIER, MIN_CELLS)) +
    theme_minimal(base_size = 10) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid=element_blank(), panel.spacing=unit(2,"pt"),
          axis.text.x=element_text(colour=INK2, size=7, angle=90, hjust=1, vjust=.5),
          axis.text.y=element_text(colour=INK2, size=7.5),
          axis.title=element_text(colour=INK2, size=9),
          strip.text.x=element_text(colour=INK, face="bold", size=8),
          strip.text.y.left=element_text(colour=INK2, face="bold", size=7.5, angle=0, hjust=1),
          strip.placement="outside",
          plot.title=element_text(colour=INK, face="bold", size=13),
          plot.subtitle=element_text(colour=INK2, size=8.5),
          plot.caption=element_text(colour=MUTED, size=7, hjust=0),
          legend.position="right", legend.key.width=unit(9,"pt"))
}

## ---- A: cell type x sample --------------------------------------------------
## Aggregate over regions: sum the COUNTS then renormalise, never average the cp10k -- a
## mean of shares weights a 36-cell group equally with a 180,000-cell one.
cat(sprintf("\nminimum cells per column: %d\n", MIN_CELLS))
A_keys <- CJ(samp = c("IB1","RB1"), celltype = CTS, unique = TRUE)
A_n <- rbindlist(lapply(seq_len(nrow(A_keys)), function(i) {
  jj <- GR[samp == A_keys$samp[i] & celltype == A_keys$celltype[i], which = TRUE]
  data.table(samp = A_keys$samp[i], celltype = A_keys$celltype[i],
             cells = if (length(jj)) sum(GR$cells[jj]) else 0L)
}))
cat("\n=== panel A columns DROPPED at this cutoff\n")
if (A_n[cells < MIN_CELLS, .N]) print(A_n[cells < MIN_CELLS][order(-cells)]) else cat("  none\n")
A <- rbindlist(lapply(seq_len(nrow(A_keys)), function(i) {
  jj <- GR[samp == A_keys$samp[i] & celltype == A_keys$celltype[i], which = TRUE]
  if (!length(jj)) return(NULL)
  n <- sum(GR$cells[jj]); if (n < MIN_CELLS) return(NULL)
  tot <- sum(GR$nCount[jj])
  data.table(samp = A_keys$samp[i], celltype = A_keys$celltype[i], cells = n,
             gene = GENES, cp10k = 1e4 * as.numeric(Matrix::rowSums(S[, jj, drop=FALSE])) / tot)
}))
sc("marker_cp10k_celltype_by_sample.csv", A[gene %in% MG$gene][order(celltype, samp, gene)])
AW <- dcast(A, gene ~ celltype + samp, value.var = "cp10k", sep = "|")
AM <- as.matrix(AW[, -1]); rownames(AM) <- AW$gene
AM <- AM[as.character(MG$gene), , drop = FALSE]
AZ <- zmat(AM)
nlab <- unique(A[, .(key = paste(celltype, samp, sep="|"), cells)])
dA <- melt(as.data.table(AZ, keep.rownames="gene"), id.vars="gene",
           variable.name="key", value.name="z")
dA[, `:=`(major = sub("\\|.*$","",key), minor = sub("^.*\\|","",key))]
dA <- merge(dA, nlab, by="key")
dA[, minor := sprintf("%s (%s)", minor, formatC(cells, big.mark=",", format="d"))]
dA[, gene := factor(gene, levels = levels(MG$gene))]
dA <- merge(dA, MG, by="gene")
dA[, major := factor(major, levels = intersect(CTS, unique(major)))]
pA <- hm(dA, "sample (cells)", "What is a \"Tfh\" cell here? Marker profile by cell type and node",
  paste0("Cell type is the marker-argmax call from step 16b. The HELD OUT block is the evidence: ",
         "PDCD1 and CTLA4 were excluded from the\ndefinition, so Tfh being high there is a ",
         "prediction the call made and passed, not a restatement of it."), length(CTS))
ggsave(file.path(FIG,"fig_16_tfh_markers_by_celltype.png"), pA, width=15, height=8.5, dpi=170, bg=SURF)
cat("wrote ", file.path(FIG,"fig_16_tfh_markers_by_celltype.png"), "\n", sep="")

## ---- B: region x cell type --------------------------------------------------
GR[, reg_key := paste(samp, bin_region, sep="|")]
Ball <- GR[reg_key %in% REGS & celltype %in% CTS]
cat("\n=== panel B columns DROPPED at this cutoff\n")
if (Ball[cells < MIN_CELLS, .N])
  print(dcast(Ball[cells < MIN_CELLS], celltype ~ reg_key, value.var = "cells", fill = NA)) else
  cat("  none\n")
## Combinations that do not appear at all -- distinct from ones dropped for small n.
miss <- CJ(celltype = CTS, reg_key = REGS)[!Ball, on = c("celltype","reg_key")]
if (nrow(miss)) { cat("\n  combinations with NO cells at all:\n"); print(miss) }
B <- Ball[cells >= MIN_CELLS]
cat(sprintf("\npanel B columns kept: %d of %d region x celltype groups\n",
            nrow(B), nrow(Ball)))
print(dcast(B[, .N, by=.(reg_key, celltype)], celltype ~ reg_key, value.var="N", fill=0L))
BM <- CP[as.character(MG$gene), B$grp_reg, drop=FALSE]
colnames(BM) <- paste(B$reg_key, B$celltype, sep="@")
BZ <- zmat(BM)
dB <- melt(as.data.table(BZ, keep.rownames="gene"), id.vars="gene",
           variable.name="key", value.name="z")
## Major = CELL TYPE, minor = node+region. So each block holds one cell type and you read
## its profile across the five regional categories side by side -- which is the arrangement
## that answers "is Tfh the same thing everywhere", rather than "what is in each region".
dB[, `:=`(major = sub("^.*@","",key), minor = sub("@.*$","",key))]
dB <- merge(dB, B[, .(key = paste(reg_key, celltype, sep="@"), cells)], by="key")
## Short region labels: 13 cell types x 5 regions is 65 columns, and the full
## "IB1|cHL_hotspot" strings do not fit on an axis at that density.
REG_SHORT <- c(`IB1|FL_focus`="IB1 FL", `IB1|cHL_hotspot`="IB1 hot",
               `IB1|background`="IB1 bg", `RB1|FL_follicular`="RB1 foll",
               `RB1|interfollicular`="RB1 inter")
dB[, minor_lab := sprintf("%s (%s)", REG_SHORT[minor], formatC(cells, big.mark=",", format="d"))]
dB[, minor_lab := factor(minor_lab,
     levels = unique(dB[order(match(minor, REGS)), minor_lab]))]
dB[, minor := minor_lab][, minor_lab := NULL]
dB[, gene := factor(gene, levels = levels(MG$gene))]
dB <- merge(dB, MG, by="gene")
dB[, major := factor(major, levels = intersect(CTS, unique(major)))]
sc("marker_cp10k_region_by_celltype.csv",
   melt(as.data.table(BM, keep.rownames="gene"), id.vars="gene",
        variable.name="region_celltype", value.name="cp10k"))
pB <- hm(dB, "node + region (cells)", "Is the Tfh profile the same in every region?",
  paste0("Major columns are cell type, minor are node+region -- so each block is one cell ",
         "type read across the five regional\ncategories. The regions encode the sample ",
         "because the two nodes share no region vocabulary (RB1 has no Hodgkin\nhotspot, and ",
         "its light-chain readout is degenerate). IB1 hot has only 36 Tfh cells: read that ",
         "column as noise."),
  length(CTS))
ggsave(file.path(FIG,"fig_16_tfh_markers_by_region.png"), pB, width=19.5, height=8.5, dpi=170, bg=SURF)
cat("wrote ", file.path(FIG,"fig_16_tfh_markers_by_region.png"), "\n", sep="")

## ---- the numbers behind the headline claim ---------------------------------
cat("\n=== Tfh vs its siblings on the HELD-OUT markers, cp10k (panel A columns)\n")
sib <- c("Tfh","Treg","CD4_conv","CD8_T","T_unresolved")
print(dcast(A[celltype %in% sib & gene %in% c(MARKERS$`PD-1 axis (held out)`,
                                              MARKERS$`other held out`)],
            gene ~ celltype + samp, value.var="cp10k")[, lapply(.SD, function(v)
              if (is.numeric(v)) round(v,3) else v)])
cat("\n=== and on the DEFINING markers (high here is definitional)\n")
print(dcast(A[celltype %in% sib & gene %in% MARKERS$`Tfh: defining`],
            gene ~ celltype + samp, value.var="cp10k")[, lapply(.SD, function(v)
              if (is.numeric(v)) round(v,3) else v)])
writeLines(format(Sys.time()), file.path(OUT,"tfh_heatmaps.done"))
cat("============================================================\nOK\n")
