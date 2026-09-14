## ============================================================================
## 19_script_b_substructure.R -- what B substructure does the 6K panel actually support?
##
## WHY. The gate passed, but 17 of 27 high-resolution clusters have modal argmax call `B`.
## That is not an error: B is 50.5% of all argmax-resolved cells in the high-confidence
## object, because RB1 is follicular lymphoma and supplies 62.6% of it. But labelling by
## modal coarse call would name most of the object "B cell" and stop, which is useless.
## So: is there panel evidence to split them, and does it vary between the B clusters?
##
## PRESENCE IS NOT SUPPORT. A module whose genes are on the panel but flat across every B
## cluster cannot subtype anything. This reports the RANGE across clusters, not just the
## level, and gates on cp10k AND detection rate moving together -- cp10k is compositional
## and fakes depression in whatever else is abundant, which in a plasma-cell-rich or
## light-chain-saturated cluster is a real risk (CLAUDE.md fact 8).
##
## cp10k HERE IS THE POOLED SHARE of the cluster's transcriptome (the step-11 / bin-table
## convention), NOT the mean of per-cell shares (the step-10 per-gene convention). The two
## are not comparable and mixing them is a documented trap.
##
## SUBTYPING WITHIN B IS A LEGITIMATE RELATIVE COMPARISON, unlike the step-09b across-cluster
## z that this project had to abandon. 09b failed because it used relative scoring to decide
## IDENTITY when one lineage dominated the object -- B markers stopped being distinctive
## between clusters precisely because everything was B. Here identity is already settled and
## the question is which B cluster is most GC-like, so comparing across B clusters is the
## right frame. Absolute cp10k is reported alongside anyway, so a module that is "highest"
## at a level indistinguishable from zero is visible as such.
##
## COUNTS COME FROM THE OBJECT, NOT FROM panel_counts_all.Rds. That file holds only 103
## genes -- the targeted subset step 16 built, which is also why panel_detection_preqc.csv
## has 103 rows and not 6,175. Using it here silently dropped 29 of the candidate markers,
## including the whole isotype axis (IGHM, IGHG1/2, IGHA1) and most of naive/mantle (IGHD,
## TCL1A, FCER2), leaving those modules to be scored on whatever fragment survived. The
## full 6,175-gene counts exist only in the BPCells matrices, so this reads them from the
## step-19a object.
##
## TWO PANEL CAVEATS, both from CLAUDE.md:
##   * CR2 and FCER2 are B markers, NOT FDC-exclusive -- that is why FDC was dropped from
##     LINEAGE_SETS_COARSE. A naive/mantle score built on them will also light up on FDC.
##   * RB1's light chain is degenerate: 178/178 neighbourhoods kappa-restricted. The
##     kappa/lambda axis cannot separate FL clone from normal B in RB1 -- the whole node is
##     clone. It is reported for IB1's sake and as a negative control for RB1.
##
## Usage: sbatch --mem=150G --cpus-per-task=2 --time=2:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_b_substructure.R [low|high]
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); setDTthreads(2)
})
args <- commandArgs(trailingOnly = TRUE)
RES  <- if (length(args) >= 1) args[1] else "high"
CL   <- paste0("hc_clusters_", RES)
OUT  <- "output/19_celltypes"; FIG <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
HC_F <- file.path(OUT, "hc_cells.csv.gz")
OBJ_F<- file.path(OUT, "hc_02_clustered.Rds")
CT_F <- "output/16_celltypes/cell_celltypes.csv.gz"
for (f in c(HC_F, OBJ_F, CT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}

MODULES <- list(
  B_core       = c("MS4A1","CD79A","CD79B","BANK1"),
  germinal_ctr = c("BCL6","AICDA","RGS13","LMO2","MME","EZH2","BACH2"),
  dark_zone    = c("MKI67","TOP2A","PCNA","CCNB1","TYMS","RRM2","CDK1","UBE2C","BIRC5","CXCR4"),
  light_zone   = c("CD83","CD40","CXCR5","ICAM1","SLAMF1","PSMB10","NFKB1"),
  ## naive/mantle holds ONLY B-restricted genes. CR2 and FCER2 are split out below: this
  ## project already records that "CR2/FCER2 are B markers; nothing FDC-exclusive survives
  ## on this panel", and a naive/mantle score containing them labelled an FDC cluster as
  ## naive B -- cluster 17, CR2 22.9 cp10k against 1.9 in the clonal core, while its own
  ## MS4A1 and CD79A sat BELOW core. CCR7/IL4R/YBX3 are dropped as not B-restricted either.
  naive_mantle = c("IGHD","TCL1A","SELL"),
  fdc_shared   = c("CR2","FCER2"),
  memory       = c("CD27","TNFRSF13B"),
  plasma       = c("PRDM1","XBP1","IRF4","MZB1","DERL3","SDC1","TNFRSF17","SEC11C","JCHAIN","POU2AF1"),
  FL_BCL2      = c("BCL2"),
  isotype_M    = c("IGHM"),
  isotype_G    = c("IGHG1/2"),
  isotype_A    = c("IGHA1"),
  kappa        = c("IGKC"),
  lambda       = c("IGLC1/2"),
  ## ---- supplementary modules, scored for EVERY cluster, not just the B ones ----
  ## Two of the largest unnamed clusters turned out to be real populations that
  ## LINEAGE_SETS_COARSE simply has no name for, so argmax could never call them and they
  ## sat in `mixed`. These modules let 19d name them at the cluster level.
  ##
  ## pDC: argmax splits these cells evenly between the PC and B sets (purity 0.227 vs 0.225,
  ## a margin of 0.002) because plasmacytoid DCs express MZB1 and IRF4 -- the same confusion
  ## CLAUDE.md records for HierAnnot, "pDC called plasma cell". CLEC4C (CD303), LILRA4 (ILT7)
  ## and IL3RA (CD123) are pDC-restricted and settle it.
  pDC          = c("CLEC4C","LILRA4","IL3RA","IRF7","PLD4","GZMB","IRF8","SERPINF1",
                   "PTPRS","ITM2C"),
  ## Hodgkin niche: NOT a cell type. CCL17/CCL22 are the project's established cHL readout
  ## (tissue-wide field plus 36 focal hotspots in IB1, zero RB1 equivalents). Per-cell HRS
  ## calling is a CLOSED avenue here -- TNFRSF8 maxes at 6 counts with 78% singletons -- so
  ## a cluster scoring high on this is a REGION, and 19d labels it as one.
  hodgkin_niche= c("CCL17","CCL22","BATF3","IL22RA2"))

cat("============================================================\n")
cat("step 19 -- B substructure on the 6K panel (", CL, ")\n", sep="")
cat("============================================================\n")

hc <- fread(HC_F)
ct <- fread(CT_F, select = c("cell_id","lineage16"))
hc <- ct[hc, on = .(cell_id)]
setnames(hc, CL, "cl")
## split() names its output with characters, so per-gene rows key on a character cluster id.
## Keep every cluster id character from here on, or the joins fail on type.
hc[, cl := as.character(cl)]

cat("reading ", OBJ_F, "\n", sep = "")
seu <- readRDS(OBJ_F)
ca  <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
cat(sprintf("object counts: %d genes x %d cells (%s)\n", nrow(ca), ncol(ca), class(ca)[1]))

miss <- setdiff(unlist(MODULES, use.names = FALSE), rownames(ca))
if (length(miss)) { cat("NOT on the panel, dropped: ", paste(miss, collapse=", "), "\n", sep="")
  MODULES <- lapply(MODULES, intersect, rownames(ca))
  MODULES <- MODULES[lengths(MODULES) > 0] }
GEN <- unique(unlist(MODULES, use.names = FALSE))
cat(sprintf("module genes found: %d\n", length(GEN)))

## Order hc to the object's cells, then materialise only the module genes.
setkey(hc, cell_id); hc <- hc[J(colnames(ca))]
if (anyNA(hc$cl)) stop(sum(is.na(hc$cl)), " object cells missing from hc_cells.csv.gz", call.=FALSE)
tot <- seu$nCount_RNA
sub <- as(as.matrix(ca[GEN, , drop = FALSE]), "dgCMatrix")
rownames(sub) <- GEN
rm(seu, ca); invisible(gc())
cat(sprintf("materialised %d x %d\n", nrow(sub), ncol(sub)))

## Which clusters are B? modal argmax call among resolved cells.
FAMS <- c("B","T","NK","myeloid","endo","fibro","PC")
mode_call <- hc[lineage16 %in% FAMS, .N, by = .(cl, lineage16)][
  order(-N), .(top = lineage16[1], purity = round(N[1]/sum(N),3)), by = cl]
cat("\n=== clusters by modal argmax call ===\n"); print(mode_call[order(top, -purity)])
BCL <- mode_call[top == "B", cl]
cat(sprintf("\nB-dominant clusters: %d of %d\n", length(BCL), nrow(mode_call)))

## Pooled cp10k per cluster per gene, and detection rate. Pooled: sum of the gene's counts
## over the cluster divided by the cluster's total transcripts.
idx <- split(seq_len(ncol(sub)), hc$cl)
per_gene <- rbindlist(lapply(names(idx), function(k) {
  j <- idx[[k]]; s <- sub[, j, drop = FALSE]
  data.table(cl = k, gene = rownames(sub),
             cp10k = as.numeric(Matrix::rowSums(s)) / sum(tot[j]) * 1e4,
             det   = as.numeric(Matrix::rowSums(s > 0)) / length(j))
}))
per_gene[, module := {
  m <- rep(NA_character_, .N)
  for (nm in names(MODULES)) m[gene %in% MODULES[[nm]]] <- nm
  m }]
sc(sprintf("b_substructure_genes_%s.csv", RES), per_gene)

mod <- per_gene[!is.na(module), .(cp10k = mean(cp10k), det = mean(det)), by = .(cl, module)]
mod <- mode_call[mod, on = .(cl)]
## supplementary modules are reported for every cluster, keyed off the object-wide median
sup <- mod[module %in% c("pDC","hodgkin_niche")]
sup[, ratio := cp10k / median(cp10k), by = module]
cat("\n=== supplementary modules: clusters scoring >3x the object median ===\n")
print(sup[ratio > 3][order(module, -ratio),
          .(module, cl, top, cp10k = round(cp10k,2), ratio = round(ratio,2))])
mod[, cl := factor(cl, levels = c(mode_call[order(-purity)][top == "B", cl],
                                  setdiff(mode_call$cl, BCL)))]

## Does the module actually VARY across the B clusters? A module flat across them cannot
## subtype anything, however respectable its gene list.
cat("\n=== module range ACROSS B-dominant clusters (pooled cp10k) ===\n")
rng <- mod[cl %in% BCL, .(lo = min(cp10k), hi = max(cp10k),
                          fold = round(max(cp10k)/pmax(min(cp10k),1e-6), 2),
                          det_lo = round(min(det),3), det_hi = round(max(det),3),
                          det_fold = round(max(det)/pmax(min(det),1e-6), 2)), by = module]
rng[, `:=`(lo = round(lo,3), hi = round(hi,3))]
rng[, usable := fold >= 2 & det_fold >= 1.5 & hi >= 0.5]
print(rng[order(-fold)])
cat("\n`usable` = varies at least 2-fold in cp10k AND 1.5-fold in detection AND reaches\n")
cat("cp10k 0.5 somewhere. cp10k and detection must move together: cp10k alone is\n")
cat("compositional and fakes depression wherever something else is abundant.\n")
sc(sprintf("b_substructure_modules_%s.csv", RES), mod)
sc(sprintf("b_substructure_range_%s.csv", RES), rng)

cat("\n=== module cp10k per B-dominant cluster ===\n")
w <- dcast(mod[cl %in% BCL], cl ~ module, value.var = "cp10k")
num <- setdiff(names(w), "cl"); w[, (num) := lapply(.SD, round, 2), .SDcols = num]
print(w)

## kappa/lambda, reported per cluster. Degenerate in RB1 by construction; the point is IB1.
lc <- per_gene[gene %in% c("IGKC","IGLC1/2"), .(cl, gene, cp10k)]
lc <- dcast(lc, cl ~ gene, value.var = "cp10k")
if (all(c("IGKC","IGLC1/2") %in% names(lc))) {
  lc[, kappa_frac := IGKC / pmax(IGKC + `IGLC1/2`, 1e-9)]
  cat("\n=== light chain per cluster (RB1 is 178/178 kappa-restricted; IB1 is the point) ===\n")
  print(lc[cl %in% BCL][order(-kappa_frac)])
  sc(sprintf("b_substructure_lightchain_%s.csv", RES), lc)
}

## ---- heatmap ---------------------------------------------------------------
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"
hm <- mod[cl %in% BCL]
hm[, z := (cp10k - mean(cp10k)) / pmax(sd(cp10k), 1e-9), by = module]
p <- ggplot(hm, aes(factor(cl), module, fill = z)) + geom_tile(colour = SURF, linewidth = .4) +
  scale_fill_gradient2(low = "#184f95", mid = "#f0efec", high = "#d03b3b", midpoint = 0,
                       limits = c(-2,2), oob = scales::squish, name = "z across B clusters") +
  labs(title = sprintf("B substructure across the %d B-dominant clusters (%s)", length(BCL), CL),
       subtitle = "Pooled cp10k per module, z-scored across B clusters only. Identity is already settled here,\nso a relative comparison is legitimate -- see the range table for whether each module varies at all.",
       x = "cluster", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.background=element_rect(fill=SURF,colour=NA), panel.grid=element_blank(),
        axis.text=element_text(colour=MUTED), axis.title=element_text(colour=INK2),
        plot.title=element_text(colour=INK,face="bold"),
        plot.subtitle=element_text(colour=INK2,size=9), legend.position="top")
sf(sprintf("b_substructure_%s.png", RES), p, 12, 5.5)
cat("\ndone\n")
