## ============================================================================
## 19_script_halo_qc.R -- is any step-19 label a neighbourhood artefact?
##
## Step 18 established that cluster assignment is contaminated by what a cell sits NEXT TO,
## and that the clustering is the MORE contaminated of the two callers. Holding a cell's own
## PC-marker count at zero, the cluster called it PC 30.1% of the time within 10 um of a
## plasma cell and 3.6% beyond 80 um -- an 8.4x gradient on proximity alone. Across that
## range the cluster's PC rate moved 6.8x against argmax's 2.5x. The mechanism is specific:
## the leaked mass is light chains, which whole-transcriptome clustering sees and argmax's
## narrow PC marker set excludes.
##
## Every label written by 19d therefore has to answer the same question: does membership
## depend on the NEIGHBOUR's marker load once the cell's OWN marker count is held fixed? A
## real cell type does not. A halo does.
##
## THIS IS THE ONLY LEVEL AT WHICH LEAKAGE CAN BE CORRECTED. Per-cell leakage scores are all
## at chance -- positional displacement AUC 0.446-0.526, donor load 0.354-0.605, compartment
## discrepancy 0.383-0.554. The myeloid effect is +0.127 log2 across 5,913 cells and is
## invisible in any single one of them. So there is no per-cell filter to apply and none is
## attempted here; the diagnostic is per label, where thousands of cells average the noise
## away and the 8.4x signal is unmissable.
##
## DONORS ARE DEFINED PHYSICALLY, as in 18_script_select_cells.R: the nearest cell carrying
## >= MIN_DONOR transcripts of family F, searched over ALL non-defect cells rather than the
## labelled subset. Only 2.05M of 4.4M cells are in the high-confidence object, and
## restricting the donor pool to those would put "nearest" three cell diameters away.
## Geometry is canvas px from the *_centroids.csv.gz tables, so the fitted mm->px affine
## never enters.
##
## Usage: sbatch --mem=150G --cpus-per-task=2 --time=3:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_halo_qc.R [high|low]
## ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(ggplot2); library(RANN); setDTthreads(2)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/lineage.R")

args <- commandArgs(trailingOnly = TRUE)
RES  <- if (length(args) >= 1) args[1] else "high"
OUT  <- "output/19_celltypes"; FIG <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
LAB_F <- file.path(OUT, "hc_cells_final.csv.gz")
if (!file.exists(LAB_F)) LAB_F <- file.path(OUT, sprintf("hc_cells_labelled_%s.csv.gz", RES))
MD_F  <- "output/16_regions/cells_meta_all.csv.gz"
CNT_F <- "output/16_regions/panel_counts_all.Rds"
for (f in c(LAB_F, MD_F, CNT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}

UMPX      <- 120.280945e-6 * 1000
MIN_DONOR <- 2L
D_MAX_UM  <- 20
FAMS      <- names(LINEAGE_SETS_COARSE)
LAB2FAM   <- c("B cell" = "B", "Plasmablast" = "PC",
               "B cell, naive/mantle" = "B", "T cell" = "T", "NK cell" = "NK",
               "Myeloid" = "myeloid", "Endothelial" = "endo",
               "Stromal/FDC" = "fibro", "Plasma cell" = "PC")

cat("============================================================\n")
cat("step 19e -- halo QC: is any label driven by the neighbour? (", RES, ")\n", sep="")
cat("============================================================\n")

lab <- fread(LAB_F, na.strings = c("NA",""))
md  <- fread(MD_F, select = c("cell_id","samp","is_defect"))
cells <- md[is_defect == FALSE]
cat(sprintf("labelled HC cells %d | non-defect donor universe %d\n", nrow(lab), nrow(cells)))

## centroids, canvas px -> physical um (y negated: canvas px runs opposite to slide mm)
ce <- rbindlist(lapply(c("IB1","RB1"), function(S) {
  f <- sprintf("output/18_leakage/%s_centroids.csv.gz", S)
  if (!file.exists(f)) stop("missing ", f, " -- run scripts/18_script_centroids.sh ", S, call.=FALSE)
  fread(f, select = c("cell_id","cx_px","cy_px"))[, samp := S][] }))
cells <- ce[cells, on = .(cell_id), nomatch = NULL]
cells[, `:=`(X_um = cx_px * UMPX, Y_um = -cy_px * UMPX)]
cat(sprintf("with centroids: %d\n", nrow(cells)))

## per-cell marker load per family
M <- readRDS(CNT_F)
gi <- match(cells$cell_id, colnames(M))
if (anyNA(gi)) stop(sum(is.na(gi)), " cells absent from the count matrix", call. = FALSE)
M <- M[, gi, drop = FALSE]
for (f in FAMS) {
  g <- intersect(LINEAGE_SETS_COARSE[[f]], rownames(M))
  cells[[paste0("k_", f)]] <- as.integer(Matrix::colSums(M[g, , drop = FALSE]))
}
rm(M); invisible(gc())

## nearest F-donor, per sample (a cell can only receive from a cell on its own slide)
cat("\n=== nearest-donor search ===\n")
for (f in FAMS) {
  kc <- paste0("k_", f)
  cells[[paste0("dd_", f)]] <- NA_real_; cells[[paste0("dk_", f)]] <- NA_integer_
  for (S in c("IB1","RB1")) {
    qi <- which(cells$samp == S)
    di <- which(cells$samp == S & cells[[kc]] >= MIN_DONOR)
    if (length(di) < 2L) next
    r <- RANN::nn2(cells[di, .(X_um, Y_um)], cells[qi, .(X_um, Y_um)], k = 2L)
    self <- cells$cell_id[di][r$nn.idx[,1]] == cells$cell_id[qi]
    idx  <- ifelse(self, r$nn.idx[,2], r$nn.idx[,1])
    dst  <- ifelse(self, r$nn.dists[,2], r$nn.dists[,1])
    set(cells, qi, paste0("dd_", f), dst)
    set(cells, qi, paste0("dk_", f), cells[[kc]][di][idx])
  }
  cat(sprintf("  %-8s donors %8d  median nn %.1f um\n", f,
              cells[get(kc) >= MIN_DONOR, .N], median(cells[[paste0("dd_", f)]], na.rm=TRUE)))
}

d <- cells[lab[, .(cell_id, celltype_cluster, cluster_status)], on = .(cell_id), nomatch = NULL]
cat(sprintf("\nHC cells with geometry and counts: %d\n", nrow(d)))

## ---- the diagnostic --------------------------------------------------------
## For each label L of family F: P(cell carries label L) as the nearest F-donor's load
## rises, with the cell's OWN F-marker count held fixed. Flat = a real cell type.
res <- rbindlist(lapply(names(LAB2FAM), function(L) {
  f <- LAB2FAM[[L]]; if (!L %in% d$celltype_cluster) return(NULL)
  x <- d[!is.na(get(paste0("dd_", f))) & get(paste0("dd_", f)) <= D_MAX_UM]
  if (!nrow(x)) return(NULL)
  x[, own := pmin(get(paste0("k_", f)), 3L)]
  x[, dq  := cut(get(paste0("dk_", f)), c(-Inf,2,3,5,9,Inf),
                 labels = c("2","3","4-5","6-9","10+"))]
  x[!is.na(dq), .(label = L, family = f, n = .N, p = 100*mean(celltype_cluster == L, na.rm=TRUE)),
    by = .(own, dq)]
}))
res <- res[n >= 200]
sc(sprintf("halo_qc_%s.csv", RES), res)

cat("\n=== P(cell carries the label) %, by the DONOR's marker load, at fixed own load ===\n")
for (L in unique(res$label)) {
  cat("\n--", L, "\n")
  print(dcast(res[label == L], own ~ dq, value.var = "p"))
}

## halo score: fold change from the lowest to the highest donor-load bin, at own = 0
h <- res[own == 0, .(lo = p[dq == "2"], hi = p[dq == "10+"],
                     n_lo = n[dq == "2"], n_hi = n[dq == "10+"]), by = .(label, family)]
h[, fold := round(hi / pmax(lo, 1e-6), 2)]
## `needs_conditioning`, NOT `halo_suspect`. A high fold here identifies a label to
## INVESTIGATE, not one to distrust: this test cannot separate leakage from genuine spatial
## aggregation, and on the pooled Plasma cell label the 12.69x flag turned out to be
## aggregation. Emitting a boolean called halo_suspect straight off the fold is what let that
## false positive propagate into the per-cell tables in the first place.
h[, needs_conditioning := fold >= 3 & n_lo >= 500 & n_hi >= 200]
## Verdicts recorded once a label has actually been conditioned on independent evidence.
COND <- c("Plasma cell" = "conditioned: spatial aggregation, not leakage (pc_halo_test.csv)",
          "Plasmablast" = "conditioned: spatial aggregation, not leakage (pc_halo_test.csv)")
h[, conditioned_verdict := COND[label]]
h[, halo_suspect := needs_conditioning & is.na(conditioned_verdict)]
cat("\n=== halo score: cells with NO marker of their own, lowest vs highest donor load ===\n")
print(h[order(-fold), .(label, family, fold, n_lo, n_hi, needs_conditioning,
                        halo_suspect, conditioned_verdict)])
sc(sprintf("halo_score_%s.csv", RES), h)
cat("\nFor reference, step 18 measured 8.4x for the InSituType PC call. A label with a LOW fold\n")
cat("here is not being assembled out of its neighbours.\n")
cat("\n*** A HIGH FOLD IS NOT A VERDICT. *** This test cannot separate leakage from genuine\n")
cat("spatial aggregation, and for plasma cells the difference was decisive: the 12.69x flag\n")
cat("on the pooled Plasma cell label was shown to be aggregation, not leakage --\n")
cat("PC-labelled cells with no plasma marker of their own carry HALF the B markers and 5-10x\n")
cat("the light chain of non-PC-labelled cells in the SAME donor-load bins, and conditioning\n")
cat("on their own light chain collapsed the gradient from 12.62x to 1.62x. Any label flagged\n")
cat("here needs that conditioning before the flag is acted on. See pc_halo_test.csv and\n")
cat("19_script_pc_refine.R.\n")

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
## POINT SIZE CARRIES n. The fold statistic is hi/lo and is therefore dominated by the
## rightmost bin, whose size differs by more than an order of magnitude between labels --
## 1,569 cells for plasma against 29,470 for myeloid. Without n on the figure, a rare-and-
## extreme tail and a common gradient look like the same result.
p <- ggplot(res[own <= 2], aes(dq, p, colour = factor(own), group = own)) +
  geom_line(linewidth=.6) + geom_point(aes(size=n)) +
  facet_wrap(~ label, scales="free_y") +
  scale_colour_manual(values=c("0"="#184f95","1"="#7b9fc9","2"="#d03b3b"),
                      name="the cell's OWN markers") +
  scale_size_continuous(trans="log10", range=c(.6,3.4), name="cells in bin",
                        labels=scales::label_number(scale_cut=scales::cut_short_scale())) +
  labs(title = "How spatially aggregated is each label? (this is NOT a leakage detector)",
       subtitle = paste0(
  "RISING LINES ARE EXPECTED, not an alarm. Cell types cluster in tissue, so a cell with no markers of its own sitting beside a\n",
  "marker-rich neighbour genuinely IS more likely to be that type with dropout. Proven here: the plasma-cell rise looked like a\n",
  "12.6x halo until conditioned on the cell's own light chain, which collapsed it to 1.6x -- aggregation, not leakage.\n",
  "READ INSTEAD: (1) the VERTICAL GAP between the lines -- own evidence beats neighbour evidence in every panel;\n",
  "(2) the SHAPE -- a steady rise (myeloid, stromal, endothelial) differs from a flat line that hooks only in the last bin (plasma);\n",
  "(3) POINT SIZE -- the fold is hi/lo, so it is set by the rightmost bin, which holds 1,569 cells for plasma and 29,470 for myeloid."),
       x = "marker load of the nearest same-family cell within 20um", y = "% carrying the label",
       caption = paste0(
  "A high fold marks a label to INVESTIGATE (needs_conditioning), never one to distrust. Converting a rise into a leakage claim ",
  "requires conditioning on\nindependent evidence of the cell's own identity -- see 19_script_pc_refine.R and pc_halo_test.csv for the worked example.")) +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.minor=element_blank(), panel.grid.major=element_line(colour=GRID,linewidth=.25),
        axis.text=element_text(colour=MUTED), axis.title=element_text(colour=INK2),
        strip.text=element_text(colour=INK2,face="bold"), plot.title=element_text(colour=INK,face="bold"),
        plot.subtitle=element_text(colour=INK2,size=8),
        plot.caption=element_text(colour=MUTED,size=7,hjust=0), legend.position="top")
sf(sprintf("halo_qc_%s.png", RES), p, 14, 8.5)
cat("\ndone\n")
