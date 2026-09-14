## ============================================================================
## 07_script_lightchain_neighborhoods.R -- light-chain neighbourhoods from two
## genes only: IGKC and IGLC1/2.
##
##   1. find cells with HIGH light-chain expression
##   2. build a spatial graph connecting nearby light-chain cells
##   3. find neighbourhoods (graph communities)
##   4. score each neighbourhood's kappa:lambda ratio
##
## Context: IB1 is an FL composite with HRS; RB1 is a recurrence of the same FL and
## is almost entirely IGKC. So kappa-restricted neighbourhoods in IB1 = FL;
## polytypic neighbourhoods = reactive germinal centres. Graph communities are used
## instead of square bins because follicles are discrete structures the graph can
## find itself.
##
## Run this BEFORE scripts/07_script_lightchain_views.R, which renders image views
## of what this produces. Submit with scripts/07_submit_lightchain_neighborhoods.sh.
##
## ---------------------------------------------------------------------------
## MOVED from scripts/sbatch/lightchain_neighborhoods.R, which produced the tables
## still sitting in output/04_merge/. Differences, all deliberate:
##
##   * INPUT is now step 1's PRE-QC-FILTER object,
##       output/01_import_v1/seuratObject_<sample>_BPCells_full.Rds
##     instead of output/04_merge/merged_processed.Rds. For IB1 that is 2,540,021
##     cells rather than 400,111 -- a 6.3x larger pool, because none of the
##     step-02/03 QC filtering has been applied. Every cell the vendor segmented
##     can now carry a light-chain call and a neighbourhood.
##     The trade-off is real and is NOT hidden: the pool now includes low-count and
##     QC-failing cells, so `qcCellsPassed` is carried through to the cells csv and
##     `artifact_suspect` (pct_cluster_i >= 50 | med_nCount < 300) is unchanged and
##     still worth filtering on. Compare against output/04_merge/ to see what the
##     QC filter was doing.
##   * OUTPUT goes to output/07_lightchain/ so the step-04 tables are preserved
##     for that comparison.
##   * sample is a command-line argument, and the sample column is
##     Run_Tissue_name (the step-1 object has no merge_sample).
##   * cell ids come from the meta.data ROWNAMES. The step-1 object also has a
##     vendor `cell_id` column, which is a different identifier -- the original
##     script's as.data.table(keep.rownames = "cell_id") would silently collide
##     with it here, so the vendor column is dropped first.
##   * only the needed meta.data columns are pulled before the data.table
##     conversion; converting all 102 columns of 2.54M cells is gratuitous.
##
## Usage: Rscript scripts/07_script_lightchain_neighborhoods.R [sample_id]
##   sbatch --mem=200G --cpus-per-task=2 --time=4:00:00 \
##     scripts/sbatch/run_rscript.sh scripts/07_script_lightchain_neighborhoods.R IB1
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(RANN); library(igraph)
  library(ggplot2); library(data.table); setDTthreads(1)
})

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"

IN_RDS <- sprintf("output/01_import_v1/seuratObject_%s_BPCells_full.Rds", SAMPLE)
## Keyed BY SAMPLE. The table filenames carry only the radius, so a second sample
## written to a shared directory would silently overwrite the first.
OUT    <- file.path("output/07_lightchain", SAMPLE); FIG <- file.path(OUT, "figures")
if (!file.exists(IN_RDS))
  stop("missing ", IN_RDS, " -- run scripts/01_import_cosmx_seurat.R first", call. = FALSE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

SAMPLE_COL <- "Run_Tissue_name"
KAPPA   <- "IGKC"; LAMBDA <- "IGLC1/2"
## CELL CALLING. The original MIN_LC = 2 with "whichever chain is larger wins" was
## far too permissive: measured on the pre-QC pool, 89,693 of 174,007 admitted cells
## had exactly 2 light-chain counts, 76,565 cells were a bare 2/0, and 85,532 of the
## 158,918 kappa/lambda calls (53.8%) rested on a dominant count of <= 2, with 10,396
## decided by a single-transcript margin. A neighbourhood's kappa fraction cannot
## mean much when half its votes are that thin.
## Now a cell enters the pool on total counts, and is only CALLED kappa/lambda if the
## dominant chain clears MIN_DOM and beats the other by MIN_MARGIN. Cells in the pool
## that fail the dominance test stay as "tie": they still anchor the spatial graph and
## count in n_cells, but they are excluded from n_called and so from kappa_frac.
MIN_LC     <- 4       # total light-chain counts to enter the pool at all
MIN_DOM    <- 3       # dominant chain must reach this many counts to be called
MIN_MARGIN <- 2       # ...and must beat the other chain by at least this much
# RADIUS SCAN. A single radius cannot be right: at 30 um the graph fragments
# (median community 1 cell) but small clonal foci survive; at 60 um connectivity
# reaches 99.6% yet the three kappa-restricted foci found at 30 um are ABSORBED
# into larger polytypic neighbourhoods and vanish. So scan, and report which foci
# persist across scales -- a focus seen at several radii is real, one seen at a
# single radius is a scale artefact. Outputs are per-radius; nothing overwrites.
RADII   <- c(0.02, 0.03, 0.04, 0.06)
K_NN    <- 40         # candidate neighbours before the radius filter
MIN_CELLS <- 20       # smallest neighbourhood we will score (was 30; see PLASMA note)

## PLASMA CELLS vs B CELLS. The pool contains two different cell types. Plasma cells
## carry light chain at 5-14% of their whole transcriptome (25-273 counts, essentially
## pure single chain, cytoplasmic); B cells carry surface BCR at ~1% (4-12 counts).
## Every plasma cell is monotypic by allelic/isotype exclusion, so a mixture of
## very-high kappa and very-high lambda cells is POLYTYPIC = reactive, not clonal --
## measured pooled kappa fraction of these cells is 0.720, right at the ~2:1 norm.
##
## They also contaminate their neighbours, chain-specifically. Measured on the r20um
## pool, for non-plasma called cells the call is dictated by the NEAREST plasma cell's
## chain and decays to background over one cell diameter:
##       <10um   92.8% kappa if that cell is kappa, 11.4% if it is lambda
##    10-20um   77.2% / 36.8%      20-40um  64.5% / 63.6%      none<40um  68.9%
## Symmetric, monotonic, chain-specific: segmentation spillover from an abundant
## cytoplasmic source. Since plasma cells cluster, it manufactures locally
## kappa-homogeneous patches -- exactly what this script is built to find.
##
## So plasma cells AND their halo are dropped before the graph is built, making
## kappa_frac a B-cell statistic by construction. MIN_CELLS drops 30 -> 20 because
## excluding ~31% of the pool made 30 a stricter bar than when it was chosen; scored
## counts are reported at both so the two changes stay separable.
##
## NOTE the plasma definition below is deliberately circular -- it types a cell by the
## same quantity whose bias we are removing. scripts/07_script_plasma_marker_check.R is
## the mitigation, not an afterthought: it ties this definition to the 11 independent
## plasma-cell markers in the panel. If that agreement fails, this exclusion is unjustified.
PC_MIN_PCT   <- 4     # light chain as % of the cell's nCount_RNA
PC_MIN_TOTAL <- 25    # ...and absolute light-chain counts
HALO_UM      <- 15    # ~one cell diameter, from the decay table above

## Independent plasma markers in the panel (only SSR4 is absent), pulled in the same
## counts read as the light chains and exported for the validation script.
PC_MARKERS <- c("MZB1","XBP1","DERL3","JCHAIN","PRDM1","SDC1","CD38","IRF4",
                "SEC11C","TNFRSF17","POU2AF1")
NULL_K  <- 0.62       # polytypic kappa fraction (normal ~2:1 kappa:lambda)

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
## Value scales run blue -> red (CLAUDE.md): kappa is the high end, so kappa is red.
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
          legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left", legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF);cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

## ---- 1. light-chain-high cells --------------------------------------------
cat("reading ", IN_RDS, "\n", sep="")
obj <- readRDS(IN_RDS)
rn  <- rownames(obj@meta.data)          # "c_2_100_1002" -- the id every join uses

NEED <- c(SAMPLE_COL, "x_slide_mm", "y_slide_mm", "inSitu_clusters",
          "nCount_RNA", "qcCellsPassed")
miss <- setdiff(NEED, colnames(obj@meta.data))
if (length(miss)) stop("meta.data is missing: ", paste(miss, collapse=", "), call.=FALSE)
md <- as.data.table(obj@meta.data[, NEED, drop = FALSE])
md[, cell_id := rn]                     # NOT the vendor cell_id column, see header

cnt <- LayerData(obj[["RNA"]], layer="counts")
stopifnot(all(c(KAPPA,LAMBDA) %in% rownames(cnt)))
stopifnot(identical(colnames(cnt), rn))
PC_HAVE <- intersect(PC_MARKERS, rownames(cnt))
if (length(PC_HAVE) < length(PC_MARKERS))
  message("plasma markers absent from the panel: ",
          paste(setdiff(PC_MARKERS, PC_HAVE), collapse=", "))
## one materialisation covering the light chains and the marker panel together, so the
## validation script never has to re-read the 1.9 GB object
lc <- as(cnt[c(KAPPA, LAMBDA, PC_HAVE), ], "dgCMatrix")
rm(obj,cnt); invisible(gc())
md[, `:=`(k_ct = as.vector(lc[KAPPA, ]), l_ct = as.vector(lc[LAMBDA, ]))]
md[, lc_total := k_ct + l_ct]
pc_mat <- Matrix::t(lc[PC_HAVE, , drop = FALSE])
md[, pc_n     := as.vector(Matrix::rowSums(pc_mat > 0))]
md[, pc_score := as.vector(Matrix::rowSums(pc_mat))]

d <- md[get(SAMPLE_COL) == SAMPLE]
cat("=== ", SAMPLE, ": ", nrow(d), " cells (pre-QC-filter pool) ===\n", sep="")
cat(sprintf("of which qcCellsPassed: %d (%.1f%%)\n",
            sum(d$qcCellsPassed %in% TRUE), 100*mean(d$qcCellsPassed %in% TRUE)))
cat("light-chain total counts per cell, quantiles:\n")
print(quantile(d$lc_total, c(.5,.75,.9,.95,.99,1)))
cat("\ncells by light-chain total:\n")
print(d[, .N, by=.(lc_total = pmin(lc_total, 5L))][order(lc_total)])

hi <- d[lc_total >= MIN_LC]
cat(sprintf("\nlight-chain-high cells (total >= %d): %d (%.1f%% of %s)\n",
            MIN_LC, nrow(hi), 100*nrow(hi)/nrow(d), SAMPLE))
cat(sprintf("of those, qcCellsPassed: %d (%.1f%%) -- the rest are what step 04 never saw\n",
            sum(hi$qcCellsPassed %in% TRUE), 100*mean(hi$qcCellsPassed %in% TRUE)))
hi[, dom := pmax(k_ct, l_ct)]
hi[, oth := pmin(k_ct, l_ct)]
hi[, called := dom >= MIN_DOM & (dom - oth) >= MIN_MARGIN]
hi[, call := ifelse(!called, "tie", ifelse(k_ct > l_ct, "kappa", "lambda"))]
cat(sprintf("dominance filter (dom >= %d and margin >= %d): %d of %d pool cells called, %d held as tie\n",
            MIN_DOM, MIN_MARGIN, sum(hi$called), nrow(hi), sum(!hi$called)))
cat("their per-cell calls:\n"); print(hi[, .N, by=call][order(-N)])
cat("dominant-count distribution of CALLED cells:\n")
print(hi[called == TRUE, .N, by = .(dom = pmin(dom, 8L))][order(dom)])
cat(sprintf("overall kappa fraction (excluding ties): %.3f\n",
            hi[call != "tie", mean(call == "kappa")]))

## ---- 1b. plasma / halo / B_clean ------------------------------------------
hi[, pct_lc := 100 * lc_total / nCount_RNA]
hi[, plasma := pct_lc >= PC_MIN_PCT & lc_total >= PC_MIN_TOTAL]
if (!any(hi$plasma)) stop("no plasma-like cells found -- check PC_MIN_PCT/PC_MIN_TOTAL",
                          call. = FALSE)
nnp <- RANN::nn2(as.matrix(hi[plasma == TRUE, .(x_slide_mm, y_slide_mm)]),
                 as.matrix(hi[, .(x_slide_mm, y_slide_mm)]), k = 1)
hi[, d_plasma_um := nnp$nn.dists[, 1] * 1000]
hi[, tier := ifelse(plasma, "plasma",
             ifelse(d_plasma_um < HALO_UM, "halo", "B_clean"))]
cat(sprintf("\nplasma (>=%g%% of RNA and >=%d counts): %d | halo (<%g um of one): %d | B_clean: %d\n",
            PC_MIN_PCT, PC_MIN_TOTAL, sum(hi$tier=="plasma"),
            HALO_UM, sum(hi$tier=="halo"), sum(hi$tier=="B_clean")))
cat(sprintf("light chain as %% of RNA -- plasma median %.1f%%, B_clean median %.1f%%\n",
            median(hi[tier=="plasma"]$pct_lc), median(hi[tier=="B_clean"]$pct_lc)))
cat(sprintf("kappa fraction -- plasma %.3f (polytypic ~%.2f expected), B_clean %.3f\n",
            hi[tier=="plasma"  & call!="tie", mean(call=="kappa")], NULL_K,
            hi[tier=="B_clean" & call!="tie", mean(call=="kappa")]))
cat("plasma cells by inSitu_clusters (expect one or two dominant):\n")
print(hi[tier=="plasma", .N, by=inSitu_clusters][order(-N)][1:5])

## export for scripts/07_script_plasma_marker_check.R -- the independent gate on the
## share-based definition above
## per-marker counts are materialised for POOL cells only -- dense over all 2.54M
## cells would be gratuitous, and only the pool is ever validated
pool_idx   <- match(hi$cell_id, md$cell_id)
stopifnot(!anyNA(pool_idx))
mk_markers <- as.data.table(as.matrix(pc_mat[pool_idx, , drop = FALSE]))
setnames(mk_markers, PC_HAVE)
mk_markers[, cell_id := hi$cell_id]
mk <- merge(hi[, .(cell_id, k_ct, l_ct, lc_total, pct_lc, nCount_RNA, call, tier,
                   d_plasma_um, inSitu_clusters, qcCellsPassed)],
            md[, .(cell_id, pc_n, pc_score)], by = "cell_id")
mk <- merge(mk, mk_markers, by = "cell_id")
fwrite(mk, file.path(OUT, "plasma_markers.csv.gz"))
cat("wrote ", file.path(OUT, "plasma_markers.csv.gz"), "\n", sep="")
rm(mk); invisible(gc())

## The graph is built on B_clean only. Plasma and halo cells stay in `hi` so they can
## still be written out and drawn, but they get no neighbourhood.
hi_all <- copy(hi)
hi     <- hi[tier == "B_clean"]
cat(sprintf("graph pool after dropping plasma+halo: %d of %d cells (%.0f%% dropped)\n",
            nrow(hi), nrow(hi_all), 100*(1 - nrow(hi)/nrow(hi_all))))

## ---- 2-4. per-radius: graph -> neighbourhoods -> scores --------------------
ALL <- list()
for (RADIUS in RADII) {
cat(sprintf('\n================ RADIUS %.0f um ================\n', RADIUS*1000))
xy <- as.matrix(hi[, .(x_slide_mm, y_slide_mm)])
nn <- RANN::nn2(xy, k = min(K_NN + 1, nrow(xy)))
i <- rep(seq_len(nrow(xy)), ncol(nn$nn.idx) - 1)
j <- as.vector(nn$nn.idx[, -1]); dist <- as.vector(nn$nn.dists[, -1])
keep <- dist <= RADIUS & i != j
cat(sprintf("\ngraph: %d B_clean cells, %d edges within %.0f um (mean degree %.1f)\n",
            nrow(xy), sum(keep), RADIUS*1000, 2*sum(keep)/nrow(xy)))
gr <- igraph::graph_from_edgelist(cbind(i[keep], j[keep]), directed = FALSE)
gr <- igraph::simplify(igraph::add_vertices(gr, max(0, nrow(xy) - igraph::vcount(gr))))
deg <- igraph::degree(gr)
cat(sprintf("connectivity: %.1f%% of B_clean cells have >=1 neighbour, %.1f%% have >=5\n",
            100*mean(deg >= 1), 100*mean(deg >= 5)))
rm(nn, i, j, dist, keep); invisible(gc())

## ---- 3. neighbourhoods ----------------------------------------------------
set.seed(1)
cl <- igraph::cluster_louvain(gr)
hi[, nb := as.integer(igraph::membership(cl))]
hi[, cc := as.integer(igraph::components(gr)$membership)]
cat(sprintf("neighbourhoods: %d Louvain communities, %d connected components\n",
            length(unique(hi$nb)), length(unique(hi$cc))))
cat("community size quantiles:\n"); print(quantile(hi[, .N, by=nb]$N, c(.5,.75,.9,.95,1)))

## ---- 4. score each neighbourhood ------------------------------------------
nbs <- hi[, .(n_cells = .N,
              n_kappa = sum(call=="kappa"), n_lambda = sum(call=="lambda"),
              n_tie = sum(call=="tie"),
              lc_counts = sum(lc_total), med_lc = as.numeric(median(lc_total)),
              x = mean(x_slide_mm), y = mean(y_slide_mm),
              span_mm = round(max(sqrt((x_slide_mm-mean(x_slide_mm))^2 +
                                       (y_slide_mm-mean(y_slide_mm))^2)), 3),
              pct_cluster_i = round(100*mean(inSitu_clusters=="i"), 1),
              pct_qc_passed = round(100*mean(qcCellsPassed %in% TRUE), 1),
              med_nCount = as.numeric(median(nCount_RNA))), by = nb]
nbs[, n_called := n_kappa + n_lambda]

## How close is a neighbourhood to the excluded plasma compartment? Plasma cells have
## no nb of their own now, so this is the only way a plasma-rich region stays visible.
ctx <- hi[, .(n_plasma_near = sum(d_plasma_um < HALO_UM * 2),
              med_d_plasma_um = round(median(d_plasma_um), 1)), by = nb]
nbs <- merge(nbs, ctx, by = "nb")

## Report at BOTH size floors so the MIN_CELLS change stays separable from the effect
## of dropping plasma cells.
cat(sprintf("scored at MIN_CELLS=20: %d | at 30: %d (of %d communities)\n",
            sum(nbs$n_called >= 20), sum(nbs$n_called >= 30), nrow(nbs)))
nbs <- nbs[n_called >= MIN_CELLS]
nbs[, kappa_frac := round(n_kappa/n_called, 3)]
nbs[, p := mapply(function(a,b) binom.test(a,b,p=NULL_K)$p.value, n_kappa, n_called)]
nbs[, q := p.adjust(p, "BH")]
nbs[, class := ifelse(q < 0.05 & kappa_frac >= 0.85, "kappa-restricted (FL-like)",
              ifelse(q < 0.05 & kappa_frac <= 0.35, "lambda-skewed",
              ifelse(q >= 0.05, "polytypic (reactive-like)", "intermediate")))]
nbs[, artifact_suspect := pct_cluster_i >= 50 | med_nCount < 300]
setorder(nbs, -kappa_frac)
cat(sprintf("\nscored neighbourhoods (>= %d called cells): %d\n", MIN_CELLS, nrow(nbs)))
cat("\n=== neighbourhood classes ===\n")
print(nbs[, .(n = .N, cells = sum(n_cells), med_kappa = median(kappa_frac)),
          by = .(class, artifact_suspect)][order(class, artifact_suspect)])
cat("\n=== top kappa-restricted neighbourhoods (clean) ===\n")
print(head(nbs[artifact_suspect == FALSE & class == "kappa-restricted (FL-like)",
  .(nb, n_cells, n_kappa, n_lambda, kappa_frac, q = signif(q,3), span_mm, x = round(x,2),
    y = round(y,2), pct_cluster_i, pct_qc_passed)], 15))
cat("\n=== polytypic neighbourhoods (clean) ===\n")
print(head(nbs[artifact_suspect == FALSE & class == "polytypic (reactive-like)",
  .(nb, n_cells, n_kappa, n_lambda, kappa_frac, q = signif(q,3), span_mm, x = round(x,2),
    y = round(y,2), pct_cluster_i, pct_qc_passed)], 15))
nbs[, radius_um := RADIUS*1000]
ALL[[as.character(RADIUS)]] <- copy(nbs)
sc(sprintf("lightchain_neighborhoods_r%02.0fum.csv", RADIUS*1000), nbs)
## Every pool cell, not just the graphed ones: the views need to draw plasma and halo
## cells as their own tiers rather than showing them as plain background. nb/cc are
## populated for B_clean cells only.
out_cells <- merge(hi_all[, .(cell_id, k_ct, l_ct, lc_total, pct_lc, dom, oth, called,
                              call, tier, d_plasma_um, x_slide_mm, y_slide_mm,
                              inSitu_clusters, nCount_RNA, qcCellsPassed)],
                   hi[, .(cell_id, nb, cc)], by = "cell_id", all.x = TRUE)
sc(sprintf("lightchain_cells_r%02.0fum.csv", RADIUS*1000), out_cells)

## ---- figures ---------------------------------------------------------------
p1 <- ggplot(nbs, aes(kappa_frac)) +
  geom_histogram(binwidth = 0.05, fill = DIV_HI, colour = SURF, linewidth = 0.3) +
  geom_vline(xintercept = NULL_K, colour = INK, linetype = "22", linewidth = 0.4) +
  scale_x_continuous(limits = c(0,1)) +
  labs(title = sprintf("Light-chain ratio per neighbourhood (%s)", SAMPLE),
       subtitle = sprintf(paste0("%d spatial neighbourhoods of light-chain-high cells (>=%d called ",
         "cells each), pre-QC-filter pool.\nDashed line = polytypic expectation (%.2f). A clonal FL ",
         "neighbourhood sits near 1.0."), nrow(nbs), MIN_CELLS, NULL_K),
       x = "kappa fraction", y = "neighbourhoods",
       caption = "Data: output/07_lightchain/lightchain_neighborhoods_r*um.csv") + theme_viz()
sf(sprintf("fig_07_lc_nb_ratio_hist_r%02.0fum.png", RADIUS*1000), p1, 8, 4.4)

p2 <- ggplot(nbs, aes(x, y, fill = kappa_frac)) +
  geom_point(shape = 21, aes(size = n_cells), colour = INK, stroke = 0.3) +
  scale_fill_gradient2(low = DIV_LO, mid = DIV_MID, high = DIV_HI, midpoint = NULL_K,
                       limits = c(0,1), name = "kappa fraction") +
  scale_size_area(max_size = 9, name = "cells") +
  coord_equal() +
  labs(title = sprintf("Where the light-chain neighbourhoods sit in %s", SAMPLE),
       subtitle = "Red = kappa-restricted (FL-like), grey = polytypic, blue = lambda-skewed.",
       x = "x (mm)", y = "y (mm)",
       caption = "Data: output/07_lightchain/lightchain_neighborhoods_r*um.csv") +
  theme_viz() + theme(panel.grid.major.y = element_blank())
sf(sprintf("fig_07_lc_nb_map_r%02.0fum.png", RADIUS*1000), p2, 7.6, 8)

p3 <- ggplot(hi[nb %in% nbs$nb], aes(x_slide_mm, y_slide_mm, colour = call)) +
  geom_point(size = 0.25) +
  scale_colour_manual(values = c(kappa = DIV_HI, lambda = DIV_LO, tie = MUTED), name = NULL) +
  coord_equal() +
  labs(title = sprintf("Light-chain-high cells in %s, coloured by call", SAMPLE),
       subtitle = sprintf("Cells with >=%d light-chain counts, in scored neighbourhoods.", MIN_LC),
       x = "x (mm)", y = "y (mm)",
       caption = "Data: output/07_lightchain/lightchain_cells_r*um.csv") +
  theme_viz() + theme(panel.grid.major.y = element_blank())
sf(sprintf("fig_07_lc_cells_map_r%02.0fum.png", RADIUS*1000), p3, 7.6, 8)
} # end radius loop
sc("lightchain_call_params.csv",
   data.table(sample = SAMPLE, min_lc = MIN_LC, min_dom = MIN_DOM,
              min_margin = MIN_MARGIN, min_cells = MIN_CELLS, null_k = NULL_K,
              pc_min_pct = PC_MIN_PCT, pc_min_total = PC_MIN_TOTAL, halo_um = HALO_UM,
              pool_cells = nrow(hi_all), called_cells = sum(hi_all$called),
              n_plasma = sum(hi_all$tier == "plasma"),
              n_halo = sum(hi_all$tier == "halo"),
              n_bclean = sum(hi_all$tier == "B_clean"),
              graph_cells = nrow(hi)))

## ---- cross-radius comparison ----------------------------------------------
A <- rbindlist(ALL, use.names = TRUE)
sc("lightchain_neighborhoods_all_radii.csv", A)
cat("\n================ CROSS-RADIUS SUMMARY ================\n")
print(A[, .(scored = .N, flagged = sum(artifact_suspect),
            kappa_restricted = sum(class == "kappa-restricted (FL-like)" & !artifact_suspect),
            lambda_skewed = sum(class == "lambda-skewed" & !artifact_suspect),
            polytypic = sum(class == "polytypic (reactive-like)" & !artifact_suspect),
            med_size = as.numeric(median(n_cells)), max_size = max(n_cells),
            med_kappa = round(median(kappa_frac), 3)), by = radius_um])

# do restricted foci recur at the same place across radii? (0.15 mm grid)
R <- A[class == "kappa-restricted (FL-like)" & !artifact_suspect]
if (nrow(R)) {
  R[, `:=`(gx = round(x/0.15), gy = round(y/0.15))]
  per <- R[, .(radii = paste(sort(unique(radius_um)), collapse=","), n_radii = uniqueN(radius_um),
               x = round(mean(x),2), y = round(mean(y),2),
               kappa = round(mean(kappa_frac),3), cells = max(n_cells)), by=.(gx,gy)][order(-n_radii)]
  cat("\n=== kappa-restricted foci, persistence across radii ===\n"); print(per)
  cat("(a focus found at multiple radii is robust; one radius only = scale artefact)\n")
  sc("lightchain_restricted_foci_persistence.csv", per)
} else cat("\nno clean kappa-restricted neighbourhoods at any radius\n")
cat("\nOK\n")
