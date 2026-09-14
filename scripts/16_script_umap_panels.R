## ============================================================================
## 16_script_umap_panels.R -- the four cell-type calls, side by side in UMAP space.
##
## Four panels, one per way of naming a cell in this project, on the same embedding:
##   1. inSitu_clusters          the vendor's InSituType call. `i` is the over-fixation defect.
##   2. pearson_clusters_low     step 03's scPearsonPCA clustering, unlabelled
##   3. audited label            step 09b's marker audit of those clusters (via lib/cluster_labels.R)
##   4. celltype                 step 16b's per-cell marker argmax (the section-10.1 call)
##
## THREE FILTERINGS, and the point of the series is the comparison between them:
##   none      every cell in the object
##   assigned  only cells the step-16b argmax actually resolved (sentinels dropped)
##   nodefect  only cells outside InSituType cluster `i`
##
## and two embeddings: the ORIGINAL `pearsonumap` (so a filter is a subset of a fixed map), or
## a RECOMPUTED one (so the filter changes the embedding itself). The recomputation replicates
## step 03's block exactly -- gene_frequency -> 900 HVGs -> sparse_quasipoisson_pca_seurat ->
## make_umap (03_script_filter_and_process.R:244-263). It does NOT re-run FindClusters:
## panels 2 and 3 must keep the cluster ids that step 09b audited, or panel 3 becomes a join
## against ids that no longer mean anything. So the labels are carried, and only the embedding
## moves -- which is the comparison being asked for.
##
## PANEL 3 HAS A SAMPLE-SPECIFIC CAVEAT, surfaced on the figure rather than buried. Step 09b
## scores each cluster by z ACROSS clusters within the sample, so when one lineage dominates a
## sample its markers stop being distinctive between clusters and the z collapses. RB1 is ~60%
## B cells and its largest cluster (836,909 cells) audits as "low-information" at only 3%
## defect -- yet its own table records best_supported = B cell on BANK1/FCER2/CD19/CD79A/CD79B
## with a 0.75 margin over T. So where the status is not confident, the panel shows the
## best-supported family in brackets. Painting 60% of RB1 as uninformative would be worse than
## the small dishonesty of a bracket.
##
## Usage: Rscript scripts/16_script_umap_panels.R <IB1|RB1|merged> <none|assigned|nodefect> <0|1> [replot]
##
## A 4th argument `replot` redraws the figure from the per-cell CSV this script already writes
## (cell_id, U1, U2 and all four labels), without touching the Seurat object. That matters
## because the recomputed embeddings cost hours of quasipoisson PCA: a change to the palette or
## the point size must never require redoing them.
## ============================================================================

if (!file.exists("scripts/lib/cluster_labels.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/cluster_labels.R")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(patchwork)
  library(ggrepel); library(data.table); setDTthreads(1)
})

args  <- commandArgs(trailingOnly = TRUE)
OBJ   <- if (length(args) >= 1) args[1] else "IB1"
FILT  <- if (length(args) >= 2) args[2] else "none"
REPRO <- if (length(args) >= 3) as.integer(args[3]) else 0L
REPLOT<- length(args) >= 4 && args[4] %in% c("replot","1","TRUE")
if (!OBJ  %in% c("IB1","RB1","merged")) stop("object must be IB1, RB1 or merged", call.=FALSE)
if (!FILT %in% c("none","assigned","nodefect")) stop("filter must be none, assigned or nodefect", call.=FALSE)
if (REPRO == 1L && FILT == "none")
  stop("reprocessing an unfiltered object just reproduces step 03; not run", call. = FALSE)

IN_RDS <- if (OBJ == "merged") "output/04_merge/merged_processed.Rds" else
          sprintf("output/03_filter_and_process/%s_02_processed.Rds", OBJ)
CT_F   <- "output/16_celltypes/cell_celltypes.csv.gz"
OUT    <- "output/16_umaps"
TAG    <- sprintf("%s_%s_%s", OBJ, FILT, ifelse(REPRO == 1L, "reumap", "origumap"))
for (f in c(IN_RDS, CT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

N_DRAW <- 250000L      # points drawn per panel; the full object is far past overplotting
LAB_MIN_N <- 200L      # a group needs this many DRAWN points to earn an in-plot label; below
                       # that the median position is not a location, it is a coin flip
SENT   <- c("ambiguous","below_background","all_markers_0")
SEED   <- 1L
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GREY<-"#dedcd4"
th <- theme_minimal(base_size = 10) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA), panel.grid=element_blank(),
        axis.text=element_blank(), axis.ticks=element_blank(),
        axis.title=element_text(colour=INK2, size=8),
        plot.title=element_text(colour=INK, face="bold", size=10),
        plot.subtitle=element_text(colour=INK2, size=7.5),
        legend.text=element_text(colour=INK2, size=6.5), legend.title=element_blank(),
        legend.position="right", legend.key.height=unit(7,"pt"))

cat("============================================================\n")
cat("step 16 UMAPs -- ", TAG, "\n", sep=""); cat("============================================================\n")

CSV_F <- file.path(OUT, sprintf("umap_cells_%s.csv.gz", TAG))
if (REPLOT) {
  if (!file.exists(CSV_F)) stop("replot needs ", CSV_F, " -- run without `replot` first",
                                call. = FALSE)
  cat("REPLOT: redrawing from ", CSV_F, " -- the object and the PCA are not touched\n", sep="")
  d <- fread(CSV_F)
  n0 <- nrow(d)
  RED <- if (REPRO == 1L) "recomputed pearsonumap" else "original pearsonumap (step 03/04)"
} else {

## ---- 1. object and labels ----------------------------------------------------
cat("reading ", IN_RDS, "\n", sep=""); seu <- readRDS(IN_RDS)
## The vendor `cell_id` / `cell_ID` META COLUMNS are a DIFFERENT identifier from the rownames,
## which are the canonical ids (c_<slide>_<fov>_<cell>) everything else in this project joins
## on. Dropping them before keep.rownames is what step 07 does; leaving them in makes
## as.data.table() emit a duplicated `cell_id` and the join below dies.
mdf <- seu@meta.data
for (bad in c("cell_id","cell_ID")) if (bad %in% colnames(mdf)) mdf[[bad]] <- NULL
md <- as.data.table(mdf, keep.rownames = "cell_id")
cat(sprintf("cells %d | reductions %s\n", ncol(seu), paste(names(seu@reductions), collapse=",")))

ct <- fread(CT_F, select = c("cell_id","celltype","lineage16"))
n0 <- nrow(md); md <- merge(md, ct, by = "cell_id", all.x = TRUE, sort = FALSE)
if (nrow(md) != n0) stop("celltype join changed the row count", call. = FALSE)
cat(sprintf("step-16b celltype present for %d of %d cells (%.1f%%)\n",
            md[, sum(!is.na(celltype))], nrow(md), 100*md[, mean(!is.na(celltype))]))
md[is.na(celltype), celltype := "not in step 16"]

## Panel 3: audited labels. For a single sample the join is cluster_id -> label. For the merged
## object the clusters are JOINT and step 09b never audited them, so the per-sample audit is
## mapped through `sample_pearson_clusters_low` ("IB1_3"), which is the only honest route.
lab_of <- function(S) {
  L <- tryCatch(cluster_labels(S, "pearson_clusters_low"), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  L[, disp := ifelse(status == "confident", label,
              ifelse(is.na(best_supported) | best_supported == "", label,
                     sprintf("%s [%s?]", label, best_supported)))]
  L[, .(cluster_id = as.character(cluster_id), disp)]
}
if (OBJ == "merged") {
  if (!"sample_pearson_clusters_low" %in% names(md))
    stop("merged object lacks sample_pearson_clusters_low", call. = FALSE)
  md[, .spl := as.character(sample_pearson_clusters_low)]
  md[, .samp := sub("_.*$", "", .spl)][, .cl := sub("^[^_]*_", "", .spl)]
  LL <- rbindlist(lapply(c("IB1","RB1"), function(S) {
    x <- lab_of(S); if (is.null(x)) return(NULL); x[, .samp := S][]
  }), fill = TRUE)
  md <- merge(md, LL, by.x = c(".samp",".cl"), by.y = c(".samp","cluster_id"),
              all.x = TRUE, sort = FALSE)
  md[, audited := ifelse(is.na(disp), "not audited", disp)]
} else {
  x <- lab_of(OBJ)
  if (is.null(x)) { md[, audited := "not audited (09b not run)"] } else {
    md[, .cl := as.character(pearson_clusters_low)]
    md <- merge(md, x, by.x = ".cl", by.y = "cluster_id", all.x = TRUE, sort = FALSE)
    md[, audited := ifelse(is.na(disp), "not audited", disp)]
  }
}
cat("\npanel 3 (audited) levels:\n"); print(md[, .N, by = audited][order(-N)])

## ---- 2. filter ---------------------------------------------------------------
keep <- switch(FILT,
  none     = rep(TRUE, nrow(md)),
  assigned = !md$celltype %in% c(SENT, "not in step 16"),
  nodefect = md$inSitu_clusters != "i")
cat(sprintf("\nfilter '%s': %d of %d cells kept (%.1f%%)\n",
            FILT, sum(keep), nrow(md), 100*mean(keep)))
if (sum(keep) < 5000) stop("filter left too few cells", call. = FALSE)

## ---- 3. embedding ------------------------------------------------------------
if (REPRO == 1L) {
  ## Replicates 03_script_filter_and_process.R:244-263 on the filtered cells.
  suppressPackageStartupMessages(library(scPearsonPCA))
  cat("\nsubsetting and recomputing the pearson UMAP on the filtered cells\n")
  seu <- subset(seu, cells = md$cell_id[keep])
  md  <- md[keep]
  ca  <- if (inherits(seu[["RNA"]], "Assay5")) seu[["RNA"]]$counts else seu[["RNA"]]@counts
  tc  <- Matrix::colSums(ca)
  gf  <- scPearsonPCA::gene_frequency(ca)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 900)
  hv  <- VariableFeatures(seu)
  ch  <- as(ca[hv, ], "dgCMatrix")
  cat(sprintf("  %d HVGs x %d cells, %.1f%% nonzero\n", nrow(ch), ncol(ch),
              100*length(ch@x)/(as.numeric(nrow(ch))*ncol(ch))))
  po  <- sparse_quasipoisson_pca_seurat(ch, totalcounts = tc, grate = gf[hv],
                                        scale.max = 10, do.scale = TRUE, do.center = TRUE)
  uo  <- scPearsonPCA::make_umap(po)
  EMB <- as.data.table(uo$ump@cell.embeddings, keep.rownames = "cell_id")
  RED <- "recomputed pearsonumap"
} else {
  if (!"pearsonumap" %in% names(seu@reductions))
    stop("no pearsonumap reduction in ", IN_RDS, call. = FALSE)
  EMB <- as.data.table(Embeddings(seu, "pearsonumap"), keep.rownames = "cell_id")
  md  <- md[keep]
  RED <- "original pearsonumap (step 03/04)"
}
setnames(EMB, 2:3, c("U1","U2"))
d <- merge(md, EMB, by = "cell_id", sort = FALSE)
if (!nrow(d)) stop("no cells survived the embedding join", call. = FALSE)
}   # end of the non-replot branch
cat(sprintf("plotting %d cells on the %s\n", nrow(d), RED))

## ---- 4. panels ---------------------------------------------------------------
set.seed(SEED)
pl <- d[sample.int(.N, min(.N, N_DRAW))]
## Three things this has to get right, each of which the first version got wrong.
##
## 1. A HIGHLIGHT MUST ACTUALLY STAND OUT. Colouring cluster `i` red while the other nine
##    InSituType levels came from the default hue palette put `a` in almost the same red, so
##    the one level the figure exists to show was indistinguishable from a neighbour. When a
##    highlight is requested, everything else is desaturated to blue-greys.
## 2. SENTINELS MUST RECEDE, not dominate. Unassigned cells are ~50% of IB1, and drawn at the
##    same size and alpha as real calls they simply grey out the panel. They get their own
##    layer, smaller and more transparent, underneath.
## 3. LEGEND ORDER: real classes first, sentinels last, so the eye is not led into the greys.
## A highlight needs the OTHER levels to stay readable. The first attempt desaturated them to
## blue-greys, which made cluster `i` pop but turned the remaining nine InSituType levels into
## an indistinguishable mush -- trading one unreadable panel for another. This is a fixed
## qualitative palette that simply CONTAINS NO CRIMSON, so `i` is unique by hue rather than by
## everything else being drained of colour.
QUAL <- c("#0072B2","#009E73","#E69F00","#56B4E9","#CC79A7","#8B5A2B",
          "#7B3294","#66A61E","#01665E","#B8860B","#4B0082","#2F4F4F")
mkpal <- function(lv, sentinels = character(), hi = NULL) {
  real <- setdiff(lv, c(sentinels, hi))
  p <- setNames(rep(GREY, length(lv)), lv)
  if (length(real))
    p[real] <- if (is.null(hi)) scales::hue_pal(l = 55, c = 100)(length(real)) else
               rep(QUAL, length.out = length(real))
  if (!is.null(hi) && hi %in% lv) p[hi] <- "#c1121f"
  p
}
mk <- function(col, ttl, sent = character(), hi = NULL) {
  v   <- as.character(pl[[col]])
  real<- sort(setdiff(unique(v), sent))
  lv  <- c(real, sort(intersect(unique(v), sent)))   # sentinels last in the legend
  pal <- mkpal(lv, sent, hi)
  dd  <- copy(pl)[, .lab := factor(v, levels = lv)][, .sent := v %in% sent]
  ## In-plot labels at each group's MEDIAN position, in addition to the legend. The median is
  ## used rather than the mean because a UMAP group is routinely non-convex or split across
  ## arms, and a mean lands in empty space between them. Sentinel classes are NOT labelled:
  ## they are diffuse by definition, so a label for them marks nothing and only adds clutter.
  ## ggrepel with a white halo rather than a filled box -- a box large enough to read would
  ## occlude the points it is describing.
  labd <- dd[.sent == FALSE, .(U1 = median(U1), U2 = median(U2), n = .N), by = .lab]
  labd <- labd[n >= LAB_MIN_N]
  labd[, col := pal[as.character(.lab)]]
  ggplot() +
    ## sentinels underneath, smaller and fainter
    { if (any(dd$.sent)) geom_point(data = dd[.sent == TRUE], aes(U1, U2, colour = .lab),
                                    size = .20, stroke = 0, alpha = .30) } +
    geom_point(data = dd[.sent == FALSE], aes(U1, U2, colour = .lab),
               size = .38, stroke = 0, alpha = .85) +
    scale_colour_manual(values = pal, breaks = lv, drop = FALSE) +
    { if (nrow(labd)) ggrepel::geom_text_repel(
        data = labd, aes(U1, U2, label = .lab), colour = labd$col,
        size = 2.7, fontface = "bold", bg.color = "white", bg.r = .14,
        box.padding = .18, point.padding = 0, min.segment.length = .3,
        segment.colour = MUTED, segment.size = .25, max.overlaps = 40, seed = SEED,
        inherit.aes = FALSE) } +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1), ncol = 1)) +
    labs(x = "UMAP 1", y = "UMAP 2", title = ttl,
         subtitle = sprintf("%d levels%s", length(lv),
           if (any(dd$.sent)) sprintf("; %.0f%% unassigned (faint, unlabelled)",
                                      100*mean(dd$.sent)) else "")) + th
}
p1 <- mk("inSitu_clusters",      "1. InSituType clusters (vendor)", hi = "i")
p2 <- mk("pearson_clusters_low", "2. pearson_clusters_low (step 03)")
p3 <- mk("audited",              "3. audited label (step 09b)",
         sent = c("not audited","not audited (09b not run)"))
p4 <- mk("celltype",             "4. celltype, marker argmax (step 16b)",
         sent = c(SENT, "not in step 16"))

sub <- sprintf(paste0("%s | filter '%s' -> %d cells of %d (%.1f%%) | %s | %d of them drawn\n",
  "Panels 2 and 3 keep the ORIGINAL cluster ids even when the embedding is recomputed, ",
  "because step 09b audited those ids.\nRed in panel 1 is over-fixation cluster `i`. ",
  "Grey in panels 3-4 is unassigned. A bracketed [family?] in panel 3 is 09b's ",
  "best-supported\nfamily where its z threshold withheld a confident call -- see the script ",
  "header on why that matters for RB1."),
  OBJ, FILT, nrow(d), n0, 100*nrow(d)/n0, RED, nrow(pl))
pg <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = sprintf("%s — four ways of naming a cell, %s", OBJ,
                    ifelse(REPRO == 1L, "embedding recomputed after filtering",
                           "original embedding")),
                  subtitle = sub,
                  theme = theme(plot.title=element_text(face="bold", colour=INK, size=13),
                                plot.subtitle=element_text(colour=INK2, size=7.5),
                                plot.background=element_rect(fill=SURF, colour=NA)))
f <- file.path(OUT, sprintf("fig_16_umap_%s.png", TAG))
ggsave(f, pg, width = 15, height = 11, dpi = 165, bg = SURF, limitsize = FALSE)
cat("wrote ", f, "\n", sep="")
if (!REPLOT)
  fwrite(d[, .(cell_id, U1, U2, inSitu_clusters, pearson_clusters_low, audited, celltype)],
         CSV_F, compress = "gzip")
writeLines(format(Sys.time()), file.path(OUT, sprintf("%s.done", TAG)))
cat("============================================================\nOK\n")
