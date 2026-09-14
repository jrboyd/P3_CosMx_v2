## ============================================================================
## 14_script_lightchain_in_out.R -- what is IN a light-chain neighborhood?
##
## Step 07 built light-chain neighborhoods from IGKC/IGLC-called B and plasma cells, and
## everything since has described those neighborhoods by the CALLED cells only. This step
## characterises ALL cells occupying the same tissue, against matched tissue elsewhere:
##
##   1. cluster prevalence IN vs OUT, for pearson_clusters_low and _high
##   2. is composition consistent across neighborhoods, and do any stick out
##   3. expression IN vs OUT within each cluster -> two heatmaps
##
## THE CIRCULARITY, STATED UP FRONT. A neighborhood IS a cluster of light-chain-called
## B/plasma cells. B and plasma clusters are therefore enriched inside it BY
## CONSTRUCTION, and any "B cells are enriched" result here is a restatement of the
## definition, not a finding. Every prevalence number is reported twice: once over all
## cells, and once with the 14,641 CALLED MEMBER CELLS REMOVED. Only the second is
## evidence about anything. The same applies to IGKC/IGLC in the expression heatmaps --
## they are the selection variable and are marked as such.
##
## THE DENSITY TRAP, which cost step 08 a whole analysis. A neighborhood is dense by
## construction, so a fraction-of-neighbours statistic is biased before any biology.
## Here the unit is a REGION, and control regions are matched on CELL COUNT rather than
## on radius: for each real neighborhood footprint holding n cells, controls take the n
## nearest cells to an out-of-neighborhood seed. Same n, same compactness, different
## tissue.
##
## THE PSEUDO-REPLICATION TRAP, which cost step 11 a Tfh result. 400k cells inside ~154
## neighborhoods are not 400k independent observations. Every p-value here is computed
## with the REGION as the unit (154 neighborhoods vs 462 matched controls), never the
## cell. Cell-level counts are reported as effect sizes only.
##
## Usage: Rscript scripts/14_script_lightchain_in_out.R [SAMPLE]
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(RANN)
  library(ggplot2); library(data.table); setDTthreads(1)
})

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"

IN_RDS  <- sprintf("output/03_filter_and_process/%s_03_HCA_cell_typed.Rds", SAMPLE)
CELLS_F <- sprintf("output/07_lightchain/%s/lightchain_cells_r40um.csv", SAMPLE)
NBHD_F  <- sprintf("output/07_lightchain/%s/lightchain_neighborhoods_r40um.csv", SAMPLE)
## Labels come from the step-09b marker audit, never straight from HierAnnot -- four of
## nine pearson_clusters_low calls are wrong there, including stroma called macrophage.
source("scripts/lib/cluster_labels.R")
OUT     <- file.path("output/14_lightchain_in_out", SAMPLE)
FIG     <- file.path(OUT, "figures")

R_IN      <- 40      # um from a called member: the neighborhood footprint. Matches the
                     # radius step 07 used to LINK members, so the footprint is the
                     # region the neighborhood was actually defined over.
R_BUF     <- 100     # um: cells between R_IN and R_BUF are a boundary zone, excluded
                     # from both groups and reported. Without it, "out" is contaminated
                     # by tissue immediately adjacent to a neighborhood.
CTRL_PER  <- 3L      # matched control regions per real neighborhood
CAND_MULT <- 8L      # candidate control seeds per slot, to select on depth from
DEFECT_CL <- "i"     # the over-fixation cluster, dropped as in steps 07/11
MIN_CELLS_RC <- 5L   # min cells for a (region, cluster) pseudobulk to be usable. A
                     # synthetic run with 20 neighborhoods showed 10 silently dropped 7 of
                     # 9 clusters -- including the one carrying a planted 2.5x effect --
                     # because a rare cluster contributes few cells to any one footprint.
MIN_REG      <- 12L  # min regions per side before a cluster is testable at all
TOP_DE    <- 40L     # genes shown in the differential heatmap
SEED      <- 1L

## Canonical lymph-node markers, shown so the reader can check the step-09b labels
## against the evidence rather than taking them on trust. 09b scores a superset of this
## panel; this is the display subset.
MARKERS <- list(
  `B`           = c("MS4A1","CD19","CD79A","CD79B","PAX5","BANK1","FCER2","CR2"),
  `GC B`        = c("BCL6","AICDA","RGS13","LMO2","MEF2B"),
  `Plasma`      = c("MZB1","XBP1","PRDM1","DERL3","TNFRSF17","JCHAIN"),
  `Light chain` = c("IGKC","IGLC1","IGLC2"),
  `T`           = c("CD3D","CD3E","CD3G","TRAC","CD2","IL7R"),
  `CD4/CD8`     = c("CD4","CD8A","CD8B"),
  `Tfh`         = c("CXCL13","PDCD1","ICOS"),
  `Treg`        = c("FOXP3","IL2RA","CTLA4"),
  `NK`          = c("NKG7","GNLY","KLRD1","NCAM1"),
  `Myeloid`     = c("CD68","CD163","LYZ","ITGAX","CSF1R","C1QA"),
  `DC`          = c("LAMP3","BATF3","CLEC9A","CD1C","FCER1A"),
  `pDC`         = c("LILRA4","IRF7"),
  `Mast`        = c("TPSAB1","CPA3","KIT"),
  `Neutrophil`  = c("FCGR3B","CSF3R","S100A8"),
  `Stroma/FDC`  = c("COL1A1","PDGFRB","ACTA2","FDCSP"),
  `Endothelial` = c("PECAM1","VWF","CDH5"),
  `Cycling`     = c("MKI67","TOP2A"))
## The selection variable. Enrichment of these inside a neighborhood is definitional.
CIRCULAR <- c("IGKC","IGLC1","IGLC2")

dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#2a78d6"; DIV_LO<-"#eb6834"; ACC<-"#1baf7a"
## Separate names for the VALUE scale. DIV_HI/DIV_LO above are also used categorically
## (in-vs-out groups, neighbourhood class), and those are identities rather than
## magnitudes -- flipping them along with the scale would be meaningless. Value scales run
## blue -> red per CLAUDE.md: red is the high end, so on the heatmaps red is up in the
## in-vs-out comparison and blue is down.
SCALE_HI<-"#eb6834"; SCALE_LO<-"#2a78d6"
GRP_COLS <- c(`in` = DIV_HI, out = MUTED)
theme_viz <- function(base=11) theme_minimal(base_size=base) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour=GRID,linewidth=.25),
        axis.text=element_text(colour=MUTED,size=base-2),
        axis.title=element_text(colour=INK2,size=base-1),
        plot.title=element_text(colour=INK,face="bold",size=base+2),
        plot.subtitle=element_text(colour=INK2,size=base-1),
        plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
        strip.text=element_text(colour=INK2,face="bold",size=base-1),
        legend.text=element_text(colour=INK2,size=base-2),
        legend.title=element_text(colour=INK2,size=base-2),
        legend.position="top", legend.justification="left")
theme_hm <- function(base=10) theme_viz(base) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=base-3),
        axis.text.y=element_text(size=base-4), axis.ticks=element_blank())
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF,limitsize=FALSE)
                       cat("wrote ",file.path(FIG,n),"\n",sep="")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote ",file.path(OUT,n),"\n",sep="")}

cat("============================================================\n")
cat("step 14 -- cells inside light-chain neighborhoods, sample ", SAMPLE, "\n", sep="")
cat("============================================================\n")

## ---- 1. cells, clusters, coordinates ----------------------------------------
for (f in c(IN_RDS, CELLS_F, NBHD_F)) if (!file.exists(f)) stop("missing ", f, call.=FALSE)
cat("reading ", IN_RDS, "\n", sep="")
obj <- readRDS(IN_RDS); rn <- rownames(obj@meta.data)
NEED <- c("x_slide_mm","y_slide_mm","nCount_RNA","inSitu_clusters",
          "pearson_clusters_low","pearson_clusters_high")
miss <- setdiff(NEED, colnames(obj@meta.data))
if (length(miss)) stop("meta.data missing: ", paste(miss, collapse=", "), call.=FALSE)
md <- as.data.table(obj@meta.data[, NEED, drop = FALSE]); md[, cell_id := rn]
md[, `:=`(cl_low  = as.character(pearson_clusters_low),
          cl_high = as.character(pearson_clusters_high))]

cnt <- LayerData(obj[["RNA"]], layer = "counts")
stopifnot(identical(colnames(cnt), rn))
cat("materialising counts ", nrow(cnt), " x ", ncol(cnt), " ...\n", sep="")
M <- as(cnt, "dgCMatrix")
rm(obj, cnt); invisible(gc())

## Cluster `i` (the over-fixation defect) is 37% of IB1 and is NOT dropped before the
## footprint is drawn. A first run did drop it and the cost was severe: 14,641 called
## member cells fell to 6,127, and 154 neighborhoods to 68 -- including 11 of the 17
## kappa-restricted FL-like ones, which is the class the whole question is about. Step 07
## defined these neighborhoods WITH those cells, so removing them first redefines the
## object of study. Cluster `i` cells are instead flagged and excluded from every
## composition and expression statistic, while still anchoring the geometry.
n0 <- nrow(md)
keep <- md[, !is.na(x_slide_mm) & !is.na(y_slide_mm) & nCount_RNA > 0]
cat(sprintf("\ncells %d -> %d after dropping coordinate/count NAs (%d removed)\n",
            n0, sum(keep), n0 - sum(keep)))
md <- md[keep]; M <- M[, keep, drop = FALSE]
md[, is_defect := inSitu_clusters == DEFECT_CL]
cat(sprintf("cluster '%s' (over-fixation defect): %d cells (%.1f%%) -- kept for geometry, excluded from all statistics\n",
            DEFECT_CL, md[, sum(is_defect)], 100*md[, mean(is_defect)]))
stopifnot(identical(colnames(M), md$cell_id))

## ---- 2. IN / BUFFER / OUT ----------------------------------------------------
lc  <- fread(CELLS_F)
nbh <- fread(NBHD_F)
mem <- lc[nb %in% nbh$nb]                     # cells belonging to a CALLED neighborhood
cat(sprintf("\nlight-chain called cells %d | called neighborhoods %d | member cells %d\n",
            nrow(lc), nrow(nbh), nrow(mem)))
mem <- mem[cell_id %in% md$cell_id]
cat(sprintf("member cells present in this object after filtering: %d\n", nrow(mem)))
stopifnot(nrow(mem) > 100)

XY  <- as.matrix(md[, .(x = x_slide_mm*1000, y = y_slide_mm*1000)])   # um
MXY <- as.matrix(mem[, .(x = x_slide_mm*1000, y = y_slide_mm*1000)])
q   <- RANN::nn2(MXY, XY, k = 1)
md[, d_mem := as.numeric(q$nn.dists[,1])]
md[, nb_near := mem$nb[q$nn.idx[,1]]]
md[, grp := fifelse(d_mem <= R_IN, "in",
             fifelse(d_mem >  R_BUF, "out", "buffer"))]
md[grp != "in", nb_near := NA_integer_]
md[, is_member := cell_id %in% mem$cell_id]

cat("\n=== cells by region group\n")
print(md[, .(cells = .N, members = sum(is_member),
             med_nCount = as.numeric(median(nCount_RNA))), by = grp][order(-cells)])
cat(sprintf("\nIN region: %d cells across %d neighborhoods (%d are called members, %.1f%%)\n",
            md[grp=="in", .N], md[grp=="in", uniqueN(nb_near)],
            md[grp=="in", sum(is_member)], 100*md[grp=="in", mean(is_member)]))

## ---- 3. matched control regions ----------------------------------------------
## Each real neighborhood footprint holds n cells. Its controls take the n nearest OUT
## cells to a random OUT seed: same cell count, same compactness, different tissue.
## Matching on COUNT rather than radius is deliberate -- neighborhoods are dense, so a
## fixed-radius control would hold fewer cells and every fraction would be noisier.
nb_sizes <- md[grp == "in" & !(is_defect), .(n = .N), by = nb_near][!is.na(nb_near)]
out_idx  <- md[, which(grp == "out" & !is_defect)]
## When neighborhoods tile the whole section there is no control tissue and the in/out
## contrast is undefined. That is the expected result for a sample whose light-chain
## compartment IS the node (RB1, the FL recurrence), so say so and stop rather than
## fabricating controls out of a handful of leftover cells.
if (length(out_idx) < 20 * max(nb_sizes$n)) {
  cat(sprintf(paste0("\n*** STOPPING: only %d usable out-of-neighborhood cells against a ",
    "largest footprint of %d.\n*** Light-chain neighborhoods cover %.1f%% of this ",
    "section, so there is no control\n*** tissue and the in-vs-out contrast is ",
    "undefined for %s. This is a result, not a failure.\n"),
    length(out_idx), max(nb_sizes$n), 100*md[, mean(grp == "in")], SAMPLE))
  fwrite(data.table(sample = SAMPLE, cells = nrow(md), in_cells = md[, sum(grp=="in")],
                    out_cells = length(out_idx), pct_in = 100*md[, mean(grp=="in")],
                    verdict = "no out-territory; in/out contrast undefined"),
         file.path(OUT, "coverage_verdict.csv"))
  quit(save = "no", status = 0)
}
OXY      <- XY[out_idx, , drop = FALSE]
cat(sprintf("\nbuilding %d matched control regions from %d out-of-neighborhood cells\n",
            nrow(nb_sizes)*CTRL_PER, length(out_idx)))
## Cells inside neighborhoods are DEEPER than cells outside (see the group table above --
## dense, better-preserved follicular tissue). cp10k removes a uniform depth shift, but
## a shallow cell's cp10k is also noisier and inflates whatever it does detect, so an
## unmatched control set biases the comparison. Controls are therefore drawn from an
## oversampled candidate pool and chosen for the closest median depth to their target
## neighborhood, matching on depth AND cell count.
kmax  <- max(nb_sizes$n)
nseed <- min(nrow(OXY), nrow(nb_sizes) * CTRL_PER * CAND_MULT)
seeds <- sample.int(nrow(OXY), nseed, replace = FALSE)
knn   <- RANN::nn2(OXY, OXY[seeds, , drop = FALSE], k = min(kmax, nrow(OXY)))
odep  <- md$nCount_RNA[out_idx]
nb_dep <- md[grp == "in" & !(is_defect), .(dep = median(nCount_RNA)), by = nb_near]
nb_sizes <- merge(nb_sizes, nb_dep, by = "nb_near")
setnames(nb_sizes, c("n","dep"), c("n_nb","dep_nb"))
nb_sizes[, `:=`(n = n_nb, dep = dep_nb)]
used <- rep(FALSE, nseed)
ctrl <- rbindlist(lapply(seq_len(nrow(nb_sizes)), function(i) {
  n <- nb_sizes$n[i]; tgt <- nb_sizes$dep[i]
  cd <- vapply(seq_len(nseed), function(s)
         if (used[s]) Inf else abs(median(odep[knn$nn.idx[s, seq_len(n)]]) - tgt), 0)
  pick <- head(order(cd), CTRL_PER)
  used[pick] <<- TRUE
  rbindlist(lapply(seq_along(pick), function(j)
    data.table(region = sprintf("ctrl_%d_%d", nb_sizes$nb_near[i], j),
               cell_row = out_idx[knn$nn.idx[pick[j], seq_len(n)]])))
}))
## A cell may fall in more than one control region; that is fine for a region-level
## statistic but it is reported so the overlap is visible rather than assumed away.
cat(sprintf("control regions %d | cell slots %d | distinct cells %d (%.2fx reuse)\n",
            uniqueN(ctrl$region), nrow(ctrl), uniqueN(ctrl$cell_row),
            nrow(ctrl)/uniqueN(ctrl$cell_row)))
## Report the achieved depth match, because a matching step that silently failed would
## look identical to one that worked.
ctrl_dep <- ctrl[, .(dep = median(md$nCount_RNA[cell_row]), n = .N), by = region]
ctrl_dep[, nb_near := as.integer(sub("^ctrl_([0-9]+)_.*$", "\\1", region))]
setnames(ctrl_dep, c("dep","n"), c("dep_ctrl","n_ctrl"))
mm <- merge(ctrl_dep, nb_sizes[, .(nb_near, n_nb, dep_nb)], by = "nb_near")
cat(sprintf("\n=== depth match achieved (median nCount per region)\n  neighborhoods %.0f | controls %.0f | ratio %.3f\n  per-pair |log2 ratio|: median %.3f, 90th pct %.3f\n",
            median(mm$dep_nb), median(mm$dep_ctrl), median(mm$dep_ctrl)/median(mm$dep_nb),
            median(abs(log2(mm$dep_ctrl/mm$dep_nb))),
            quantile(abs(log2(mm$dep_ctrl/mm$dep_nb)), .9)))
cat(sprintf("  cell-count match: median %d vs %d per region\n",
            as.integer(median(mm$n_nb)), as.integer(median(mm$n_ctrl))))
sc("control_region_matching.csv", mm)

in_rows <- md[, which(grp == "in")]
reg <- rbind(
  data.table(region = paste0("nb_", md$nb_near[in_rows]), cell_row = in_rows, set = "in"),
  data.table(region = ctrl$region, cell_row = ctrl$cell_row, set = "out"))
reg[, cl_low := md$cl_low[cell_row]][, cl_high := md$cl_high[cell_row]]
reg[, is_member := md$is_member[cell_row]][, is_defect := md$is_defect[cell_row]]

## ============================================================================
## PART 1 -- cluster prevalence IN vs OUT
## ============================================================================
labs <- tryCatch(cluster_labels(SAMPLE), error = function(e) {
  cat("\nNOTE: ", conditionMessage(e), "\n", sep=""); NULL })
if (!is.null(labs)) {
  cat("\n=== step-09b marker-verified labels in use\n")
  print(labs[, .(clustering, cluster_id, label, status, top_z,
                 hierannot = hierannot_immune_core)])
}

prev_one <- function(clcol, drop_members) {
  d <- md[grp %in% c("in","out") & !(is_defect)]
  if (drop_members) d <- d[!(is_member)]
  tab <- d[, .N, by = c("grp", clcol)]
  tot <- tab[, .(tot = sum(N)), by = grp]
  tab <- merge(tab, tot, by = "grp")[, pct := 100*N/tot]
  setnames(tab, clcol, "cluster_id")
  w <- dcast(tab, cluster_id ~ grp, value.var = c("N","pct"))
  for (cn in c("N_in","N_out","pct_in","pct_out"))
    if (!cn %in% names(w)) w[, (cn) := 0]
  w[, log2_enr := log2((pct_in + .01) / (pct_out + .01))]
  w[, clustering := clcol][, members_dropped := drop_members][]
}
PREV <- rbindlist(lapply(c("cl_low","cl_high"), function(cc)
  rbindlist(lapply(c(FALSE, TRUE), function(dm) prev_one(cc, dm)))), fill = TRUE)
PREV[, clustering := fifelse(clustering=="cl_low","pearson_clusters_low","pearson_clusters_high")]
if (!is.null(labs)) {
  PREV[, cluster_id := as.integer(cluster_id)]
  PREV <- merge(PREV, labs[, .(cluster_id, clustering, label, status,
                               hierannot_immune_core)],
                by = c("cluster_id","clustering"), all.x = TRUE, sort = FALSE)
}
setorder(PREV, clustering, members_dropped, -log2_enr)
cat("\n============================================================\n")
cat("PART 1 -- cluster prevalence, IN vs OUT\n")
cat("============================================================\n")
for (cc in unique(PREV$clustering)) for (dm in c(FALSE, TRUE)) {
  cat(sprintf("\n=== %s | called member cells %s\n", cc,
              if (dm) "REMOVED (this is the interpretable one)" else "included (circular for B/plasma)"))
  print(PREV[clustering==cc & members_dropped==dm,
             .(cluster_id, pct_in = round(pct_in,2), pct_out = round(pct_out,2),
               log2_enr = round(log2_enr,2), n_in = N_in,
               label = if ("label" %in% names(PREV)) label else NA)])
}
sc("cluster_prevalence_in_out.csv", PREV)

## region-level significance: cluster fraction per region, neighborhoods vs controls
frac_reg <- function(clcol) {
  d <- reg[!(is_member) & !(is_defect)]
  f <- d[, .N, by = c("region","set",clcol)]
  t <- f[, .(tot = sum(N)), by = .(region, set)]
  f <- merge(f, t, by = c("region","set"))[, frac := N/tot]
  full <- CJ(region = unique(f$region), cl = unique(f[[clcol]]), unique = TRUE)
  setnames(full, "cl", clcol)
  f <- merge(full, f, by = c("region", clcol), all.x = TRUE)
  f[is.na(frac), frac := 0]
  f[, set := fifelse(grepl("^nb_", region), "in", "out")]
  f[]
}
test_prev <- function(clcol) {
  f <- frac_reg(clcol)
  f[, {
    a <- frac[set=="in"]; b <- frac[set=="out"]
    w <- suppressWarnings(wilcox.test(a, b))
    .(median_in = round(100*median(a),2), median_out = round(100*median(b),2),
      log2_enr = round(log2((median(a)+1e-4)/(median(b)+1e-4)),2),
      n_regions_in = length(a), p = w$p.value)
  }, by = c(clcol)][, q := p.adjust(p, "BH")][order(p)]
}
cat("\n=== region-level test (unit = region, NOT cell), members removed\n")
RT <- rbindlist(lapply(c("cl_low","cl_high"), function(cc) {
  r <- test_prev(cc); setnames(r, cc, "cluster_id"); r[, clustering := cc][] }), fill=TRUE)
print(RT)
sc("cluster_prevalence_region_test.csv", RT)

## ============================================================================
## PART 2 -- are neighborhoods consistent, and do any stick out?
## ============================================================================
cat("\n============================================================\n")
cat("PART 2 -- composition consistency across neighborhoods\n")
cat("============================================================\n")
fl <- frac_reg("cl_low")
CM <- dcast(fl, region + set ~ cl_low, value.var = "frac", fill = 0)
cmat <- as.matrix(CM[, setdiff(names(CM), c("region","set")), with = FALSE])
rownames(cmat) <- CM$region

## Bray-Curtis to the group centroid. Lower = more typical.
bc <- function(x, y) sum(abs(x-y)) / sum(x+y)
cent_in  <- colMeans(cmat[CM$set=="in", , drop=FALSE])
cent_out <- colMeans(cmat[CM$set=="out", , drop=FALSE])
CM[, d_centroid := vapply(seq_len(.N), function(i)
      bc(cmat[i,], if (set[i]=="in") cent_in else cent_out), 0)]
cat("\n=== dispersion around the group centroid (Bray-Curtis)\n")
print(CM[, .(regions = .N, median_d = round(median(d_centroid),3),
             iqr = round(IQR(d_centroid),3), max_d = round(max(d_centroid),3)), by = set])
dw <- suppressWarnings(wilcox.test(CM[set=="in", d_centroid], CM[set=="out", d_centroid]))
cat(sprintf("neighborhoods %s consistent than matched control regions (Wilcoxon p = %.3g)\n",
            if (median(CM[set=="in",d_centroid]) < median(CM[set=="out",d_centroid]))
              "MORE" else "LESS", dw$p.value))

NBC <- merge(CM[set=="in"], nbh[, .(region = paste0("nb_", nb), class, artifact_suspect,
                                    n_cells_called = n_cells, kappa_frac, span_mm,
                                    pct_cluster_i, med_nCount)],
             by = "region", all.x = TRUE)
thr <- quantile(NBC$d_centroid, .95)
NBC[, outlier := d_centroid >= thr]
cat(sprintf("\n=== %d neighborhoods flagged as compositional outliers (>=95th pct, BC >= %.3f)\n",
            sum(NBC$outlier), thr))
print(NBC[outlier == TRUE][order(-d_centroid)][,
      .(region, d_centroid = round(d_centroid,3), class, kappa_frac,
        n_cells_called, artifact_suspect)])
cat("\n=== composition dispersion by step-07 class\n")
print(NBC[, .(n = .N, median_d = round(median(d_centroid),3),
              outliers = sum(outlier)), by = class][order(-n)])
sc("neighborhood_composition.csv", NBC)

## does composition differ between FL-like and reactive-like neighborhoods?
cls <- NBC[, .N, by = class][N >= 5, class]
if (length(cls) >= 2 && all(c("kappa-restricted (FL-like)","polytypic (reactive-like)") %in% cls)) {
  cat("\n=== kappa-restricted (FL-like) vs polytypic (reactive-like) composition\n")
  cmp <- rbindlist(lapply(colnames(cmat), function(k) {
    a <- NBC[class=="kappa-restricted (FL-like)"][[k]]
    b <- NBC[class=="polytypic (reactive-like)"][[k]]
    data.table(cluster_id = k, median_FL = round(100*median(a),2),
               median_reactive = round(100*median(b),2),
               p = suppressWarnings(wilcox.test(a,b))$p.value)
  }))[, q := p.adjust(p,"BH")][order(p)]
  print(cmp); sc("neighborhood_class_composition.csv", cmp)
}

p_disp <- ggplot(CM, aes(d_centroid, fill = set)) +
  geom_density(alpha = .5, colour = NA) +
  scale_fill_manual(values = GRP_COLS, name = NULL,
                    labels = c(`in`="light-chain neighborhoods", out="matched control regions")) +
  labs(title = "Are light-chain neighborhoods a consistent kind of tissue?",
       subtitle = sprintf(paste0("Bray-Curtis distance from each region to its own group's ",
         "mean composition (%s, members removed).\nA tighter distribution means the ",
         "neighborhoods resemble each other more than random tissue does."),
         "pearson_clusters_low"),
       x = "Bray-Curtis distance to group centroid", y = "density",
       caption = sprintf("Data: %s/neighborhood_composition.csv", OUT)) + theme_viz()
sf("fig_14_composition_dispersion.png", p_disp, 8.5, 5)

if (nrow(NBC) > 3) {
  mds <- cmdscale(dist(cmat[CM$set=="in", , drop=FALSE], method="manhattan"), k = 2)
  MD <- data.table(region = CM[set=="in", region], MDS1 = mds[,1], MDS2 = mds[,2])
  MD <- merge(MD, NBC[, .(region, class, outlier, n_cells_called)], by="region")
  p_mds <- ggplot(MD, aes(MDS1, MDS2, colour = class, size = n_cells_called)) +
    geom_point(alpha = .8, stroke = 0) +
    geom_point(data = MD[outlier==TRUE], shape = 21, colour = INK, fill = NA,
               stroke = .5, show.legend = FALSE) +
    scale_colour_manual(values = c(`kappa-restricted (FL-like)`=DIV_HI,
                                   `polytypic (reactive-like)`=DIV_LO,
                                   `intermediate`=ACC, `lambda-skewed`=MUTED), name=NULL) +
    scale_size_continuous(range = c(1,6), name = "called cells") +
    labs(title = "Neighborhood composition does not separate by light-chain class",
         subtitle = paste0("MDS of pearson_clusters_low composition, called members removed. ",
                           "Ringed points are the 95th-percentile\ncompositional outliers. ",
                           "If FL-like and reactive-like neighborhoods sat in different ",
                           "regions, the clone would have a distinct microenvironment."),
         caption = sprintf("Data: %s/neighborhood_composition.csv", OUT)) + theme_viz()
  sf("fig_14_composition_mds.png", p_mds, 9, 6)
}

## ============================================================================
## PART 3 -- expression IN vs OUT, per cluster
## ============================================================================
cat("\n============================================================\n")
cat("PART 3 -- expression IN vs OUT within each cluster\n")
cat("============================================================\n")
CLCOL <- "cl_low"     # 9 clusters -> 18 heatmap columns, legible
CLS   <- sort(unique(md[[CLCOL]]))

## pooled pseudobulk per (cluster, group), members removed
dd <- md[grp %in% c("in","out") & !(is_member) & !(is_defect)]
dd[, key := paste0("c", get(CLCOL), "__", grp)]
ridx <- match(dd$cell_id, colnames(M))
KF <- factor(dd$key)
IND <- sparseMatrix(i = ridx, j = as.integer(KF), x = 1,
                    dims = c(ncol(M), nlevels(KF)),
                    dimnames = list(NULL, levels(KF)))
PB  <- as.matrix(M %*% IND)                       # genes x (cluster,group)
tot <- colSums(PB)
CP  <- sweep(PB, 2, pmax(tot,1), "/") * 1e4       # counts per 10k, per group
## DETECTION RATE, which cp10k cannot substitute for. cp10k is COMPOSITIONAL: it sums to
## 10,000, so a real rise in one set of genes forces an apparent fall in everything else,
## and the fall lands hardest on the highest-abundance genes. That artifact is exactly
## what a first pass produced -- PHKG2, A1BG and DUSP2 (the top of the panel at 91-114
## cp10k) came out "down inside neighborhoods" in all 9 clusters at ~2x, while their
## DETECTION rate moved only 0.87-0.96x. Their cp10k fell because CXCL13, CLU, CD79A and
## MS4A1 genuinely rose. Detection rate is composition-free, so both are reported and a
## direction is only trusted when the two agree.
DET <- as.matrix((M > 0) %*% IND)
DET <- sweep(DET, 2, pmax(as.numeric(table(KF)[colnames(DET)]),1), "/")
ncell <- as.numeric(table(KF)[colnames(CP)])
cat("\npseudobulk groups (cluster x in/out), cells and total counts:\n")
print(data.table(group = colnames(CP), cells = ncell, counts = tot)[order(group)])
sc("pseudobulk_cp10k.csv", data.table(gene = rownames(CP), as.data.table(round(CP,4))))

## region-level pseudobulk, for significance with the region as the unit
rr <- reg[!(is_member) & !(is_defect)]
rr[, rkey := paste0(region, "__c", get(CLCOL))]
rn_ok <- rr[, .N, by = .(rkey, region, set, cl = get(CLCOL))][N >= MIN_CELLS_RC]
rr <- rr[rkey %in% rn_ok$rkey]
RK  <- factor(rr$rkey)
RIND <- sparseMatrix(i = rr$cell_row, j = as.integer(RK), x = 1,
                     dims = c(ncol(M), nlevels(RK)),
                     dimnames = list(NULL, levels(RK)))
RPB <- as.matrix(M %*% RIND)
RCP <- sweep(RPB, 2, pmax(colSums(RPB),1), "/") * 1e4
cat(sprintf("\nregion x cluster pseudobulks with >= %d cells: %d (of %d possible)\n",
            MIN_CELLS_RC, ncol(RCP), uniqueN(reg$region)*length(CLS)))
print(rn_ok[, .(pseudobulks = .N), by = .(cl, set)][order(cl, set)])

## vectorised Mann-Whitney per gene, within each cluster
mw <- function(A, B) {
  ## A, B: genes x regions. Rank across the pooled regions per gene row.
  X <- cbind(A, B); nA <- ncol(A); nB <- ncol(B)
  R <- t(apply(X, 1, rank))
  U  <- rowSums(R[, seq_len(nA), drop=FALSE]) - nA*(nA+1)/2
  mu <- nA*nB/2; sg <- sqrt(nA*nB*(nA+nB+1)/12)
  z  <- (U - mu)/sg
  list(z = z, p = 2*pnorm(-abs(z)))
}
## Pooled over all clusters first. Per-cluster tests below are the interesting ones but
## they lose rare clusters to power; this one is always available and always reported.
## The unit here is the REGION, not the region x cluster: pooling region x cluster keys
## would count the same neighborhood up to 9 times and reintroduce exactly the
## pseudo-replication this design exists to avoid.
gr   <- reg[!(is_member) & !(is_defect)]
GK   <- factor(gr$region)
GIND <- sparseMatrix(i = gr$cell_row, j = as.integer(GK), x = 1,
                     dims = c(ncol(M), nlevels(GK)), dimnames = list(NULL, levels(GK)))
GPB  <- as.matrix(M %*% GIND)
GCP  <- sweep(GPB, 2, pmax(colSums(GPB),1), "/") * 1e4
gset <- ifelse(grepl("^nb_", colnames(GCP)), "in", "out")
pi_ <- colnames(GCP)[gset=="in"]; po_ <- colnames(GCP)[gset=="out"]
DE_POOL <- NULL
if (length(pi_) >= MIN_REG && length(po_) >= MIN_REG) {
  st <- mw(GCP[, pi_, drop=FALSE], GCP[, po_, drop=FALSE])
  DE_POOL <- data.table(gene = rownames(GCP),
                        mean_in = rowMeans(GCP[, pi_, drop=FALSE]),
                        mean_out = rowMeans(GCP[, po_, drop=FALSE]),
                        z = st$z, p = st$p)[, q := p.adjust(p,"BH")]
  DE_POOL[, log2fc := log2((mean_in + .1)/(mean_out + .1))]
  DE_POOL[, circular := gene %in% CIRCULAR]
  setorder(DE_POOL, p)
  cat(sprintf("\n=== POOLED over all cell types, unit = region (%d neighborhoods vs %d controls)\n    %d genes q<0.05 & |log2FC|>0.5 of %d\n",
              length(pi_), length(po_), DE_POOL[q<.05 & abs(log2fc)>.5, .N], nrow(DE_POOL)))
  print(head(DE_POOL[!(circular)][, .(gene, mean_in=round(mean_in,2),
        mean_out=round(mean_out,2), log2fc=round(log2fc,2), q=signif(q,3))], 20))
  sc("expression_in_out_pooled.csv", DE_POOL)
}

testable <- rn_ok[, .N, by=.(cl,set)][, dcast(.SD, cl ~ set, value.var="N", fill=0)]
testable[, testable := (get("in") >= MIN_REG) & (out >= MIN_REG)]
cat("\n=== which clusters are testable per-cluster (>= ", MIN_REG,
    " regions each side)\n", sep="")
print(testable)
sc("de_cluster_testability.csv", testable)

DE <- rbindlist(lapply(CLS, function(CLI) {
  cl <- CLI
  ki <- rn_ok$rkey[rn_ok$cl == CLI & rn_ok$set == "in"]
  ko <- rn_ok$rkey[rn_ok$cl == CLI & rn_ok$set == "out"]
  if (length(ki) < MIN_REG || length(ko) < MIN_REG) {
    cat(sprintf("  cluster %s skipped: %d in / %d out regions\n", CLI, length(ki), length(ko)))
    return(NULL)
  }
  if (!all(c(paste0("c",CLI,"__in"), paste0("c",CLI,"__out")) %in% colnames(CP))) return(NULL)
  A <- RCP[, ki, drop=FALSE]; B <- RCP[, ko, drop=FALSE]
  st <- mw(A, B)
  data.table(cluster_id = cl, gene = rownames(RCP),
             cp10k_in  = CP[, paste0("c",cl,"__in")],
             cp10k_out = CP[, paste0("c",cl,"__out")],
             det_in    = DET[, paste0("c",cl,"__in")],
             det_out   = DET[, paste0("c",cl,"__out")],
             med_in = rowMeans(A), med_out = rowMeans(B),
             z = st$z, p = st$p, n_in = length(ki), n_out = length(ko))
}))
DE[, log2fc := log2((cp10k_in + .1)/(cp10k_out + .1))]
DE[, log2fc_det := log2((det_in + 1e-4)/(det_out + 1e-4))]
DE[, q := p.adjust(p, "BH"), by = cluster_id]
DE[, circular := gene %in% CIRCULAR]
## Concordant = cp10k and detection agree in direction and detection moves materially.
## Anything else is most likely the compositional squeeze described above.
## 0.2 was too permissive: A1BG cleared it in 2 of 9 clusters at -0.27/-0.22 while
## failing in the other 7, which was enough to pull it into the figure. Genes with a real
## change move detection by well over 1 log2 unit (CXCL13 +2.1, CLU +1.6), so 0.3 keeps
## every genuine signal and drops the marginal compositional passes.
DE[, concordant := sign(log2fc) == sign(log2fc_det) & abs(log2fc_det) > .3]
cat(sprintf("\nconcordance: %d of %d significant gene x cluster results have cp10k and detection agreeing\n",
            DE[q < .05 & abs(log2fc) > .5 & concordant, .N],
            DE[q < .05 & abs(log2fc) > .5, .N]))
cat("=== significant but DISCORDANT (cp10k moves, detection does not) -- top 10\n")
print(head(DE[q < .05 & abs(log2fc) > .5 & !concordant][order(-abs(log2fc)),
      .(cluster_id, gene, cp10k_in = round(cp10k_in,1), cp10k_out = round(cp10k_out,1),
        log2fc = round(log2fc,2), det_log2 = round(log2fc_det,2))], 10))
setorder(DE, cluster_id, p)
cat(sprintf("\ndifferential tests: %d gene x cluster combinations, %d clusters\n",
            nrow(DE), uniqueN(DE$cluster_id)))
cat("\n=== significant genes per cluster (q < 0.05 AND |log2FC| > 0.5)\n")
print(DE[q < .05 & abs(log2fc) > .5, .(n_sig = .N,
      up_in = sum(log2fc>0), down_in = sum(log2fc<0)), by = cluster_id][order(-n_sig)])
cat("\n=== top 25 by q, across all clusters\n")
print(head(DE[order(q)][, .(cluster_id, gene, cp10k_in = round(cp10k_in,2),
       cp10k_out = round(cp10k_out,2), log2fc = round(log2fc,2),
       q = signif(q,3), circular)], 25))
sc("expression_in_out_de.csv.gz", DE)

## ---- heatmaps ----------------------------------------------------------------
lab_low <- if (!is.null(labs))
  cluster_label_vec(SAMPLE, "pearson_clusters_low") else setNames(paste0("c", CLS), CLS)
col_lab <- function(cl, g) sprintf("%s %s", cl, g)
ORD <- unlist(lapply(CLS, function(cl) c(paste0("c",cl,"__in"), paste0("c",cl,"__out"))))
ORD <- ORD[ORD %in% colnames(CP)]

hm <- function(genes, title, sub, file, h) {
  genes <- intersect(genes, rownames(CP)); if (!length(genes)) return(invisible(NULL))
  X <- CP[genes, ORD, drop=FALSE]
  Z <- t(scale(t(log1p(X))))                       # z per gene across all columns
  Z[!is.finite(Z)] <- 0
  ord <- if (nrow(Z) > 2) rownames(Z)[hclust(dist(Z))$order] else rownames(Z)
  D <- as.data.table(as.table(Z)); setnames(D, c("gene","col","z"))
  D[, cl  := sub("^c(.*)__.*$", "\\1", col)]
  D[, grp := sub("^.*__", "", col)]
  D[, collab := factor(sprintf("%s | %s", cl, grp),
                       levels = unlist(lapply(CLS, function(k) sprintf("%s | %s", k, c("in","out")))))]
  D[, gene := factor(gene, levels = ord)]
  D[, panel := factor(sprintf("%s\n%s", cl, ifelse(is.na(lab_low[cl]),"",lab_low[cl])),
                      levels = unique(sprintf("%s\n%s", CLS, ifelse(is.na(lab_low[CLS]),"",lab_low[CLS]))))]
  lim <- max(abs(D$z), na.rm=TRUE)
  p <- ggplot(D, aes(grp, gene, fill = z)) + geom_tile() +
    facet_grid(~ panel, switch = "x") +
    scale_fill_gradient2(low = SCALE_LO, mid = "#f7f6f2", high = SCALE_HI, midpoint = 0,
                         limits = c(-lim, lim), name = "z of log1p(cp10k), per gene") +
    scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
    labs(title = title, subtitle = sub, x = NULL, y = NULL,
         caption = sprintf("Data: %s/pseudobulk_cp10k.csv | columns are pearson_clusters_low, in/out paired",
                           OUT)) +
    theme_hm() + theme(panel.spacing = unit(1.5,"pt"),
                       strip.placement = "outside",
                       strip.text.x = element_text(size = 7, lineheight = 1.05))
  sf(file, p, 11, h)
}

mk <- data.table(gene = unlist(MARKERS),
                 family = rep(names(MARKERS), lengths(MARKERS)))
mk <- mk[gene %in% rownames(CP)]
cat(sprintf("\nmarker panel: %d of %d genes present on the 6K panel\n",
            nrow(mk), length(unlist(MARKERS))))
sc("marker_panel_used.csv", mk)
hm(mk$gene,
   "Marker expression inside vs outside light-chain neighborhoods",
   paste0("Each cluster contributes two columns: cells IN a neighborhood footprint and ",
          "cells in matched OUT tissue.\nCalled member cells are REMOVED, so B/plasma ",
          "signal here is not the selection variable restated. IGKC/IGLC are the ",
          "selection\nvariable and are shown only for reference. Cluster labels are the ",
          "step-09b marker-verified calls;\nany label not confidently resolved carries its ",
          "status in parentheses."),
   "fig_14_heatmap_markers.png", max(6, .16*nrow(mk)))

## Rank by the DETECTION effect, not cp10k. cp10k ranking is itself compositional, so it
## systematically promotes the panel's highest-abundance genes; detection rate cannot be
## pushed around by what other genes are doing.
top <- DE[!(circular) & q < .05 & concordant == TRUE][order(-abs(log2fc_det))]
sel <- unique(c(head(top[log2fc_det>0][order(-log2fc_det)]$gene, TOP_DE/2),
                head(top[log2fc_det<0][order(log2fc_det)]$gene,  TOP_DE/2)))
if (length(sel) < 4) sel <- head(DE[!(circular) & concordant == TRUE][order(q)]$gene, TOP_DE)
sel <- unique(sel)
cat(sprintf("\ndifferential heatmap: %d genes selected (q<0.05, largest |log2FC|, circular excluded)\n",
            length(sel)))
hm(sel,
   "Genes differing inside vs outside light-chain neighborhoods",
   sprintf(paste0("Top %d genes by |log2 detection ratio| among those significant with the REGION as ",
     "the unit (%d neighborhoods vs %d matched\ncontrol regions, BH q < 0.05) -- not the ",
     "cell, which would be pseudo-replicated. Called members removed; IGKC/IGLC excluded\n",
     "as the selection variable. Restricted to genes where cp10k and DETECTION RATE agree, ",
     "which removes the\ncompositional artifact that made the panel's highest-abundance ",
     "genes look uniformly depressed."),
     length(sel), uniqueN(reg[set=="in"]$region), uniqueN(reg[set=="out"]$region)),
   "fig_14_heatmap_differential.png", max(6, .16*length(sel)))

cat("\n============================================================\n")
cat("OK -- outputs in ", OUT, "\n", sep="")
cat("============================================================\n")
