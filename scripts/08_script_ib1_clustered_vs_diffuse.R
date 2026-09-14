## ============================================================================
## 08_script_ib1_clustered_vs_diffuse.R -- IB1 only. What ARE the B cells in the
## clustered regions, and what do they sit in?
##
## WHY THIS EXISTS. scripts/08_script_lightchain_sample_compare.R settled the disease
## PATTERN in IB1 using light chain alone: clustered B cells sit at kappa_frac 0.757
## while diffuse B cells sit at 0.630, i.e. on the 0.62 polytypic baseline, whereas RB1
## (known extensive FL) is skewed in both compartments. So IB1's FL reads as FOCAL.
## That argument is made entirely from IGKC/IGLC1-2. It says the clustered cells are
## clonally skewed; it says nothing about what they are or what surrounds them.
##
## This script answers that with two readouts that do not use light chain at all:
##   B. curated B-cell state / FL marker panels
##   C. composition of the 30 um microenvironment
## and a third cell set (D) as a sensitivity check on the light-chain-high gate.
##
## THE MEASUREMENT PROBLEM, and why half this script is confound control.
## IB1 carries a large over-fixation defect population -- vendor inSitu_clusters "i",
## 37% of all cells, diagnosed in output/07_overfixation/ib1_overfixation_summary.md --
## and it is distributed very unevenly across the compartments. Measured on the r40um
## table before anything else was done:
##
##     group                       cells    % cluster i    mean nCount_RNA
##     clustered (any scored nb)  14,641          10.6                 879
##     diffuse (unscored)         21,693          34.9                 780
##
## A naive expression comparison between those two sets returns fixation and depth
## genes, not biology. So: cluster "i" is excluded from the primary and every table is
## also written with it retained; every score is a SHARE of the cell's transcriptome
## rather than a count; and every headline number is additionally standardised across
## nCount_RNA deciles. A difference that does not survive depth standardisation is
## reported as depth-driven, not corrected.
##
## GROUPS (3-way, not 2-way). "Clustered" in step 08a means "member of any scored
## neighbourhood", which lumps kappa-restricted foci in with reactive germinal centres.
## Comparing that against diffuse would largely measure follicular vs interfollicular
## architecture. So the kappa-restricted foci are split out:
##     FL-like           B_clean, scored nb, class == "kappa-restricted (FL-like)"
##     clustered, other  B_clean, scored nb, any other class
##     diffuse           B_clean, not in a scored nb
##     dropped           B_clean in an artifact_suspect nb (60-70% cluster "i")
##
## Usage: Rscript scripts/08_script_ib1_clustered_vs_diffuse.R [sample]
##   sbatch --mem=200G --cpus-per-task=2 --time=4:00:00 \
##     scripts/sbatch/run_rscript.sh scripts/08_script_ib1_clustered_vs_diffuse.R IB1
## or just:  scripts/08_submit_ib1_clustered_vs_diffuse.sh
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(RANN)
  library(ggplot2); library(data.table); setDTthreads(1)
})

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"

IN_RDS  <- sprintf("output/01_import_v1/seuratObject_%s_BPCells_full.Rds", SAMPLE)
LC_DIR  <- file.path("output/07_lightchain", SAMPLE)
## The step-02 output directory was renamed 02_filter_QC -> 02_spatial_QC partway
## through this project. Resolve whichever exists rather than hardcoding either, so this
## script works before, during and after that rename.
.qc_dir <- local({
  cand <- c("output/02_spatial_QC", "output/02_filter_QC")
  hit  <- cand[dir.exists(cand)]
  if (!length(hit)) stop("no step-02 QC directory found; looked for: ",
                         paste(cand, collapse=", "), call. = FALSE)
  hit[1]
})
QC_CSV  <- file.path(.qc_dir, "qc_spatial_cells.csv.gz")
OUT     <- "output/08_ib1_clustered_vs_diffuse"
FIG     <- file.path(OUT, "figures")

RADII      <- c(20, 30, 40, 60)   # um; all reported, so the radius is not a hidden choice
FOCUS_R    <- 40                  # primary
NN_UM      <- 30                  # microenvironment radius
K_NN       <- 80                  # candidate neighbours before the radius filter
MIN_NBR    <- 5                   # neighbours required before a composition is computed
DEFECT_CL  <- "i"                 # the over-fixation cluster, see header
GC_POS_N   <- 2                   # "GC-positive" = at least this many GC markers detected
SAMPLE_COL <- "Run_Tissue_name"

if (!file.exists(IN_RDS))
  stop("missing ", IN_RDS, " -- run scripts/01_import_cosmx_seurat.R first", call. = FALSE)
if (!dir.exists(LC_DIR))
  stop("missing ", LC_DIR, " -- run scripts/07_submit_lightchain_neighborhoods.sh ",
       SAMPLE, " first", call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

## ---- marker panels ---------------------------------------------------------
## Panel membership was checked against the IB1 BPCells row_names before writing this.
## Absent from the 6K panel and therefore not used: FDCSP, SERPINA9, GCSAM, MEF2B,
## MYBL1, TOX2. CD10 is on the panel under its gene symbol MME.
##
## BCL2 IS DELIBERATELY ITS OWN PANEL. BCL6 and MME say "germinal centre", which a
## reactive follicle satisfies just as well as FL does; they are scored in GC, where
## they belong. The FL discriminator is BCL2 expression IN a GC-phenotype cell --
## t(14;18) puts BCL2 under the IGH enhancer, and a normal GC B cell is BCL2-negative.
## Scoring BCL2 inside a combined "FL" panel would have let a simple GC signal carry it.
PANELS <- list(
  BCL2       = c("BCL2"),
  GC         = c("AICDA","BCL6","MME","RGS13","LMO2","CXCR4","CD83","FAS"),
  FL_mut_prx = c("EZH2","TNFRSF14","EPHA7","CREBBP","KMT2D"),
  NBC_mantle = c("IGHD","TCL1A","FCER2","CR2","SELL"),
  MBC        = c("CD27","TNFRSF13B"),
  PROLIF     = c("MKI67","TOP2A","PCNA","STMN1"),
  BCELL      = c("MS4A1","CD79A","CD79B","BANK1","CD19","CD22"),
  PC         = c("MZB1","XBP1","DERL3","JCHAIN","PRDM1","SDC1","CD38","IRF4",
                 "SEC11C","TNFRSF17","POU2AF1"))
## FL_mut_prx are FL MUTATION targets scored by EXPRESSION. Expression is not mutation
## status; this panel is exploratory and is labelled as such wherever it is plotted.
NEG_CTRL <- c("BCELL", "PC")   # must not differ; see the verdict block

## Tfh, for the microenvironment. Scored only on neighbours already called lineage=="T",
## so it is Tfh-among-T rather than a mush of T content and Tfh programme.
TFH <- c("CXCR5","PDCD1","ICOS","IL21","BTLA","SH2D1A","CD200","CD40LG")
## NO FDC SCORE. scripts/01_import_cosmx_seurat.R:82 dropped FDC from LINEAGE_SETS with
## "FDC dropped: CR2/FCER2 are B markers; nothing FDC-exclusive survives on this panel".
## That still holds -- CR2 would score mantle B cells as FDC -- so FDC content is simply
## not measurable here and is not reported rather than reported badly.

GENES_ALL <- unique(c(unlist(PANELS, use.names = FALSE), TFH))

## ---- house style (scripts/07_script_lightchain_neighborhoods.R:126-147) -----
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#184f95"; DIV_MID<-"#f0efec"; DIV_LO<-"#d03b3b"
GRP_LEVELS <- c("FL-like", "clustered, other", "diffuse")
## FL-like takes DIV_HI, the kappa colour used throughout step 7, so the three-way split
## reads consistently against the light-chain figures.
GRP_COLS <- c("FL-like" = DIV_HI, "clustered, other" = "#7ba3d0", "diffuse" = "#a9a7a1")
LIN_COLS <- c(B = "#4e79a7", T = "#59a14f", NK = "#b07aa1", myeloid = "#e15759",
              endo = "#76b7b2", fibro = "#edc948", PC = "#9c755f",
              ambiguous = "#c3c2b7", below_background = "#dedcd4",
              all_markers_0 = "#efeee8")

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
          strip.text=element_text(colour=INK2,face="bold",size=base-1),
          legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left", legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF)
                       cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

## ---- reusable helpers ------------------------------------------------------
## Detection rate per gene (scripts/sbatch/qc_gc_detection.R:28-36).
det <- function(m, idx) {
  data.table(gene = rownames(m),
             pct_detected = 100 * Matrix::rowSums(m[, idx, drop = FALSE] > 0) / length(idx),
             mean_count   = Matrix::rowMeans(m[, idx, drop = FALSE]))
}

## DEPTH STANDARDISATION -- the rigorous version of the depth control.
## Normalising by nCount_RNA corrects the confound; comparing within depth strata
## removes it (scripts/07_script_plasma_marker_check.R:187-200). Here that per-decile
## comparison is collapsed back to one number per group by re-weighting each group's
## decile means to the POOLED decile sizes, i.e. direct standardisation. Deciles in
## which a group has no cells drop out and the weights renormalise over the rest.
std_by_depth <- function(dt, valcol, grpcol = "group") {
  w <- dt[, .(w = .N), by = depth_dec]
  w[, w := w / sum(w)]
  s <- dt[, .(v = mean(get(valcol)), n = .N), by = c("depth_dec", grpcol)]
  s <- merge(s, w, by = "depth_dec")
  s[, .(std = sum(v * w) / sum(w), n = sum(n), deciles = .N), by = c(grpcol)]
}

## Wilcoxon that never aborts the run: tiny or degenerate groups return NA rather than
## throwing, so one empty panel cannot take the whole script down.
safe_w <- function(a, b) {
  if (length(a) < 3 || length(b) < 3) return(NA_real_)
  out <- tryCatch(stats::wilcox.test(a, b)$p.value, error = function(e) NA_real_)
  as.numeric(out)
}

cat("============================================================\n")
cat("step 08 -- ", SAMPLE, ": clustered vs diffuse B cells\n", sep = "")
cat("============================================================\n")

## ---- A1. group assignment from the step-07 tables --------------------------
## Cheap: these are csvs. Done for every radius before the object is opened, so a
## missing table fails in seconds rather than after a 1.9 GB read.
assign_groups <- function(R) {
  fc <- file.path(LC_DIR, sprintf("lightchain_cells_r%dum.csv", R))
  fn <- file.path(LC_DIR, sprintf("lightchain_neighborhoods_r%dum.csv", R))
  if (!all(file.exists(fc, fn)))
    stop("missing step-07 table(s):\n  ", paste(c(fc, fn)[!file.exists(c(fc, fn))],
         collapse = "\n  "), call. = FALSE)
  lc  <- fread(fc)
  nbt <- fread(fn)[, .(nb, class, artifact_suspect)]
  x   <- merge(lc, nbt, by = "nb", all.x = TRUE, sort = FALSE)
  ## nb is NA for plasma/halo cells and non-NA-but-unscored for B_clean cells that
  ## landed in a community below MIN_CELLS; both leave `class` NA after the merge.
  x[, group := fifelse(tier != "B_clean", NA_character_,
              fifelse(is.na(class),                          "diffuse",
              fifelse(artifact_suspect %in% TRUE,            "dropped",
              fifelse(class == "kappa-restricted (FL-like)", "FL-like",
                                                             "clustered, other"))))]
  x[, radius_um := R]
  x[]
}

GRP <- rbindlist(lapply(RADII, assign_groups), fill = TRUE)
asum <- GRP[!is.na(group), .(cells = .N,
                             pct_defect = round(100*mean(inSitu_clusters == DEFECT_CL), 1),
                             mean_nCount = round(mean(nCount_RNA)),
                             kappa_frac = round(mean(call[call != "tie"] == "kappa"), 3)),
            by = .(radius_um, group)][order(radius_um, group)]
cat("\n=== group assignment, all radii\n"); print(asum)
sc("group_assignment_summary.csv", asum)

pool <- GRP[radius_um == FOCUS_R & group %in% GRP_LEVELS]
pool[, group := factor(group, levels = GRP_LEVELS)]
cat(sprintf("\nprimary radius %d um: %d B_clean cells in the three groups (%d dropped as artifact_suspect)\n",
            FOCUS_R, nrow(pool),
            GRP[radius_um == FOCUS_R & group == "dropped", .N]))

## ---- A2. the object: metadata + marker counts ------------------------------
cat("\nreading ", IN_RDS, "\n", sep = "")
obj <- readRDS(IN_RDS)
rn  <- rownames(obj@meta.data)          # "c_2_100_1002" -- the id every join uses

NEED <- c(SAMPLE_COL, "x_slide_mm", "y_slide_mm", "inSitu_clusters", "nCount_RNA",
          "qcCellsPassed", "lineage", "strongB", "neighbor_leiden",
          "inSitu_neighborhood_assignment")
miss <- setdiff(NEED, colnames(obj@meta.data))
if (length(miss)) stop("meta.data is missing: ", paste(miss, collapse = ", "), call. = FALSE)
NICHE_COLS <- grep("^inSitu_neighbours_", colnames(obj@meta.data), value = TRUE)
md <- as.data.table(obj@meta.data[, c(NEED, NICHE_COLS), drop = FALSE])
md[, cell_id := rn]

cnt <- LayerData(obj[["RNA"]], layer = "counts")
stopifnot(identical(colnames(cnt), rn))
HAVE <- intersect(GENES_ALL, rownames(cnt))
if (length(HAVE) < length(GENES_ALL))
  message("genes absent from the panel and dropped: ",
          paste(setdiff(GENES_ALL, HAVE), collapse = ", "))
## one materialisation of every marker this script needs, over all cells -- the
## microenvironment stage scores NEIGHBOUR cells, which are not restricted to the pool
M <- as(cnt[HAVE, , drop = FALSE], "dgCMatrix")
rm(obj, cnt); invisible(gc())
colnames(M) <- rn

PANELS <- lapply(PANELS, intersect, HAVE)
PANELS <- PANELS[lengths(PANELS) > 0]
TFH_H  <- intersect(TFH, HAVE)
## BCL2 and GC carry the headline read (B2) and the verdict; without them this script
## has nothing to say, so refuse rather than emit a table with silent holes in it.
if (!all(c("BCL2","GC") %in% names(PANELS)) || !"BCL2" %in% HAVE)
  stop("BCL2 and/or the GC panel are absent from this panel -- the FL read cannot be made",
       call. = FALSE)
if (!length(TFH_H)) message("no Tfh markers on the panel -- section C Tfh will be empty")

md <- md[get(SAMPLE_COL) == SAMPLE]
md <- md[!is.na(x_slide_mm) & !is.na(y_slide_mm) & nCount_RNA > 0]
cat(sprintf("%s cells with usable coordinates and counts: %d\n", SAMPLE, nrow(md)))

## qc_passing_all -- the step-02/03 filter that took IB1 from 2.54M to 400k. Carried as
## a reported covariate, not a filter, since step 07 deliberately worked pre-filter.
if (file.exists(QC_CSV)) {
  qc <- fread(QC_CSV, select = c("cell", "qc_passing_all"))
  md <- merge(md, qc, by.x = "cell_id", by.y = "cell", all.x = TRUE, sort = FALSE)
} else {
  message("missing ", QC_CSV, " -- qc_passing_all will be reported as NA")
  md[, qc_passing_all := NA]
}

## per-cell panel shares, as % of the cell's transcriptome (never raw counts; see the
## depth-confound note at scripts/07_script_plasma_marker_check.R:46-57)
gidx <- match(md$cell_id, colnames(M))
stopifnot(!anyNA(gidx))
for (p in names(PANELS)) {
  v <- Matrix::colSums(M[PANELS[[p]], gidx, drop = FALSE])
  set(md, j = paste0("s_", p), value = 100 * as.vector(v) / md$nCount_RNA)
  set(md, j = paste0("n_", p),
      value = as.vector(Matrix::colSums(M[PANELS[[p]], gidx, drop = FALSE] > 0)))
}
md[, gc_pos   := n_GC  >= GC_POS_N]
md[, bcl2_pos := as.vector(M["BCL2", gidx]) > 0]
md[, tfh_pct  := 100 * as.vector(Matrix::colSums(M[TFH_H, gidx, drop = FALSE])) / nCount_RNA]

## ---- A3. join the groups onto the metadata ---------------------------------
P <- merge(pool[, .(cell_id, group, nb, call, pct_lc, d_plasma_um)],
           md, by = "cell_id", sort = FALSE)
cat(sprintf("pool cells matched into the object: %d of %d\n", nrow(P), nrow(pool)))
if (nrow(P) < 0.99 * nrow(pool))
  stop("more than 1% of pool cells failed to join -- cell_id convention mismatch",
       call. = FALSE)
P[, group := factor(group, levels = GRP_LEVELS)]
P[, defect := inSitu_clusters == DEFECT_CL]

cat("\n=== confounders across the three groups (primary radius)\n")
print(P[, .(cells = .N,
            pct_defect_i   = round(100*mean(defect), 1),
            pct_qc_passing = round(100*mean(qc_passing_all %in% TRUE), 1),
            med_nCount     = as.numeric(median(nCount_RNA)),
            kappa_frac     = round(mean(call[call != "tie"] == "kappa"), 3)), by = group][order(group)])

## THE PRIMARY ANALYSIS SET: cluster "i" excluded. Everything below is computed twice,
## once on Q (primary) and once on P (with the defect population retained), and both
## are written out -- the gap between them IS the size of the fixation confound.
Q <- P[defect == FALSE]
cat(sprintf("\nprimary set after excluding cluster \"%s\": %d cells (%.1f%% removed)\n",
            DEFECT_CL, nrow(Q), 100*(1 - nrow(Q)/nrow(P))))
print(Q[, .N, by = group][order(group)])
## Everything downstream reshapes on all three group names being present. Fail here
## with a readable message rather than inside a dcast.
empty <- setdiff(GRP_LEVELS, unique(as.character(Q$group)))
if (length(empty))
  stop("no cells left in group(s): ", paste(empty, collapse = ", "),
       " -- the 3-way split cannot be built", call. = FALSE)

## ============================================================================
## B. MARKER PANELS
## ============================================================================
cat("\n============================================================\n")
cat("B. marker panels\n")
cat("============================================================\n")

panel_table <- function(d, tag) {
  d <- copy(d)
  ## deciles are computed WITHIN the analysis set, so the standardisation target is
  ## that set's own depth distribution
  d[, depth_dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                       include.lowest = TRUE, labels = FALSE)]
  rbindlist(lapply(names(PANELS), function(p) {
    col <- paste0("s_", p)
    raw <- d[, .(cells = .N, med = median(get(col)), mean = mean(get(col))), by = group]
    std <- std_by_depth(d, col)
    ref <- d[group == "diffuse", get(col)]
    out <- merge(raw, std[, .(group, std_mean = std)], by = "group")
    out[, `:=`(panel = p, genes = paste(PANELS[[p]], collapse = "/"), set = tag,
               neg_control = p %in% NEG_CTRL)]
    out[, p_vs_diffuse := vapply(as.character(group), function(g)
          if (g == "diffuse") NA_real_ else safe_w(d[group == g, get(col)], ref), 0)]
    out[]
  }))
}
PT <- rbind(panel_table(Q, "primary (no cluster i)"), panel_table(P, "with cluster i"))
## BH within each set, over the tests actually performed -- the diffuse rows are the
## reference and carry no p-value, so they must not inflate n
PT[, q_vs_diffuse := p.adjust(p_vs_diffuse, "BH", n = sum(!is.na(p_vs_diffuse))), by = set]
PT[, std_ratio_vs_diffuse := round(std_mean / std_mean[group == "diffuse"], 3),
   by = .(set, panel)]
setcolorder(PT, c("set","panel","group","cells","med","mean","std_mean",
                  "std_ratio_vs_diffuse","p_vs_diffuse","q_vs_diffuse","neg_control","genes"))
cat("\n=== panel scores, % of transcriptome (std_mean is depth-standardised)\n")
print(PT[set == "primary (no cluster i)"][order(panel, group)][
        , .(panel, group, cells, med = round(med,3), std_mean = round(std_mean,3),
            ratio = std_ratio_vs_diffuse, q = signif(q_vs_diffuse,3), neg_control)])
sc("marker_panel_scores.csv", PT)

## ---- B2. the FL read: BCL2 positivity inside GC-phenotype cells ------------
## This is the histologic hallmark, not a composite score. A reactive germinal centre
## is BCL2-NEGATIVE; t(14;18) makes the FL clone BCL2-positive while it keeps the GC
## programme. So the question is BCL2 positivity CONDITIONAL on being GC-phenotype.
fl <- rbindlist(lapply(list(list(Q,"primary (no cluster i)"), list(P,"with cluster i")),
  function(z) {
    d <- copy(z[[1]])
    d[, depth_dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                         include.lowest = TRUE, labels = FALSE)]
    g <- d[, .(cells = .N,
               pct_gc_pos          = 100*mean(gc_pos),
               pct_bcl2_pos        = 100*mean(bcl2_pos),
               pct_bcl2_pos_in_gc  = 100*mean(bcl2_pos[gc_pos]),
               n_gc_pos            = sum(gc_pos)), by = group]
    sd1 <- std_by_depth(d, "bcl2_pos")[, .(group, std_pct_bcl2_pos = 100*std)]
    sd2 <- std_by_depth(d[gc_pos == TRUE], "bcl2_pos")[, .(group, std_pct_bcl2_in_gc = 100*std)]
    g <- merge(merge(g, sd1, by = "group"), sd2, by = "group", all.x = TRUE)
    g[, set := z[[2]]][]
  }))
cat("\n=== BCL2 positivity, overall and within GC-phenotype cells\n")
print(fl[set == "primary (no cluster i)"][order(group)])
sc("fl_bcl2_in_gc.csv", fl)

## ---- B3. per gene, so a panel result is traceable to its genes -------------
pg <- rbindlist(lapply(GRP_LEVELS, function(g) {
  idx <- match(Q[group == g, cell_id], colnames(M))
  d <- det(M, idx); d[, `:=`(group = g, cells = length(idx))]; d[]
}))
pg[, panel := vapply(gene, function(x) {
     hit <- names(PANELS)[vapply(PANELS, function(p) x %in% p, TRUE)]
     if (length(hit)) paste(hit, collapse = "+") else "TFH" }, "")]
pgw <- dcast(pg, gene + panel ~ group, value.var = "pct_detected")
setnames(pgw, GRP_LEVELS, c("FL_like","clustered_other","diffuse"))
pgw[, ratio_FLlike_vs_diffuse := round(FL_like / pmax(diffuse, 1e-9), 2)]
cat("\n=== per-gene detection rate (%), primary set, top and bottom by ratio\n")
print(rbind(head(pgw[order(-ratio_FLlike_vs_diffuse)], 8),
            head(pgw[order( ratio_FLlike_vs_diffuse)], 8)))
sc("marker_per_gene.csv", pg)

## ---- B4. the depth-decile table itself, kept as evidence -------------------
Qd <- copy(Q)
Qd[, depth_dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                      include.lowest = TRUE, labels = FALSE)]
dm <- rbindlist(lapply(names(PANELS), function(p)
  Qd[, .(panel = p, mean_score = mean(get(paste0("s_", p))), cells = .N),
     by = .(depth_dec, group)]))
sc("marker_depth_matched.csv", dm)
## how consistent is each panel across depth? the number quoted in the verdict
consist <- dcast(dm, panel + depth_dec ~ group, value.var = "mean_score")
setnames(consist, GRP_LEVELS, c("FL_like","clustered_other","diffuse"))
cons <- consist[!is.na(FL_like) & !is.na(diffuse),
                .(deciles = .N, fl_above_diffuse = sum(FL_like > diffuse)), by = panel]
cat("\n=== depth consistency: deciles in which FL-like exceeds diffuse\n"); print(cons)

## ---- B5. neighbourhood as the unit of replication --------------------------
## n = 1 patient. 976 FL-like cells are not 976 independent observations; the
## neighbourhood is the closest thing to a replicate this design has.
nbagg <- Q[group != "diffuse",
           c(.(group = first(group), cells = .N),
             lapply(.SD, median)), by = nb, .SDcols = paste0("s_", names(PANELS))]
nbagg <- nbagg[cells >= 20]
nbres <- rbindlist(lapply(names(PANELS), function(p) {
  a <- nbagg[group == "FL-like",          get(paste0("s_", p))]
  b <- nbagg[group == "clustered, other", get(paste0("s_", p))]
  data.table(panel = p, n_nb_FLlike = length(a), n_nb_other = length(b),
             med_FLlike = median(a), med_other = median(b), p = safe_w(a, b))
}))
nbres[, q := p.adjust(p, "BH", n = max(sum(!is.na(p)), 1L))]
cat(sprintf("\n=== neighbourhood-level test (n=%d FL-like nbs vs %d other, >=20 cells each)\n",
            nbagg[group == "FL-like", .N], nbagg[group == "clustered, other", .N]))
print(nbres)
sc("marker_neighbourhood_level.csv", nbres)

## ---- B figures -------------------------------------------------------------
pp <- PT[set == "primary (no cluster i)"]
pp[, lab := ifelse(neg_control, paste0(panel, "  [control]"), panel)]
p1 <- ggplot(pp, aes(group, std_mean, fill = group)) +
  geom_col(width = .68) +
  geom_text(aes(label = sprintf("%.3f", std_mean)), vjust = -0.35, size = 2.6, colour = INK2) +
  facet_wrap(~ lab, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = GRP_COLS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(title = "Marker panels across the three B-cell groups in IB1",
       subtitle = paste0("Depth-standardised share of each cell's transcriptome, cluster \"",
                         DEFECT_CL, "\" excluded. [control] panels must NOT differ:\n",
                         "BCELL is cell identity, PC is residual plasma contamination. ",
                         "FL_mut_prx is expression of FL mutation targets, not mutation status."),
       x = NULL, y = "% of transcriptome (depth-standardised mean)",
       caption = sprintf("Data: %s/marker_panel_scores.csv", OUT)) +
  theme_viz() + theme(axis.text.x = element_blank())
sf("fig_08_marker_panels_by_group.png", p1, 11, 6.4)

p2 <- ggplot(dm[panel %in% c("BCL2","GC","PROLIF","BCELL")],
             aes(depth_dec, mean_score, colour = group, group = group)) +
  geom_line(linewidth = .7) + geom_point(size = 1.6) +
  facet_wrap(~ panel, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = GRP_COLS, name = NULL) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Does the difference survive at every sequencing depth?",
       subtitle = paste0("Comparing within nCount_RNA deciles removes the depth confound ",
                         "rather than correcting for it.\nA panel whose lines cross or ",
                         "converge is depth-driven. BCELL is the identity control and should be flat."),
       x = "nCount_RNA decile", y = "% of transcriptome (mean)",
       caption = sprintf("Data: %s/marker_depth_matched.csv", OUT)) + theme_viz()
sf("fig_08_marker_depth_matched.png", p2, 13, 4.4)

p3 <- ggplot(pg[panel != "TFH"], aes(reorder(gene, pct_detected), pct_detected, fill = group)) +
  geom_col(position = position_dodge(.8), width = .74) + coord_flip() +
  scale_fill_manual(values = GRP_COLS, name = NULL) +
  labs(title = "Every panel gene individually, not just the aggregate score",
       subtitle = "So a panel-level result can be traced to one gene or to the whole set",
       x = NULL, y = "% of cells with >=1 count",
       caption = sprintf("Data: %s/marker_per_gene.csv", OUT)) + theme_viz(xgrid = TRUE)
sf("fig_08_per_gene_detection.png", p3, 8.5, 10)

## BCL2 map -- the FL hallmark, measured independently of light chain
mp <- Q[order(bcl2_pos)]
p4 <- ggplot(mp, aes(x_slide_mm, y_slide_mm)) +
  geom_point(data = mp[bcl2_pos == FALSE], colour = "#e4e3dc", size = .18, stroke = 0) +
  geom_point(data = mp[bcl2_pos == TRUE], aes(colour = group), size = .45, stroke = 0) +
  scale_colour_manual(values = GRP_COLS, name = NULL,
                      guide = guide_legend(override.aes = list(size = 3))) +
  coord_equal() +
  labs(title = "Where the BCL2-positive B cells are",
       subtitle = paste0("Pale grey = BCL2-negative B cells. BCL2 is the FL hallmark ",
                         "(t(14;18)) and is measured\nindependently of light chain, so this ",
                         "is not a restatement of the step-07 map."),
       x = "x (mm)", y = "y (mm)",
       caption = sprintf("Data: %s/marker_panel_scores.csv", OUT)) +
  theme_viz() + theme(panel.grid.major.y = element_blank())
sf("fig_08_bcl2_map.png", p4, 8.5, 9)

## ============================================================================
## C. MICROENVIRONMENT
## ============================================================================
cat("\n============================================================\n")
cat("C. microenvironment within ", NN_UM, " um\n", sep = "")
cat("============================================================\n")

## C1 -- computed. The kd-tree is built over ALL cells (that is what physically
## surrounds an index cell) but only the pool cells are queried, so the returned
## matrices are 36k x K rather than 2.5M x K.
microenv <- function(idx_dt, pool_dt, tag) {
  tree <- as.matrix(pool_dt[, .(x_slide_mm, y_slide_mm)])
  qry  <- as.matrix(idx_dt[, .(x_slide_mm, y_slide_mm)])
  nn   <- RANN::nn2(tree, qry, k = min(K_NN, nrow(tree)))
  keep <- nn$nn.dists <= NN_UM / 1000
  ## The index cells are themselves in the tree, so each query returns itself at
  ## distance 0. Left in, every index B cell would contribute one guaranteed B
  ## neighbour and inflate f_B by ~1/tot. Masked by identity rather than by distance,
  ## because duplicate coordinates make the distance-0 row ambiguous.
  selfmask <- nn$nn.idx == match(idx_dt$cell_id, pool_dt$cell_id)
  selfmask[is.na(selfmask)] <- FALSE
  keep <- keep & !selfmask
  ## saturation: if the K-th neighbour is still inside the radius the count is censored
  sat  <- mean(nn$nn.dists[, ncol(nn$nn.dists)] <= NN_UM / 1000)
  cat(sprintf("  %-22s index %6d | tree %7d | K=%d saturated for %.1f%% of index cells\n",
              tag, nrow(qry), nrow(tree), ncol(nn$nn.idx), 100*sat))
  lin  <- as.character(pool_dt$lineage)[nn$nn.idx]
  dim(lin) <- dim(nn$nn.idx)
  lin[!keep] <- NA_character_
  levs <- levels(droplevels(factor(pool_dt$lineage)))
  comp <- vapply(levs, function(L) rowSums(lin == L, na.rm = TRUE), numeric(nrow(qry)))
  tot  <- rowSums(keep)
  ## Tfh scored on T-lineage neighbours only
  tfh  <- matrix(pool_dt$tfh_pct[nn$nn.idx], nrow = nrow(qry))
  isT  <- lin == "T" & !is.na(lin)
  tfh[!isT] <- NA_real_
  list(comp = comp, tot = tot, n_T = rowSums(isT),
       tfh_mean = rowMeans(tfh, na.rm = TRUE), saturation = sat, levs = levs)
}

## primary: neighbours in the defect cluster are excluded, since a cell surrounded by
## unreadable cells has an unreadable microenvironment, not an empty one
nbr_pool <- md[inSitu_clusters != DEFECT_CL]
me <- microenv(Q, nbr_pool, "primary (no cluster i)")
me_all <- microenv(Q, md, "with cluster i")

ME <- data.table(cell_id = Q$cell_id, group = Q$group,
                 nb = Q$nb, nCount_RNA = Q$nCount_RNA,
                 n_neighbours = me$tot, n_T = me$n_T, tfh_mean = me$tfh_mean)
## A composition is a fraction OF the neighbours, so it is undefined when there are
## none. Dividing by pmax(tot, 1) instead would score an isolated cell as 0% of every
## lineage, which drags every fraction down in whichever group is sparsest -- and the
## diffuse group is sparsest BY DEFINITION, since not being in a dense patch is what
## put it there. That would have manufactured "clustered B cells sit among more T
## cells" out of pure geometry. So fractions are NA below MIN_NBR neighbours, and
## local density is reported as its own outcome instead of leaking into the others.
for (L in me$levs)     ME[, (paste0("f_",  L)) := me$comp[, L]     / me$tot]
ME[, n_neighbours_withi := me_all$tot]
for (L in me_all$levs) ME[, (paste0("fi_", L)) := me_all$comp[, L] / me_all$tot]

dens <- ME[, .(cells = .N,
               med_neighbours  = as.numeric(median(n_neighbours)),
               mean_neighbours = round(mean(n_neighbours), 2),
               pct_isolated    = round(100*mean(n_neighbours < MIN_NBR), 1)),
           by = group][order(group)]
cat("\n=== local density within ", NN_UM, " um -- an outcome in its own right, and the\n", sep = "")
cat("    reason composition fractions below ", MIN_NBR, " neighbours are dropped\n", sep = "")
print(dens)
sc("microenv_density.csv", dens)

MEok <- ME[n_neighbours >= MIN_NBR]
cat(sprintf("cells with >= %d neighbours, used for composition: %d of %d (%.1f%%)\n",
            MIN_NBR, nrow(MEok), nrow(ME), 100*nrow(MEok)/nrow(ME)))

fcols <- grep("^f_", names(ME), value = TRUE)
comp_sum <- MEok[, c(.(cells = .N, med_neighbours = as.numeric(median(n_neighbours))),
                     lapply(.SD, function(v) round(100*mean(v, na.rm = TRUE), 2))),
                 by = group, .SDcols = fcols][order(group)]
cat("\n=== lineage composition of the ", NN_UM, " um neighbourhood (% of neighbours)\n", sep = "")
print(comp_sum)

fi <- grep("^fi_", names(ME), value = TRUE)
comp_i <- ME[n_neighbours_withi >= MIN_NBR,
             c(.(cells = .N), lapply(.SD, function(v) round(100*mean(v, na.rm = TRUE), 2))),
             by = group, .SDcols = fi][order(group)]
cat("\n=== same, with cluster \"", DEFECT_CL, "\" neighbours retained (sensitivity)\n", sep = "")
print(comp_i)
sc("microenv_lineage_composition_withdefect.csv", comp_i)

## deciles over all index cells, so the composition and Tfh standardisations share one
## depth scale even though they apply different sufficiency filters
MEd <- copy(ME)
MEd[, depth_dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                       include.lowest = TRUE, labels = FALSE)]
comp_std <- rbindlist(lapply(fcols, function(f) {
  s <- std_by_depth(MEd[n_neighbours >= MIN_NBR], f); s[, feature := f][] }))
comp_std[, std_pct := round(100*std, 2)]
sc("microenv_lineage_composition.csv",
   merge(melt(comp_sum, id.vars = c("group","cells","med_neighbours"),
              variable.name = "feature", value.name = "pct_raw"),
         comp_std[, .(group, feature, pct_depth_std = std_pct)],
         by = c("group","feature"), all.x = TRUE))

## Tfh among T neighbours
tfh_sum <- ME[n_T >= 3, .(cells = .N, med_T_neighbours = as.numeric(median(n_T)),
                          mean_tfh_pct = round(mean(tfh_mean, na.rm = TRUE), 4)), by = group][order(group)]
tfh_std <- std_by_depth(MEd[n_T >= 3 & is.finite(tfh_mean)], "tfh_mean")
cat("\n=== Tfh programme in T-lineage neighbours (% of their transcriptome)\n")
print(merge(tfh_sum, tfh_std[, .(group, depth_std = round(std, 4))], by = "group"))
tfh_p <- safe_w(ME[group == "FL-like" & n_T >= 3, tfh_mean],
                ME[group == "diffuse" & n_T >= 3, tfh_mean])
cat(sprintf("FL-like vs diffuse, Wilcoxon p = %s\n", signif(tfh_p, 3)))
sc("microenv_tfh.csv", merge(tfh_sum, tfh_std[, .(group, depth_std_tfh_pct = std)], by = "group"))

## C2 -- vendor niches. Free, and independent of everything above EXCEPT that they were
## computed on inSitu_clusters, which contains the defect cluster. Corroboration only.
vn <- Q[, .N, by = .(group, niche = inSitu_neighborhood_assignment)]
vn[, pct_of_group := round(100*N/sum(N), 1), by = group]
vnw <- dcast(vn, niche ~ group, value.var = "pct_of_group", fill = 0)
cat("\n=== vendor inSitu_neighborhood_assignment, % of each group\n")
cat("    NOTE: computed on inSitu_clusters, which includes the defect cluster --\n")
cat("    corroboration only, not ground truth\n")
setorderv(vnw, "FL-like", -1)
print(vnw)
sc("microenv_vendor_niche.csv", vn)

## ---- C figures -------------------------------------------------------------
cl <- melt(comp_sum, id.vars = c("group","cells","med_neighbours"),
           variable.name = "lineage", value.name = "pct")
cl[, lineage := sub("^f_", "", lineage)]
cl <- merge(cl, dens[, .(group, med_neighbours_all = med_neighbours)], by = "group")
cl[, xlab := sprintf("%s\n(median %g neighbours)", group, med_neighbours_all)]
p5 <- ggplot(cl, aes(reorder(xlab, as.integer(group)), pct, fill = lineage)) +
  geom_col(width = .7) +
  scale_fill_manual(values = LIN_COLS, name = NULL) +
  labs(title = paste0("What surrounds each group of B cells, within ", NN_UM, " um"),
       subtitle = paste0("Marker-based `lineage` from step 01, over index cells with >= ",
                         MIN_NBR, " neighbours. Density differs between the\ngroups and is ",
                         "shown on the axis, because a fraction of nothing is not 0% -- see ",
                         "microenv_density.csv.\nDefect-cluster neighbours excluded. FDC is ",
                         "absent by design: nothing FDC-exclusive survives on this panel."),
       x = NULL, y = "% of neighbours",
       caption = sprintf("Data: %s/microenv_lineage_composition.csv", OUT)) + theme_viz()
sf("fig_08_niche_composition.png", p5, 8.5, 6)

tp <- ME[n_T >= 3 & is.finite(tfh_mean)]
p6 <- ggplot(tp, aes(group, tfh_mean + 0.01, fill = group)) +
  geom_violin(scale = "width", colour = NA, alpha = .85) +
  geom_boxplot(width = .12, outlier.shape = NA, fill = SURF, colour = INK2, linewidth = .3) +
  scale_y_log10() + scale_fill_manual(values = GRP_COLS, guide = "none") +
  labs(title = "Tfh programme in the T cells surrounding each group",
       subtitle = paste0("Scored on T-lineage neighbours only, so this is Tfh-among-T rather ",
                         "than T content.\nlog scale, +0.01 offset. Index cells with >=3 T ",
                         "neighbours within ", NN_UM, " um."),
       x = NULL, y = "Tfh markers, % of neighbour transcriptome",
       caption = sprintf("Data: %s/microenv_tfh.csv", OUT)) + theme_viz()
sf("fig_08_tfh_scores.png", p6, 6.4, 5)

p7 <- ggplot(vn, aes(niche, pct_of_group, fill = group)) +
  geom_col(position = position_dodge(.8), width = .74) +
  scale_fill_manual(values = GRP_COLS, name = NULL) +
  labs(title = "Vendor niche assignment across the three groups",
       subtitle = paste0("Independent of everything computed here, EXCEPT that the vendor ",
                         "niches were built on\ninSitu_clusters -- which contains the defect ",
                         "cluster. Corroboration, not ground truth."),
       x = NULL, y = "% of group",
       caption = sprintf("Data: %s/microenv_vendor_niche.csv", OUT)) +
  theme_viz() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
sf("fig_08_vendor_niche_enrichment.png", p7, 9, 5)

## ============================================================================
## D. SENSITIVITY -- all B cells in each region, not just light-chain-high ones
## ============================================================================
cat("\n============================================================\n")
cat("D. sensitivity: all B cells in each region\n")
cat("============================================================\n")
## The primary analysis runs on the step-07 pool, which requires lc_total >= 4 -- a
## selection for cells with detectable Ig. This repeats the headline panels on ALL B
## cells, assigning region membership spatially instead.

exclude <- GRP[radius_um == FOCUS_R & tier %in% c("plasma","halo"), unique(cell_id)]
allB <- md[(lineage == "B" | strongB %in% TRUE) & !cell_id %in% exclude]
cat(sprintf("all B cells (lineage==B or strongB), minus plasma/halo: %d\n", nrow(allB)))

mem <- Q[group != "diffuse", .(x_slide_mm, y_slide_mm, group)]
nnm <- RANN::nn2(as.matrix(mem[, .(x_slide_mm, y_slide_mm)]),
                 as.matrix(allB[, .(x_slide_mm, y_slide_mm)]), k = 1)
allB[, d_mem_um := nnm$nn.dists[, 1] * 1000]
allB[, group := factor(fifelse(d_mem_um <= FOCUS_R,
                               as.character(mem$group)[nnm$nn.idx[, 1]], "diffuse"),
                       levels = GRP_LEVELS)]
allB <- allB[inSitu_clusters != DEFECT_CL]
cat(sprintf("after excluding cluster \"%s\": %d\n", DEFECT_CL, nrow(allB)))
print(allB[, .(cells = .N, med_nCount = as.numeric(median(nCount_RNA)),
               pct_in_lc_pool = round(100*mean(cell_id %in% Q$cell_id), 1)), by = group][order(group)])

SENS <- rbind(panel_table(Q, "pool only (primary)")[, .(set, panel, group, cells, std_mean)],
              panel_table(allB, "all region B cells")[, .(set, panel, group, cells, std_mean)])
SENS[, ratio_vs_diffuse := round(std_mean / std_mean[group == "diffuse"], 3),
     by = .(set, panel)]
cat("\n=== pool-only vs all-region-B, depth-standardised panel means\n")
print(dcast(SENS, panel + group ~ set, value.var = "ratio_vs_diffuse")[order(panel, group)])
sc("sensitivity_allB_vs_pool.csv", SENS)

p8 <- ggplot(SENS[group != "diffuse"], aes(panel, ratio_vs_diffuse, fill = group)) +
  geom_hline(yintercept = 1, colour = INK, linetype = "22", linewidth = .4) +
  geom_col(position = position_dodge(.8), width = .74) +
  facet_wrap(~ set, nrow = 1) +
  scale_fill_manual(values = GRP_COLS, name = NULL) +
  labs(title = "Does the light-chain-high gate create the result?",
       subtitle = paste0("Panel score relative to diffuse. Left: the step-07 pool ",
                         "(lc_total >= 4). Right: every B cell in the\nsame regions, ",
                         "assigned spatially. Agreement means the gate is not driving it."),
       x = NULL, y = "ratio to diffuse (depth-standardised)",
       caption = sprintf("Data: %s/sensitivity_allB_vs_pool.csv", OUT)) +
  theme_viz() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
sf("fig_08_sensitivity_allB.png", p8, 11, 5.2)

## ============================================================================
## E. SUMMARY AND VERDICT
## ============================================================================
## every lookup below returns exactly one number or NA, so a panel or group that did not
## survive filtering degrades the summary row rather than breaking the vector length
.one <- function(v) if (length(v) == 1) as.numeric(v) else NA_real_
gv <- function(p, g, col = "std_ratio_vs_diffuse")
  .one(PT[set == "primary (no cluster i)" & panel == p & group == g, get(col)])
fv <- function(g, col) .one(fl[set == "primary (no cluster i)" & group == g, get(col)])
cv <- function(g, col) .one(comp_sum[group == g, if (col %in% names(comp_sum)) get(col) else NULL])

ctrl_max <- max(abs(1 - c(gv("BCELL","FL-like"), gv("BCELL","clustered, other"),
                          gv("PC","FL-like"),    gv("PC","clustered, other"))), na.rm = TRUE)
ok_ctrl  <- is.finite(ctrl_max) && ctrl_max <= 0.25
bcl2_fl  <- fv("FL-like",          "std_pct_bcl2_in_gc")
bcl2_df  <- fv("diffuse",          "std_pct_bcl2_in_gc")
bcl2_co  <- fv("clustered, other", "std_pct_bcl2_in_gc")
cons_bcl2 <- cons[panel == "BCL2"]

SUM <- data.table(
  metric = c("cells_FL_like","cells_clustered_other","cells_diffuse",
             "pct_defect_removed",
             "ratio_BCELL_FLlike_vs_diffuse","ratio_PC_FLlike_vs_diffuse",
             "ratio_GC_FLlike_vs_diffuse","ratio_BCL2_FLlike_vs_diffuse",
             "ratio_PROLIF_FLlike_vs_diffuse",
             "pct_bcl2_in_gc_FLlike","pct_bcl2_in_gc_clustered_other","pct_bcl2_in_gc_diffuse",
             "bcl2_deciles_FLlike_above_diffuse",
             "n_neighbourhoods_FLlike",
             "med_neighbours_FLlike","med_neighbours_diffuse",
             "microenv_pct_T_FLlike","microenv_pct_T_diffuse",
             "microenv_saturation_pct"),
  value = c(Q[group=="FL-like", .N], Q[group=="clustered, other", .N], Q[group=="diffuse", .N],
            round(100*(1 - nrow(Q)/nrow(P)), 1),
            gv("BCELL","FL-like"), gv("PC","FL-like"),
            gv("GC","FL-like"), gv("BCL2","FL-like"), gv("PROLIF","FL-like"),
            round(bcl2_fl,2), round(bcl2_co,2), round(bcl2_df,2),
            .one(cons_bcl2$fl_above_diffuse),
            nbagg[group=="FL-like", .N],
            .one(dens[group=="FL-like", med_neighbours]),
            .one(dens[group=="diffuse", med_neighbours]),
            cv("FL-like", "f_T"), cv("diffuse", "f_T"),
            round(100*me$saturation, 1)))
sc("summary.csv", SUM)

cat("\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")
print(SUM)
cat("\n")
cat(sprintf("[%s] negative controls: largest deviation from diffuse across BCELL and PC is %.2fx\n",
            ifelse(ok_ctrl, "ok", "NO"), ctrl_max))
if (!ok_ctrl)
  cat("     >25%% -- the groups differ in cell IDENTITY or plasma content, not state.\n",
      "     Every other panel above is uninterpretable. Fix the grouping before reading on.\n", sep = "")
cat(sprintf("[%s] BCL2 within GC-phenotype cells: FL-like %.1f%% vs diffuse %.1f%% (clustered-other %.1f%%)\n",
            ifelse(is.finite(bcl2_fl) && is.finite(bcl2_df) && bcl2_fl > bcl2_df, "ok", "--"),
            bcl2_fl, bcl2_df, bcl2_co))
cat(sprintf("[%s] BCL2 depth consistency: FL-like above diffuse in %d of %d deciles\n",
            ifelse(nrow(cons_bcl2) && cons_bcl2$fl_above_diffuse >= 7, "ok", "--"),
            if (nrow(cons_bcl2)) cons_bcl2$fl_above_diffuse else NA,
            if (nrow(cons_bcl2)) cons_bcl2$deciles else NA))
cat("\nRead the panel table, the depth-decile figure and the neighbourhood-level test\n")
cat("together. A panel that separates only in the pooled comparison, or only in a few\n")
cat("depth deciles, or not at all once neighbourhoods rather than cells are the unit,\n")
cat("has not shown anything. If the kappa-restricted foci match the rest of the\n")
cat("clustered pool on every panel, that is the finding: they are not transcriptionally\n")
cat("separable at 6K-panel depth in this sample. It would not contradict step 08a,\n")
cat("which located them by light chain and never claimed they were.\n")
cat("============================================================\n")
cat("OK\n")
