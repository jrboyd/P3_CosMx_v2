## ============================================================================
## 08_script_ib1_hierannot_bcells.R -- the clustered-vs-diffuse comparison, redone on
## a POSITIVE B-cell definition.
##
## WHY THIS SUPERSEDES scripts/08_script_ib1_clustered_vs_diffuse.R.
## That script asked the right question and failed its own negative-control gate: its
## "diffuse" compartment scored 2.9x lower on B-identity markers than its clustered
## one, because step 07's `B_clean` tier is defined ENTIRELY BY EXCLUSION -- light-
## chain-high, not plasma, not within 15 um of plasma -- and never tests for B identity.
## The marker `lineage` call put that compartment at 13.4% B cells versus 48.9% for
## clustered. Nothing contrasted against it could be read.
##
## Step 09 supplies the missing piece: a cluster-level B call from HierAnnot
## (immune_core on pearson_clusters_low), 31,651 cells whose B markers run at 0.642%
## of transcriptome against 0.152% for the next lineage -- a 4.2x margin -- and which
## 6 of 6 hierarchy x resolution combinations place in B lineage for 95% of the cells.
## This script rebuilds the comparison on those cells.
##
## WHAT CHANGES, AND WHAT THAT COSTS.
##   * The cell set is now HierAnnot B cells, not the step-07 light-chain pool. Group
##     membership becomes SPATIAL: nearest scored-neighbourhood member within the
##     graph radius. The light-chain gate (lc_total >= 4) no longer selects the cells.
##   * BCELL IS NO LONGER AN INDEPENDENT CONTROL. Selecting cells on B identity
##     guarantees the B panel is flat; that flatness is now a sanity check on the join,
##     not evidence. PC and the depth standardisation remain live controls, and the
##     old pool-based set is carried alongside so the change is visible.
##   * COVERAGE IS THE BINDING LIMIT. The HierAnnot call comes from the POST-QC object
##     (400,111 cells). `qc_passing_all` coverage of the step-08 groups is FL-like
##     45.4% / clustered-other 44.3% / diffuse 18.0%, so the diffuse arm rests on a
##     filtered fifth of that compartment and is reported but NOT treated as decisive.
##     Fixing that needs the pre-QC object re-clustered, which is not done here.
##
## Usage: Rscript scripts/08_script_ib1_hierannot_bcells.R [sample]
##   scripts/08_submit_ib1_hierannot_bcells.sh
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(RANN)
  library(ggplot2); library(data.table); setDTthreads(1)
})

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"

IN_RDS  <- sprintf("output/01_import_v1/seuratObject_%s_BPCells_full.Rds", SAMPLE)
LC_DIR  <- file.path("output/07_lightchain", SAMPLE)
BCALL   <- file.path("output/09_hierannot", SAMPLE, "bcell_call.csv.gz")
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
OUT     <- "output/08_ib1_hierannot_bcells"
FIG     <- file.path(OUT, "figures")

RADII      <- c(20, 30, 40, 60)
FOCUS_R    <- 40                  # mm/1000; also the region-assignment radius
NN_UM      <- 30                  # microenvironment radius
K_NN       <- 80
MIN_NBR    <- 5                   # neighbours required before a composition is computed
DEFECT_CL  <- "i"
GC_POS_N   <- 2
SAMPLE_COL <- "Run_Tissue_name"

for (f in c(IN_RDS, BCALL)) if (!file.exists(f))
  stop("missing ", f, if (identical(f, BCALL))
       " -- run scripts/09_submit_hierannot_bcell_call.sh first" else "", call. = FALSE)
if (!dir.exists(LC_DIR)) stop("missing ", LC_DIR, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

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
## BCELL is now a JOIN CHECK, not a control -- see the header. PC is still a real
## control: nothing in the cell selection used plasma markers.
JOIN_CHECK <- "BCELL"; NEG_CTRL <- "PC"
TFH <- c("CXCR5","PDCD1","ICOS","IL21","BTLA","SH2D1A","CD200","CD40LG")
GENES_ALL <- unique(c(unlist(PANELS, use.names = FALSE), TFH))

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#184f95"; DIV_LO<-"#d03b3b"
GRP_LEVELS <- c("FL-like", "clustered, other", "diffuse")
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

std_by_depth <- function(dt, valcol, grpcol = "group") {
  w <- dt[, .(w = .N), by = depth_dec]; w[, w := w / sum(w)]
  s <- dt[, .(v = mean(get(valcol)), n = .N), by = c("depth_dec", grpcol)]
  s <- merge(s, w, by = "depth_dec")
  s[, .(std = sum(v * w) / sum(w), n = sum(n), deciles = .N), by = c(grpcol)]
}
safe_w <- function(a, b) {
  if (length(a) < 3 || length(b) < 3) return(NA_real_)
  as.numeric(tryCatch(stats::wilcox.test(a, b)$p.value, error = function(e) NA_real_))
}
add_dec <- function(d) { d <- copy(d)
  d[, depth_dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                       include.lowest = TRUE, labels = FALSE)]; d[] }

cat("============================================================\n")
cat("step 08 (HierAnnot B cells) -- ", SAMPLE, "\n", sep = "")
cat("============================================================\n")

## ---- A. step-07 neighbourhood members, the anchors for region assignment ---
lc  <- fread(file.path(LC_DIR, sprintf("lightchain_cells_r%dum.csv", FOCUS_R)))
nbt <- fread(file.path(LC_DIR, sprintf("lightchain_neighborhoods_r%dum.csv", FOCUS_R)))[
         , .(nb, class, artifact_suspect)]
x <- merge(lc, nbt, by = "nb", all.x = TRUE, sort = FALSE)
x[, group := fifelse(tier != "B_clean", NA_character_,
            fifelse(is.na(class),                          "diffuse",
            fifelse(artifact_suspect %in% TRUE,            "dropped",
            fifelse(class == "kappa-restricted (FL-like)", "FL-like",
                                                           "clustered, other"))))]
MEM <- x[group %in% c("FL-like","clustered, other"),
         .(x_slide_mm, y_slide_mm, group, nb)]
cat(sprintf("step-07 neighbourhood member cells at r%dum: %d (%d FL-like, %d other)\n",
            FOCUS_R, nrow(MEM), MEM[group=="FL-like", .N], MEM[group!="FL-like", .N]))

## ---- B. the object, marker shares, and the HierAnnot B call ----------------
cat("\nreading ", IN_RDS, "\n", sep = "")
obj <- readRDS(IN_RDS)
rn  <- rownames(obj@meta.data)
NEED <- c(SAMPLE_COL, "x_slide_mm", "y_slide_mm", "inSitu_clusters", "nCount_RNA",
          "qcCellsPassed", "lineage", "strongB")
miss <- setdiff(NEED, colnames(obj@meta.data))
if (length(miss)) stop("meta.data is missing: ", paste(miss, collapse=", "), call.=FALSE)
md <- as.data.table(obj@meta.data[, NEED, drop = FALSE]); md[, cell_id := rn]

cnt <- LayerData(obj[["RNA"]], layer = "counts")
stopifnot(identical(colnames(cnt), rn))
HAVE <- intersect(GENES_ALL, rownames(cnt))
M <- as(cnt[HAVE, , drop = FALSE], "dgCMatrix"); colnames(M) <- rn
rm(obj, cnt); invisible(gc())
PANELS <- lapply(PANELS, intersect, HAVE); PANELS <- PANELS[lengths(PANELS) > 0]
TFH_H  <- intersect(TFH, HAVE)

md <- md[get(SAMPLE_COL) == SAMPLE][!is.na(x_slide_mm) & !is.na(y_slide_mm) & nCount_RNA > 0]
if (file.exists(QC_CSV)) {
  qc <- fread(QC_CSV, select = c("cell","qc_passing_all"))
  md <- merge(md, qc, by.x = "cell_id", by.y = "cell", all.x = TRUE, sort = FALSE)
} else md[, qc_passing_all := NA]

gidx <- match(md$cell_id, colnames(M)); stopifnot(!anyNA(gidx))
for (p in names(PANELS)) {
  sub <- M[PANELS[[p]], gidx, drop = FALSE]
  set(md, j = paste0("s_", p), value = 100 * as.vector(Matrix::colSums(sub)) / md$nCount_RNA)
  set(md, j = paste0("n_", p), value = as.vector(Matrix::colSums(sub > 0)))
}
md[, gc_pos   := n_GC >= GC_POS_N]
md[, bcl2_pos := as.vector(M["BCL2", gidx]) > 0]
md[, tfh_pct  := 100 * as.vector(Matrix::colSums(M[TFH_H, gidx, drop=FALSE])) / nCount_RNA]

BC <- fread(BCALL, select = c("cell_id","is_B","is_plasma","cluster_id","label"))
cat(sprintf("\nHierAnnot call: %d cells, %d B, %d plasma\n",
            nrow(BC), BC[, sum(is_B)], BC[, sum(is_plasma)]))
md <- merge(md, BC, by = "cell_id", all.x = TRUE, sort = FALSE)
md[, `:=`(is_B = is_B %in% TRUE, is_plasma = is_plasma %in% TRUE)]
cat(sprintf("of %d %s cells with coordinates, %d carry a HierAnnot call (%.1f%%)\n",
            nrow(md), SAMPLE, md[!is.na(cluster_id), .N], 100*md[, mean(!is.na(cluster_id))]))

## ---- C. the analysis set: HierAnnot B cells, assigned spatially ------------
B <- md[is_B == TRUE & inSitu_clusters != DEFECT_CL]
cat(sprintf("\nHierAnnot B cells outside cluster \"%s\": %d\n", DEFECT_CL, nrow(B)))
nn1 <- RANN::nn2(as.matrix(MEM[, .(x_slide_mm, y_slide_mm)]),
                 as.matrix(B[,   .(x_slide_mm, y_slide_mm)]), k = 1)
B[, d_mem_um := nn1$nn.dists[, 1] * 1000]
B[, group := factor(fifelse(d_mem_um <= FOCUS_R, as.character(MEM$group)[nn1$nn.idx[,1]],
                            "diffuse"), levels = GRP_LEVELS)]
B[, nb := fifelse(d_mem_um <= FOCUS_R, MEM$nb[nn1$nn.idx[,1]], NA_integer_)]
cat("\n=== the analysis set\n")
print(B[, .(cells = .N, med_nCount = as.numeric(median(nCount_RNA)),
            pct_qc_passing = round(100*mean(qc_passing_all %in% TRUE), 1),
            med_d_mem_um = round(median(d_mem_um))), by = group][order(group)])
empty <- setdiff(GRP_LEVELS, unique(as.character(B$group)))
if (length(empty)) stop("no HierAnnot B cells in group(s): ",
                        paste(empty, collapse = ", "), call. = FALSE)

## the old pool-based set, carried so the change from step 08a is visible
Qold <- merge(x[group %in% GRP_LEVELS, .(cell_id, group)], md, by = "cell_id", sort = FALSE)
Qold <- Qold[inSitu_clusters != DEFECT_CL]
Qold[, group := factor(group, levels = GRP_LEVELS)]

## ---- D. panels --------------------------------------------------------------
panel_table <- function(d, tag) {
  d <- add_dec(d)
  rbindlist(lapply(names(PANELS), function(p) {
    col <- paste0("s_", p)
    raw <- d[, .(cells = .N, med = median(get(col))), by = group]
    std <- std_by_depth(d, col)
    ref <- d[group == "diffuse", get(col)]
    out <- merge(raw, std[, .(group, std_mean = std)], by = "group")
    out[, `:=`(panel = p, set = tag,
               role = fifelse(p == JOIN_CHECK, "join check",
                      fifelse(p == NEG_CTRL, "negative control", "readout")))]
    out[, p_vs_diffuse := vapply(as.character(group), function(g)
          if (g == "diffuse") NA_real_ else safe_w(d[group == g, get(col)], ref), 0)]
    out[]
  }))
}
PT <- rbind(panel_table(B, "HierAnnot B (primary)"), panel_table(Qold, "step-07 pool (08a)"))
PT[, q_vs_diffuse := p.adjust(p_vs_diffuse, "BH", n = max(sum(!is.na(p_vs_diffuse)),1L)), by = set]
PT[, ratio := round(std_mean / std_mean[group == "diffuse"], 3), by = .(set, panel)]
setcolorder(PT, c("set","panel","role","group","cells","med","std_mean","ratio",
                  "p_vs_diffuse","q_vs_diffuse"))
cat("\n=== panel scores, % of transcriptome, depth-standardised\n")
print(PT[set == "HierAnnot B (primary)"][order(panel, group)][
        , .(panel, role, group, cells, std_mean = round(std_mean,3), ratio,
            q = signif(q_vs_diffuse,3))])
sc("marker_panel_scores.csv", PT)

cat("\n=== the same contrast, old cell set vs new\n")
print(dcast(PT[group != "diffuse"], panel + group ~ set, value.var = "ratio")[order(panel, group)])

## the FL read: BCL2 conditional on GC phenotype
fl <- rbindlist(lapply(list(list(B,"HierAnnot B (primary)"), list(Qold,"step-07 pool (08a)")),
  function(z) { d <- add_dec(z[[1]])
    g <- d[, .(cells = .N, pct_gc_pos = 100*mean(gc_pos),
               pct_bcl2_pos = 100*mean(bcl2_pos),
               pct_bcl2_in_gc = 100*mean(bcl2_pos[gc_pos]), n_gc_pos = sum(gc_pos)), by = group]
    s <- std_by_depth(d[gc_pos == TRUE], "bcl2_pos")[, .(group, std_bcl2_in_gc = 100*std)]
    merge(g, s, by = "group", all.x = TRUE)[, set := z[[2]]][] }))
cat("\n=== BCL2 positivity within GC-phenotype cells\n")
print(fl[set == "HierAnnot B (primary)"][order(group)])
sc("fl_bcl2_in_gc.csv", fl)

## depth-decile consistency
Bd <- add_dec(B)
dm <- rbindlist(lapply(names(PANELS), function(p)
  Bd[, .(panel = p, mean_score = mean(get(paste0("s_", p))), cells = .N),
     by = .(depth_dec, group)]))
sc("marker_depth_matched.csv", dm)
cw <- dcast(dm, panel + depth_dec ~ group, value.var = "mean_score")
setnames(cw, GRP_LEVELS, c("FL_like","clustered_other","diffuse"), skip_absent = TRUE)
cons <- cw[!is.na(FL_like) & !is.na(diffuse),
           .(deciles = .N, fl_above_diffuse = sum(FL_like > diffuse)), by = panel]
cat("\n=== depth consistency: deciles in which FL-like exceeds diffuse\n"); print(cons)

## neighbourhood as the unit of replication
nbagg <- B[!is.na(nb), c(.(group = first(group), cells = .N),
                         lapply(.SD, median)), by = nb, .SDcols = paste0("s_", names(PANELS))]
nbagg <- nbagg[cells >= 10]
nbres <- rbindlist(lapply(names(PANELS), function(p) {
  a <- nbagg[group == "FL-like", get(paste0("s_", p))]
  b <- nbagg[group == "clustered, other", get(paste0("s_", p))]
  data.table(panel = p, n_nb_FLlike = length(a), n_nb_other = length(b),
             med_FLlike = median(a), med_other = median(b), p = safe_w(a, b))
}))
nbres[, q := p.adjust(p, "BH", n = max(sum(!is.na(p)),1L))]
cat(sprintf("\n=== neighbourhood-level test (%d FL-like nbs vs %d other, >=10 B cells each)\n",
            nbagg[group=="FL-like", .N], nbagg[group=="clustered, other", .N]))
print(nbres)
sc("marker_neighbourhood_level.csv", nbres)

## ---- E. microenvironment ----------------------------------------------------
microenv <- function(idx_dt, pool_dt, tag) {
  tree <- as.matrix(pool_dt[, .(x_slide_mm, y_slide_mm)])
  qry  <- as.matrix(idx_dt[, .(x_slide_mm, y_slide_mm)])
  nn   <- RANN::nn2(tree, qry, k = min(K_NN, nrow(tree)))
  keep <- nn$nn.dists <= NN_UM / 1000
  ## index cells are in the tree; without this every B cell donates itself one
  ## guaranteed B neighbour
  sm <- nn$nn.idx == match(idx_dt$cell_id, pool_dt$cell_id); sm[is.na(sm)] <- FALSE
  keep <- keep & !sm
  sat <- mean(nn$nn.dists[, ncol(nn$nn.dists)] <= NN_UM/1000)
  cat(sprintf("  %-24s index %6d | tree %7d | saturated %.1f%%\n",
              tag, nrow(qry), nrow(tree), 100*sat))
  lin <- as.character(pool_dt$lineage)[nn$nn.idx]; dim(lin) <- dim(nn$nn.idx)
  lin[!keep] <- NA_character_
  levs <- levels(droplevels(factor(pool_dt$lineage)))
  comp <- vapply(levs, function(L) rowSums(lin == L, na.rm=TRUE), numeric(nrow(qry)))
  tfh <- matrix(pool_dt$tfh_pct[nn$nn.idx], nrow = nrow(qry))
  isT <- lin == "T" & !is.na(lin); tfh[!isT] <- NA_real_
  list(comp = comp, tot = rowSums(keep), n_T = rowSums(isT),
       tfh_mean = rowMeans(tfh, na.rm = TRUE), levs = levs)
}
cat("\n=== microenvironment within ", NN_UM, " um\n", sep = "")
me <- microenv(B, md[inSitu_clusters != DEFECT_CL], "primary (no cluster i)")
ME <- data.table(cell_id = B$cell_id, group = B$group, nb = B$nb,
                 nCount_RNA = B$nCount_RNA,
                 n_neighbours = me$tot, n_T = me$n_T, tfh_mean = me$tfh_mean)
## a fraction of no neighbours is undefined, not zero -- and the diffuse group is
## sparsest by definition, so pmax(tot,1) would fabricate a composition difference
for (L in me$levs) ME[, (paste0("f_", L)) := me$comp[, L] / me$tot]
dens <- ME[, .(cells = .N, med_neighbours = as.numeric(median(n_neighbours)),
               pct_isolated = round(100*mean(n_neighbours < MIN_NBR), 1)), by = group][order(group)]
cat("\n=== local density (an outcome in its own right)\n"); print(dens)
sc("microenv_density.csv", dens)

fcols <- grep("^f_", names(ME), value = TRUE)
MEok  <- ME[n_neighbours >= MIN_NBR]
comp_sum <- MEok[, c(.(cells = .N), lapply(.SD, function(v) round(100*mean(v, na.rm=TRUE), 2))),
                 by = group, .SDcols = fcols][order(group)]
cat("\n=== lineage composition of the neighbourhood (% of neighbours)\n"); print(comp_sum)
MEd <- add_dec(ME)
comp_std <- rbindlist(lapply(fcols, function(f) {
  s <- std_by_depth(MEd[n_neighbours >= MIN_NBR], f); s[, feature := f][] }))
sc("microenv_lineage_composition.csv",
   merge(melt(comp_sum, id.vars = c("group","cells"), variable.name = "feature",
              value.name = "pct_raw"),
         comp_std[, .(group, feature, pct_depth_std = round(100*std, 2))],
         by = c("group","feature"), all.x = TRUE))

tfh_sum <- ME[n_T >= 3, .(cells = .N, mean_tfh_pct = round(mean(tfh_mean, na.rm=TRUE), 4)),
              by = group][order(group)]
cat("\n=== Tfh programme in T-lineage neighbours (CELL level -- see the warning below)\n")
print(tfh_sum)
cat(sprintf("cell-level Wilcoxon, FL-like vs diffuse p = %s | vs clustered-other p = %s\n",
            signif(safe_w(ME[group=="FL-like" & n_T>=3, tfh_mean],
                          ME[group=="diffuse" & n_T>=3, tfh_mean]), 3),
            signif(safe_w(ME[group=="FL-like" & n_T>=3, tfh_mean],
                          ME[group=="clustered, other" & n_T>=3, tfh_mean]), 3)))
cat("!! Those cell-level p-values are NOT usable. The cells come from a handful of\n")
cat("   neighbourhoods in one patient, so they are pseudo-replicates: n is the number\n")
cat("   of neighbourhoods, not of cells. The table below is the honest test.\n")

## NEIGHBOURHOOD-LEVEL microenvironment test. Same unit of replication as the marker
## panels use, for the same reason.
menb <- ME[!is.na(nb) & n_neighbours >= MIN_NBR,
           c(.(group = first(group), cells = .N,
               tfh = median(tfh_mean[n_T >= 3], na.rm = TRUE)),
             lapply(.SD, function(v) median(v, na.rm = TRUE))),
           by = nb, .SDcols = fcols][cells >= 10]
menbres <- rbindlist(lapply(c("tfh", fcols), function(f) {
  a <- menb[group == "FL-like",          get(f)]; a <- a[is.finite(a)]
  b <- menb[group == "clustered, other", get(f)]; b <- b[is.finite(b)]
  data.table(feature = f, n_nb_FLlike = length(a), n_nb_other = length(b),
             med_FLlike = if (length(a)) median(a) else NA_real_,
             med_other  = if (length(b)) median(b) else NA_real_,
             p = safe_w(a, b))
}))
menbres[, q := p.adjust(p, "BH", n = max(sum(!is.na(p)), 1L))]
cat(sprintf("\n=== microenvironment, NEIGHBOURHOOD level (%d FL-like vs %d other)\n",
            menb[group=="FL-like", .N], menb[group=="clustered, other", .N]))
print(menbres)
sc("microenv_tfh.csv", tfh_sum)
sc("microenv_neighbourhood_level.csv", menbres)

## ---- F. figures --------------------------------------------------------------
pp <- PT[set == "HierAnnot B (primary)"]
pp[, lab := fifelse(role == "readout", panel, paste0(panel, "  [", role, "]"))]
p1 <- ggplot(pp, aes(group, std_mean, fill = group)) +
  geom_col(width = .68) +
  geom_text(aes(label = sprintf("%.3f", std_mean)), vjust = -0.35, size = 2.6, colour = INK2) +
  facet_wrap(~ lab, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = GRP_COLS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(title = "Marker panels across HierAnnot-defined B cells in IB1",
       subtitle = paste0("Depth-standardised share of transcriptome, cluster \"", DEFECT_CL,
                         "\" excluded. BCELL is a JOIN CHECK, not a control:\nselecting cells on ",
                         "B identity guarantees it is flat. PC is the live negative control."),
       x = NULL, y = "% of transcriptome (depth-standardised)",
       caption = sprintf("Data: %s/marker_panel_scores.csv", OUT)) +
  theme_viz() + theme(axis.text.x = element_blank())
sf("fig_08b_marker_panels.png", p1, 11, 6.4)

cmp <- dcast(PT[group != "diffuse"], panel + group ~ set, value.var = "ratio")
setnames(cmp, c("HierAnnot B (primary)","step-07 pool (08a)"), c("hierannot","pool"),
         skip_absent = TRUE)
p2 <- ggplot(melt(cmp, id.vars = c("panel","group"), variable.name = "set",
                  value.name = "ratio"),
             aes(panel, ratio, fill = set)) +
  geom_hline(yintercept = 1, colour = INK, linetype = "22", linewidth = .4) +
  geom_col(position = position_dodge(.8), width = .74) +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_manual(values = c(hierannot = DIV_HI, pool = MUTED), name = NULL) +
  labs(title = "What changed when the B cells were defined positively",
       subtitle = paste0("Panel score relative to diffuse. 'pool' is step 08a's light-chain ",
                         "tier, which failed its B-identity\ncontrol; 'hierannot' is the ",
                         "cluster-level B call from step 09."),
       x = NULL, y = "ratio to diffuse (depth-standardised)",
       caption = sprintf("Data: %s/marker_panel_scores.csv", OUT)) +
  theme_viz() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
sf("fig_08b_old_vs_new.png", p2, 12, 5.2)

p3 <- ggplot(dm[panel %in% c("BCL2","GC","PROLIF","BCELL")],
             aes(depth_dec, mean_score, colour = group, group = group)) +
  geom_line(linewidth = .7) + geom_point(size = 1.6) +
  facet_wrap(~ panel, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = GRP_COLS, name = NULL) + scale_x_continuous(breaks = 1:10) +
  labs(title = "Does the difference survive at every sequencing depth?",
       subtitle = "Within-decile comparison removes the depth confound rather than correcting it.",
       x = "nCount_RNA decile", y = "% of transcriptome (mean)",
       caption = sprintf("Data: %s/marker_depth_matched.csv", OUT)) + theme_viz()
sf("fig_08b_depth_matched.png", p3, 13, 4.4)

cl <- melt(comp_sum, id.vars = c("group","cells"), variable.name = "lineage", value.name = "pct")
cl[, lineage := sub("^f_", "", lineage)]
cl <- merge(cl, dens[, .(group, md_n = med_neighbours)], by = "group")
cl[, xlab := sprintf("%s\n(median %g neighbours)", group, md_n)]
p4 <- ggplot(cl, aes(reorder(xlab, as.integer(group)), pct, fill = lineage)) +
  geom_col(width = .7) + scale_fill_manual(values = LIN_COLS, name = NULL) +
  labs(title = paste0("What surrounds each group of B cells, within ", NN_UM, " um"),
       subtitle = paste0("Index cells with >= ", MIN_NBR, " neighbours. Density differs and is ",
                         "shown on the axis -- a fraction of\nnothing is not 0%. FDC is absent ",
                         "by design: nothing FDC-exclusive survives on this panel."),
       x = NULL, y = "% of neighbours",
       caption = sprintf("Data: %s/microenv_lineage_composition.csv", OUT)) + theme_viz()
sf("fig_08b_niche_composition.png", p4, 8.5, 6)

## ---- G. verdict ---------------------------------------------------------------
.one <- function(v) if (length(v) == 1) as.numeric(v) else NA_real_
gv <- function(p, g) .one(PT[set=="HierAnnot B (primary)" & panel==p & group==g, ratio])
fv <- function(g) .one(fl[set=="HierAnnot B (primary)" & group==g, std_bcl2_in_gc])
join_ok <- abs(1 - gv(JOIN_CHECK, "FL-like")) <= 0.25 &&
           abs(1 - gv(JOIN_CHECK, "clustered, other")) <= 0.25
pc_ok   <- abs(1 - gv(NEG_CTRL, "FL-like")) <= 0.25 &&
           abs(1 - gv(NEG_CTRL, "clustered, other")) <= 0.25

SUM <- data.table(
  metric = c("cells_FL_like","cells_clustered_other","cells_diffuse",
             "ratio_BCELL_FLlike","ratio_PC_FLlike","ratio_GC_FLlike","ratio_BCL2_FLlike",
             "ratio_PROLIF_FLlike","bcl2_in_gc_FLlike","bcl2_in_gc_clustered_other",
             "bcl2_in_gc_diffuse","bcl2_deciles_FLlike_above_diffuse",
             "n_nb_FLlike","med_neighbours_FLlike","med_neighbours_diffuse"),
  value = c(B[group=="FL-like",.N], B[group=="clustered, other",.N], B[group=="diffuse",.N],
            gv("BCELL","FL-like"), gv("PC","FL-like"), gv("GC","FL-like"),
            gv("BCL2","FL-like"), gv("PROLIF","FL-like"),
            round(fv("FL-like"),2), round(fv("clustered, other"),2), round(fv("diffuse"),2),
            .one(cons[panel=="BCL2", fl_above_diffuse]),
            nbagg[group=="FL-like",.N],
            .one(dens[group=="FL-like", med_neighbours]),
            .one(dens[group=="diffuse", med_neighbours])))
sc("summary.csv", SUM)
cat("\n============================================================\nSUMMARY\n")
print(SUM)
cat(sprintf("\n[%s] join check: BCELL within 25%% of diffuse (expected -- cells were selected on B identity)\n",
            ifelse(join_ok, "ok", "NO")))
cat(sprintf("[%s] negative control: PC within 25%% of diffuse\n", ifelse(pc_ok, "ok", "NO")))
if (!pc_ok) cat("     plasma content still differs between groups -- readouts are confounded by it.\n")
cat("\nThe FL-like vs clustered-other contrast is the one this design supports cleanly.\n")
cat("The diffuse arm reaches only ~18% of that compartment (post-QC coverage) and is\n")
cat("reported for completeness, not as a decisive comparison.\n")
cat("============================================================\nOK\n")
