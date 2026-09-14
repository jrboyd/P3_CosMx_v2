## ============================================================================
## 09_script_hierannot_bcell_call.R -- a positive, cluster-level B-cell call for IB1.
##
## WHY. Step 08 failed its own negative-control gate: the "diffuse" compartment
## scored 2.9x LOWER on B-identity markers than the clustered one, and the marker
## `lineage` call puts it at 13.4% B cells versus 48.9% for clustered. The cause is
## that step 07's `B_clean` tier is defined ENTIRELY BY EXCLUSION -- light-chain-high,
## not plasma, not within 15 um of plasma -- and never tests for B identity. In
## follicles that exclusion lands on B cells; interfollicularly it admits whatever
## picked up ambient Ig. Nothing downstream can be fixed without a B call that is
## made positively.
##
## Per-cell marker calls are not the answer either: at a median ~800 counts, `lineage`
## leaves 22% ambiguous and 23% below background, and the step-03 HCA tonsil transfer
## was DISABLED (scripts/03_script_filter_and_process.R:79-119) because it called 78%
## of RB1's strongB cells "epithelial". Cluster-level annotation pools hundreds to
## thousands of cells per call, which is the level this depth actually supports.
##
## HierAnnot (NanoString/Bruker CosMx-Analysis-Scratch-Space) scores cluster mean
## profiles against a builtin marker hierarchy. Run through reticulate against the
## `hierannot_env` virtualenv. API and both known gotchas follow
## /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx/10a_HierAnnot.R from the v1 project.
##
## SCOPE LIMIT, stated up front because it constrains what step 08 can then do.
## The input is the POST-QC object (400,111 cells). Step 07's light-chain pool comes
## from the PRE-QC object (2,540,021 cells) and its `qc_passing_all` coverage is
## badly unbalanced across the step-08 groups:
##      FL-like 45.4%   clustered-other 44.3%   diffuse 18.0%
## So this call reaches ~45% of both clustered groups but only ~18% of diffuse. The
## FL-like vs clustered-other contrast stays clean; any diffuse contrast built on it
## rests on a filtered slice and must say so. `pearson_clusters_*` exist only on this
## object -- they were computed in step 03 on filtered cells -- so they cannot be
## extended to the pre-QC pool without re-clustering.
##
## Usage: Rscript scripts/09_script_hierannot_bcell_call.R [sample]
##   scripts/09_submit_hierannot_bcell_call.sh
## ============================================================================

suppressPackageStartupMessages({
  library(reticulate)
})
## use_virtualenv MUST run before anything initialises Python.
VENV <- Sys.getenv("HIERANNOT_VENV", "hierannot_env")
use_virtualenv(VENV, required = TRUE)

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); setDTthreads(1)
})

args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"

IN_RDS       <- sprintf("output/03_filter_and_process/%s_03_HCA_cell_typed.Rds", SAMPLE)
CLUSTER_COLS <- c("pearson_clusters_low", "pearson_clusters_high")
OUT          <- file.path("output/09_hierannot", SAMPLE)
FIG          <- file.path(OUT, "figures")

## Marker evidence used to VALIDATE whatever HierAnnot returns. The label strings are
## never trusted on their own -- this project has already been burned once by a
## transfer that produced confident, wrong labels.
MARKERS <- list(
  B       = c("MS4A1","CD79A","CD79B","BANK1","CD19","CD22"),
  T       = c("CD3D","CD3E","CD3G","TRAC"),
  myeloid = c("CD68","CD163","LYZ","ITGAX","C1QA","C1QB","CD14"),
  PC      = c("MZB1","XBP1","PRDM1","SDC1","DERL3","TNFRSF17"),
  epi     = c("KRT5","KRT19","EPCAM","KRT8","KRT18","KRT14"),
  fibro   = c("COL1A1","COL1A2","PDGFRB","DCN","ACTA2"),
  endo    = c("PECAM1","VWF","CLDN5","CD34"))
## A cluster is accepted as B only if B markers are the strongest lineage AND clear
## this margin over the runner-up. Mirrors LINEAGE_MARGIN in step 01 (2:1 on log2).
B_MARGIN <- 1.5

if (!file.exists(IN_RDS))
  stop("missing ", IN_RDS, " -- run scripts/03_script_filter_and_process.R first",
       call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#184f95"; DIV_LO<-"#d03b3b"
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

## SEURAT (not HierAnnot) prefixes purely numeric group names with "g" -- see
## AverageExpression's "First group.by variable ... starts with a number, appending
## `g` to ensure valid variable names". Alphabetic names are left alone. A bare
## sub("g","",x) strips the first "g" ANYWHERE -- in the v1 project that silently
## turned "singleton" into "sinleton" and broke the join. Anchor to a leading g
## followed by a digit. pearson_clusters_low ("0".."8") and _high ("0".."17") are
## both purely numeric, so every cluster id here WILL come back prefixed.
unprefix_g <- function(x) sub("^g(?=[0-9])", "", as.character(x), perl = TRUE)

cat("============================================================\n")
cat("step 09 -- HierAnnot B-cell call, ", SAMPLE, "\n", sep = "")
cat("============================================================\n")
print(py_config())

hier <- import("hierannot")

## ---- 1. which builtin hierarchies to run -----------------------------------
HIER <- unique(c(hier$match_builtin_hierarchy("immune"),
                 hier$match_builtin_hierarchy("tonsil"),
                 hier$match_builtin_hierarchy("tme")))
HIER <- as.character(HIER)
cat("\nbuiltin hierarchies matched: ", paste(HIER, collapse = ", "), "\n", sep = "")
if (!length(HIER)) stop("no builtin hierarchies matched -- check the hierannot version",
                        call. = FALSE)

## ---- 2. the object and its cluster columns ---------------------------------
cat("\nreading ", IN_RDS, "\n", sep = "")
seu <- readRDS(IN_RDS)
cat(sprintf("%s: %d cells x %d features\n", SAMPLE, ncol(seu), nrow(seu)))

have <- intersect(CLUSTER_COLS, colnames(seu@meta.data))
if (!length(have)) {
  cand <- grep("cluster", colnames(seu@meta.data), value = TRUE, ignore.case = TRUE)
  stop("none of ", paste(CLUSTER_COLS, collapse = "/"), " are in this object.\n",
       "cluster-like columns present: ", paste(cand, collapse = ", "), call. = FALSE)
}
if (length(have) < length(CLUSTER_COLS))
  message("absent and skipped: ", paste(setdiff(CLUSTER_COLS, have), collapse = ", "))
for (cc in have)
  cat(sprintf("  %-24s %d clusters, %d cells with a call\n", cc,
              uniqueN(na.omit(seu@meta.data[[cc]])), sum(!is.na(seu@meta.data[[cc]]))))

## Seurat 5 renamed AverageExpression's `slot` argument to `layer`; the v1 reference
## script predates that. Pick whichever this install actually accepts.
avg_counts <- function(obj, group_col) {
  a <- if ("layer" %in% names(formals(Seurat::AverageExpression)))
         Seurat::AverageExpression(obj, assays = "RNA", layer = "counts",
                                   group.by = group_col, verbose = FALSE)
       else
         Seurat::AverageExpression(obj, assays = "RNA", slot = "counts",
                                   group.by = group_col, verbose = FALSE)
  as.data.frame(a$RNA)
}

## ---- 3. cluster mean profiles, once per clustering -------------------------
## The v1 script recomputed these inside the hierarchy loop; they do not depend on the
## hierarchy, and over 400k cells that is the expensive part.
CM <- list()
for (cc in have) {
  safe <- as.character(seu@meta.data[[cc]])
  safe[is.na(safe)] <- "NA"
  seu$hierannot_cluster_safe <- safe
  m <- avg_counts(seu, "hierannot_cluster_safe")
  m <- m[, colnames(m) != "NA", drop = FALSE]
  CM[[cc]] <- m
  cat(sprintf("cluster means for %-24s %d genes x %d clusters\n", cc, nrow(m), ncol(m)))
}
seu$hierannot_cluster_safe <- NULL

## ---- 4. per-cell marker evidence, aggregated to clusters -------------------
## Independent of HierAnnot. This is what decides whether a label is believable.
cnt <- LayerData(seu[["RNA"]], layer = "counts")
MK  <- lapply(MARKERS, intersect, rownames(cnt))
absent <- unlist(lapply(names(MARKERS), function(k) setdiff(MARKERS[[k]], MK[[k]])))
if (length(absent)) message("markers absent from the panel: ", paste(absent, collapse = ", "))
nC <- Matrix::colSums(cnt)
ev <- data.table(cell_id = colnames(seu), nCount = as.numeric(nC))
for (k in names(MK)) {
  if (!length(MK[[k]])) { ev[, (k) := NA_real_]; next }
  ## share of the cell's transcriptome, the depth-normalised form used throughout
  ## steps 07-08; raw counts would just re-measure sequencing depth
  ev[, (k) := 100 * as.numeric(Matrix::colSums(cnt[MK[[k]], , drop = FALSE])) / pmax(nC, 1)]
}
for (cc in have) ev[, (cc) := as.character(seu@meta.data[[cc]])]
rm(cnt); invisible(gc())

## ---- 5. run HierAnnot: every hierarchy x every clustering ------------------
res_flat  <- list()      # python result objects, for compare_hierarchy_fit
exports   <- list()      # R-side cluster tables
for (hn in HIER) {
  cat("\n---- hierarchy: ", hn, " ----\n", sep = "")
  root <- hier$get_builtin_hierarchy(hn)
  for (cc in have) {
    pipe <- hier$HierAnnotPipeline(root_programs = root, preset = "auto",
                                   input_type = "raw_cluster_means")
    result <- pipe$fit_score(CM[[cc]])
    exp_df <- hier$make_cluster_annotation_export_summary(
      result,
      label_with_cluster            = TRUE,
      unknown_on_low_confidence     = TRUE,
      unknown_branch_raw_threshold  = 0.1,
      mixed_on_parent_mixing        = TRUE,
      unresolved_candidate_min_score     = 1.0,
      unresolved_candidate_min_raw_score = 0.2)
    e <- as.data.table(as.data.frame(exp_df))
    e[, cluster_id := unprefix_g(cluster_id)]
    e[, `:=`(hierarchy = hn, clustering = cc)]
    key <- paste(hn, cc)
    res_flat[[key]] <- result
    exports[[key]]  <- e
    cat(sprintf("  %-24s %d clusters annotated\n", cc, nrow(e)))
  }
}
EX <- rbindlist(exports, fill = TRUE)
## guard the g-prefix strip actually reconnected the ids
for (cc in have) {
  ids <- unique(EX[clustering == cc, cluster_id])
  obs <- unique(na.omit(as.character(seu@meta.data[[cc]])))
  if (!length(intersect(ids, obs)))
    stop("cluster ids from HierAnnot do not match ", cc, " after the g-prefix strip.\n",
         "  hierannot: ", paste(head(ids, 5), collapse = ","), "\n",
         "  seurat   : ", paste(head(obs, 5), collapse = ","), call. = FALSE)
}
sc("hierannot_cluster_labels.csv", EX)

## ---- 6. fit comparison -- MUST run while the Python session is live --------
cmp <- tryCatch({
  d <- as.data.table(as.data.frame(hier$compare_hierarchy_fit(results = res_flat)))
  if ("name" %in% names(d)) {
    d[, hierarchy  := sub(" .*$", "", name)]
    d[, clustering := sub("^\\S+ ", "", name)]
  }
  d
}, error = function(e) { message("compare_hierarchy_fit failed: ", conditionMessage(e))
                         data.table() })
if (nrow(cmp)) { cat("\n=== hierarchy fit comparison\n"); print(cmp); sc("hierannot_fit_comparison.csv", cmp) }

## ---- 7. derive B lineage from the PATH, then validate against the markers --
## COLUMN CHOICE, learned the hard way. `annot_l1_label`/`annot_l2_label` are NOT this
## row's own assignment -- they report the best-scoring program at each level
## independent of the branch the pipeline actually chose, so they disagree with the
## call and with each other. Observed in this very run: tonsil_tme cluster 8 has
## annot_path "Immune > B cell > Plasma cell" but annot_l1_label "Fibroblast", and
## immune_core cluster 0 has annot_label "Mast cell" paired with annot_l2_label "Treg".
## The authority is `annot_path` (the chosen branch) with `annot_export_label` as the
## display form, because export carries the pipeline's own hedging -- "unknown",
## "mixed_<a>.<b>", "candidate_<x>" -- which is exactly the QC signal we want to keep.
EX[, is_B_path  := grepl("(^|>\\s*)B cell(\\s*>|$)", annot_path)]
EX[, is_PC_path := grepl("Plasma cell", annot_path)]
EX[, bcat := fifelse(is_B_path & is_PC_path, "B-lineage: plasma",
            fifelse(is_B_path,               "B-lineage: B cell", "not B"))]
## Plasma cells are B-lineage but are precisely what step 07/08 must keep OUT of a
## B-cell comparison, so they get their own category rather than being folded in.
cat("\n=== path-derived B lineage per hierarchy x clustering (clusters)\n")
print(dcast(EX[, .N, by = .(hierarchy, clustering, bcat)],
            hierarchy + clustering ~ bcat, value.var = "N", fill = 0))
LAB_COL <- "annot_export_label"

lin <- names(MK)[lengths(MK) > 0]
val <- rbindlist(lapply(HIER, function(hn) rbindlist(lapply(have, function(cc) {
  map <- EX[hierarchy == hn & clustering == cc, .(cluster_id, label = get(LAB_COL))]
  d   <- merge(ev[, c("cell_id", cc, lin), with = FALSE], map,
               by.x = cc, by.y = "cluster_id", all.x = TRUE)
  d[!is.na(label), c(.(hierarchy = hn, clustering = cc, cells = .N),
                     lapply(.SD, function(v) round(mean(v), 4))),
    by = label, .SDcols = lin]
}))))
## which lineage does the marker evidence actually support, and by how much
val[, top_marker := lin[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = lin]
val[, top_val    := do.call(pmax, .SD), .SDcols = lin]
val[, second_val := apply(as.matrix(.SD), 1, function(r) sort(r, decreasing = TRUE)[2]),
    .SDcols = lin]
val[, margin := round(top_val / pmax(second_val, 1e-9), 2)]
val[, marker_supports_B := top_marker == "B" & margin >= B_MARGIN]
setcolorder(val, c("hierarchy","clustering","label","cells","top_marker","margin",
                   "marker_supports_B"))
cat("\n=== every HierAnnot label vs the independent marker evidence\n")
print(val[order(hierarchy, clustering, -cells)])
sc("hierannot_label_marker_validation.csv", val)

## ---- 8. validate the path-derived categories, then pick the B call ---------
## The category comes from the hierarchy PATH, not from pattern-matching a label
## string. It is then required to survive an independent marker test -- exactly the
## test the step-03 HCA transfer failed when it called 78% of RB1's strongB cells
## epithelial. Both halves must hold.
bval <- rbindlist(lapply(HIER, function(hn) rbindlist(lapply(have, function(cc) {
  map <- EX[hierarchy == hn & clustering == cc, .(cluster_id, bcat)]
  d   <- merge(ev[, c("cell_id", cc, lin), with = FALSE], map,
               by.x = cc, by.y = "cluster_id", all.x = TRUE)
  d[!is.na(bcat), c(.(hierarchy = hn, clustering = cc, cells = .N),
                    lapply(.SD, function(v) round(mean(v), 4))),
    by = bcat, .SDcols = lin]
}))))
bval[, top_marker := lin[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = lin]
bval[, top_val    := do.call(pmax, .SD), .SDcols = lin]
bval[, second_val := apply(as.matrix(.SD), 1, function(r) sort(r, decreasing = TRUE)[2]),
     .SDcols = lin]
bval[, margin := round(top_val / pmax(second_val, 1e-9), 2)]
bval[, supported := fcase(
  bcat == "B-lineage: B cell", top_marker == "B"  & margin >= B_MARGIN,
  bcat == "B-lineage: plasma", top_marker %in% c("PC", "B"),
  default = NA)]
setcolorder(bval, c("hierarchy","clustering","bcat","cells","top_marker","margin","supported"))
cat("\n=== path-derived categories vs the independent marker evidence\n")
print(bval[order(hierarchy, clustering, bcat)])
sc("hierannot_bcat_marker_validation.csv", bval)

## primary combination: the one recovering the most marker-supported B cells,
## tie-broken by the cleanest margin
cand <- bval[bcat == "B-lineage: B cell" & supported %in% TRUE,
             .(b_cells = sum(cells), margin = min(margin)), by = .(hierarchy, clustering)]
if (!nrow(cand))
  stop("no hierarchy x clustering yields a B-cell category the markers support -- ",
       "do not build a B call on this output", call. = FALSE)
setorder(cand, -b_cells, -margin)
cat("\n=== candidate combinations (marker-supported B-cell category)\n"); print(cand)
BEST_H <- cand$hierarchy[1]; BEST_C <- cand$clustering[1]
cat(sprintf("\nPRIMARY: '%s' on '%s' -- %d B cells, marker margin %.2fx\n",
            BEST_H, BEST_C, cand$b_cells[1], cand$margin[1]))

sup_B <- bval[hierarchy == BEST_H & clustering == BEST_C &
              bcat == "B-lineage: B cell" & supported %in% TRUE, .N] > 0
bmap  <- EX[hierarchy == BEST_H & clustering == BEST_C,
            .(cluster_id, label = get(LAB_COL), annot_path, annot_confidence, bcat)]
CALL  <- merge(ev[, c("cell_id", BEST_C, "nCount", lin), with = FALSE], bmap,
               by.x = BEST_C, by.y = "cluster_id", all.x = TRUE)
setnames(CALL, BEST_C, "cluster_id")
CALL[, `:=`(is_B      = bcat %in% "B-lineage: B cell" & sup_B,
            is_plasma = bcat %in% "B-lineage: plasma",
            hierarchy = BEST_H, clustering = BEST_C)]

## every combination's per-cell call too, so step 08 can test sensitivity to the choice
ALL <- rbindlist(lapply(HIER, function(hn) rbindlist(lapply(have, function(cc) {
  m <- merge(EX[hierarchy == hn & clustering == cc,
                .(cluster_id, label = get(LAB_COL), annot_path, bcat)],
             bval[hierarchy == hn & clustering == cc, .(bcat, supported)],
             by = "bcat", all.x = TRUE)
  d <- merge(ev[, c("cell_id", cc), with = FALSE], m, by.x = cc, by.y = "cluster_id",
             all.x = TRUE)
  data.table(cell_id = d$cell_id, hierarchy = hn, clustering = cc, label = d$label,
             annot_path = d$annot_path, bcat = d$bcat,
             is_B      = d$bcat %in% "B-lineage: B cell" & d$supported %in% TRUE,
             is_plasma = d$bcat %in% "B-lineage: plasma")
}))))

cat("\n=== yield by combination\n")
print(ALL[, .(cells = .N, n_B = sum(is_B), pct_B = round(100*mean(is_B), 1),
              n_plasma = sum(is_plasma), pct_plasma = round(100*mean(is_plasma), 1)),
          by = .(hierarchy, clustering)][order(-pct_B)])
## CONCORDANCE IS MEASURED ON B-LINEAGE, NOT ON is_B. The hierarchies disagree about
## whether a B-marker-dominant cluster's terminal node is "B cell" or "Plasma cell"
## far more than they disagree about WHICH cells are B-lineage, so counting is_B votes
## reads as total disagreement when the combinations actually agree on the cells.
ALL[, blin := bcat != "not B"]
conc <- rbindlist(lapply(c("blin","is_B"), function(v) {
  w  <- dcast(ALL, cell_id ~ hierarchy + clustering, value.var = v)
  kk <- setdiff(names(w), "cell_id")
  w[, n_votes := rowSums(as.matrix(.SD)), .SDcols = kk]
  d  <- w[, .(cells = .N), by = n_votes][order(n_votes)]
  d[, `:=`(measure = v, n_combinations = length(kk))][]
}))
cat("\n=== do the combinations agree on which cells are B-lineage?\n")
print(dcast(conc, n_votes ~ measure, value.var = "cells", fill = 0))
NC <- conc$n_combinations[1]
cat(sprintf("unanimous B-lineage across all %d: %d cells | B-lineage in at least one: %d\n",
            NC, conc[measure=="blin" & n_votes==NC, sum(cells)],
            conc[measure=="blin" & n_votes>0, sum(cells)]))
sc("hierannot_bcall_concordance.csv", conc)

fwrite(CALL[, c("cell_id","cluster_id","label","annot_path","bcat","is_B","is_plasma",
                "annot_confidence","nCount", lin), with = FALSE],
       file.path(OUT, "bcell_call.csv.gz"))
cat("wrote ", file.path(OUT, "bcell_call.csv.gz"), "\n", sep = "")
fwrite(ALL, file.path(OUT, "hierannot_cell_labels.csv.gz"))
cat("wrote ", file.path(OUT, "hierannot_cell_labels.csv.gz"), "\n", sep = "")

## ---- 9. how the new call compares with what step 08 was using --------------
## strongB is the step-01 per-cell heuristic; `lineage` is its argmax version. Both
## are per-cell and depth-limited, which is the whole reason for doing this.
if (all(c("lineage","strongB") %in% colnames(seu@meta.data))) {
  old <- data.table(cell_id = colnames(seu),
                    lineage = as.character(seu$lineage),
                    strongB = seu$strongB %in% TRUE)
  cf <- merge(CALL[, .(cell_id, is_B)], old, by = "cell_id")
  cat("\n=== new cluster-level B call vs the per-cell heuristics\n")
  print(cf[, .(cells = .N, pct_lineage_B = round(100*mean(lineage == "B"), 1),
               pct_strongB = round(100*mean(strongB), 1)), by = is_B][order(-is_B)])
  cat(sprintf("agreement with lineage==B: %.1f%%   with strongB: %.1f%%\n",
              100*mean(cf$is_B == (cf$lineage == "B")), 100*mean(cf$is_B == cf$strongB)))
  sc("bcell_call_vs_percell.csv",
     cf[, .(cells = .N, lineage_B = sum(lineage == "B"), strongB = sum(strongB)), by = is_B])
}

## ---- 10. figures -----------------------------------------------------------
vl <- melt(val, id.vars = c("hierarchy","clustering","label","cells","marker_supports_B"),
           measure.vars = lin, variable.name = "lineage", value.name = "pct")
p1 <- ggplot(vl[hierarchy == BEST_H & clustering == BEST_C],
             aes(reorder(label, cells), pct, fill = lineage)) +
  geom_col(position = "fill", width = .75) + coord_flip() +
  labs(title = "Does the marker evidence back each HierAnnot label?",
       subtitle = sprintf("%s on %s. Bars are the relative share of each lineage's marker set, per label.\nA label accepted as B must have B markers strongest AND clearing %.1fx the runner-up.",
                          BEST_H, BEST_C, B_MARGIN),
       x = NULL, y = "share of lineage marker signal",
       caption = sprintf("Data: %s/hierannot_label_marker_validation.csv", OUT)) +
  theme_viz(xgrid = TRUE)
sf("fig_09_label_marker_validation.png", p1, 10, 7)

yb <- ALL[, .(pct_B = round(100*mean(is_B), 1)), by = .(hierarchy, clustering)]
p2 <- ggplot(yb, aes(hierarchy, pct_B, fill = clustering)) +
  geom_col(position = position_dodge(.8), width = .74) +
  geom_text(aes(label = sprintf("%.1f", pct_B)), position = position_dodge(.8),
            vjust = -0.35, size = 3, colour = INK2) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(title = "B-cell fraction by hierarchy and clustering resolution",
       subtitle = "Marker-supported labels only. Wide disagreement here means the call is resolution-dependent\nand step 08 should carry it as a sensitivity rather than a single number.",
       x = NULL, y = "% of cells called B",
       caption = sprintf("Data: %s/hierannot_cell_labels.csv.gz", OUT)) +
  theme_viz() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
sf("fig_09_bcell_yield.png", p2, 8.5, 5)

cat("\n============================================================\n")
cat(sprintf("PRIMARY B CALL: %s / %s -- %d of %d cells B (%.1f%%), %d plasma (%.1f%%)\n",
            BEST_H, BEST_C, CALL[, sum(is_B)], nrow(CALL), 100*CALL[, mean(is_B)],
            CALL[, sum(is_plasma)], 100*CALL[, mean(is_plasma)]))
cat("Read hierannot_bcat_marker_validation.csv before using this. The category comes\n")
cat("from annot_path, and is accepted only where independent markers back it -- the\n")
cat("label string is never sufficient, which is the test the step-03 transfer failed.\n")
cat("Plasma cells are kept SEPARATE from the B call: they are B-lineage, but they are\n")
cat("exactly what step 07/08 must exclude from a B-cell comparison.\n")
cat("Coverage caveat: this is the POST-QC object, so it reaches ~45% of step 08's\n")
cat("clustered groups but only ~18% of its diffuse group.\n")
cat("============================================================\n")
cat("OK\n")
