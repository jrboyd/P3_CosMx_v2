## ============================================================================
## 09b_script_marker_relabel.R -- audit HierAnnot against markers, and relabel.
##
## WHY THIS STEP EXISTS. Step 09 runs HierAnnot over pearson_clusters_low/high. Reading
## its calls against a canonical marker panel showed four of nine `pearson_clusters_low`
## labels are wrong, and wrong in ways that would invert an interpretation:
##
##   cl 2  "Myeloid > Macrophage"  ->  stromal/FDC   (PDGFRB, COL1A1, ACTA2, CLU, CR2)
##   cl 5  "Neutrophil"            ->  endothelial   (CDH5, VWF, PECAM1)
##   cl 6  "B cell > Plasma cell"  ->  pDC           (LILRA4, IRF7)
##   cl 7  "Unresolved"            ->  cDC           (BATF3, FCER1A)
##
## and cluster 0 -- 45% of the node -- is called "Mast cell", which no lymph node can be.
## The other hierarchies are worse: `tme_core` and `tonsil_tme` both call that cluster
## "Epithelial > Luminal epithelial" in a lymph node, despite outranking `immune_core` on
## HierAnnot's own fit_rank_score. The fit score does not track correctness here, so the
## audit cannot be skipped by picking a better-scoring hierarchy.
##
## WHAT THIS DOES. Scores every cluster against a curated lymph-node marker panel,
## assigns the best-supported family, and records WHY -- the margin over the runner-up,
## how much of the cluster's signal any marker explains, and whether the assigned lineage
## can plausibly be that abundant. Output is the canonical label table that downstream
## scripts read through scripts/lib/cluster_labels.R.
##
## THIS IS AN AUDIT, NOT AN ORACLE. It replaces one automated call with another. What it
## adds is that every call is accompanied by its evidence and a status, so a weak call
## announces itself instead of propagating silently. Clusters that markers cannot resolve
## are labelled `low_information` rather than forced into a lineage.
##
## Usage: Rscript scripts/09b_script_marker_relabel.R [SAMPLE]
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(data.table); setDTthreads(1)
})

args    <- commandArgs(trailingOnly = TRUE)
SAMPLE  <- if (length(args) >= 1) args[1] else "IB1"
## "merged" is the step-04 object, which carries its own pearson_clusters_* (7 low, 26
## high) alongside the per-sample ids in sample_pearson_clusters_*. Note it is built on the
## STRICT tier (IB1 400,111 + RB1 1,399,189), so its cluster sizes are not comparable to
## the matched-tier counts used elsewhere in the project.
IN_RDS  <- if (SAMPLE == "merged") "output/04_merge/merged_processed.Rds" else
           sprintf("output/03_filter_and_process/%s_03_HCA_cell_typed.Rds", SAMPLE)
LAB_F   <- sprintf("output/09_hierannot/%s/hierannot_cluster_labels.csv", SAMPLE)
OUT     <- file.path("output/09b_marker_relabel", SAMPLE)
FIG     <- file.path(OUT, "figures")
## inSitu_clusters is the VENDOR InSituType call, ids "a".."j". It is the only clustering
## whose levels mean the same thing in both samples and in the merged object, because the
## vendor assigned them once rather than re-deriving them per object -- which makes it the
## only one where "does IB1's cluster c mean what RB1's cluster c means" is directly
## answerable. Its ids are LETTERS, which is why cluster_id is character throughout below.
CLUSTERINGS <- c("pearson_clusters_low", "pearson_clusters_high", "inSitu_clusters")
DEFECT_CL   <- "i"

## Cluster ids are character so letter ids survive, but a plain sort() then orders 26
## clusters "0,1,10,11,...,2,20" and silently mis-columns the heatmap. Order numerically
## where every id is a number, lexically otherwise.
cid_order <- function(x) {
  x <- as.character(x)
  if (all(grepl("^[0-9]+$", x))) x[order(as.integer(x))] else sort(x)
}

## Thresholds for the status call. Deliberately explicit rather than tuned: each one
## encodes a stated belief about what "resolved" means, and each is reported per cluster
## so a reader can disagree with the threshold and still use the numbers.
MIN_TOP_Z    <- 0.75  # top PLAUSIBLE family must stand this far above the cross-cluster mean
MIN_STRONG_Z <- 1.50  # ...and this far to count as strongly rather than weakly supported.
                      # Without this tier, cluster 0 passed as "confident" on top_z 1.04
                      # while every other confident call scored 2.0-2.6 -- a real gap that
                      # a pass/fail threshold hides. `weak_support` names it instead.
MIN_MARGIN   <- 0.30  # ...and beat the runner-up family by this much
DEFECT_MAX   <- 0.50  # a cluster more than half over-fixation defect gets no lineage label
MIN_FAM_GENES <- 2L   # a family with fewer genes than this on the panel cannot be
                      # `confident`: the z-per-gene-then-average safeguard needs >=2 genes
                      # to mean anything. Only Neutrophil (1 gene) is affected.
MIN_MK_SHARE <- 0.010 # ...and lineage markers must be >=1% of the cluster's transcripts.
                      # This is a floor for a marker-dead cluster, not a discriminator:
                      # observed values run 1.4-2.7%, so it fires only on something worse
                      # than anything currently in the data.

## Sub-lineages and their parents. When the top two families are a child and its own
## parent within MIN_MARGIN, the defensible call is the PARENT -- reporting "Tfh" over
## "T cell" on a margin of 0.141 claims a subset resolution the data does not support.
PARENT <- c(`Tfh` = "T cell", `Treg` = "T cell", `GC B cell` = "B cell")

## Curated lymph-node panel. `cap` is the largest share of a node this lineage can
## credibly occupy; exceeding it flags `implausible_abundance` rather than silently
## accepting that 45% of a lymph node is mast cells. Caps are generous on purpose --
## they exist to catch absurdity, not to enforce an expected composition.
FAMILIES <- list(
  `B cell`            = list(cap = .70, g = c("MS4A1","CD19","CD79A","CD79B","PAX5","BANK1","CD22","FCER2","CR2")),
  `GC B cell`         = list(cap = .40, g = c("BCL6","AICDA","RGS13","LMO2","MEF2B","MKI67")),
  `Plasma cell`       = list(cap = .30, g = c("MZB1","XBP1","PRDM1","DERL3","TNFRSF17","JCHAIN")),
  `T cell`            = list(cap = .70, g = c("CD3D","CD3E","CD3G","TRAC","CD2","IL7R","CD8A","CD8B","CD4")),
  `Treg`              = list(cap = .20, g = c("FOXP3","IL2RA","CTLA4")),
  `Tfh`               = list(cap = .25, g = c("CXCL13","PDCD1","ICOS")),
  `NK cell`           = list(cap = .15, g = c("NKG7","GNLY","KLRD1","NCAM1","KIR2DL1/3")),
  `Macrophage`        = list(cap = .40, g = c("CD68","CD163","CSF1R","C1QA","C1QB","MRC1","MSR1","AIF1")),
  `Monocyte/myeloid`  = list(cap = .35, g = c("LYZ","ITGAX","S100A8","IRF8","FLT3","CD14")),
  `cDC`               = list(cap = .12, g = c("BATF3","CLEC9A","XCR1","CD1C","FCER1A","LAMP3","CCR7")),
  `pDC`               = list(cap = .08, g = c("LILRA4","IRF7","IL3RA")),
  `Mast cell`         = list(cap = .06, g = c("TPSAB1","CPA3","KIT")),
  ## S100A9 removed: monocytes and macrophages express it heavily, so it is a myeloid
  ## marker rather than a neutrophil one. With it in, cluster 1 -- unambiguously
  ## macrophage by CD68/CSF1R/C1QA/CD163 -- was called Neutrophil by a margin of 0.189.
  `Neutrophil`        = list(cap = .12, g = c("FCGR3B","CSF3R")),
  `Endothelial`       = list(cap = .18, g = c("PECAM1","VWF","CDH5","CLDN5")),
  `Stromal/FDC`       = list(cap = .25, g = c("COL1A1","PDGFRB","ACTA2","FDCSP","CLU","VCAM1","C3","CCL21","SERPINE2")),
  `Epithelial`        = list(cap = .10, g = c("KRT5","KRT14","EPCAM","KRT19")))

dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
## Value scales run blue -> red (CLAUDE.md): high family score is red.
DIV_HI<-"#eb6834"; DIV_LO<-"#2a78d6"
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
        legend.position="top", legend.justification="left")
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF,limitsize=FALSE)
                       cat("wrote ",file.path(FIG,n),"\n",sep="")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote ",file.path(OUT,n),"\n",sep="")}

cat("============================================================\n")
cat("step 09b -- marker audit of HierAnnot labels, sample ", SAMPLE, "\n", sep="")
cat("============================================================\n")

if (!file.exists(IN_RDS)) stop("missing ", IN_RDS, call. = FALSE)
obj <- readRDS(IN_RDS); rn <- rownames(obj@meta.data)
need <- c(CLUSTERINGS, "inSitu_clusters", "nCount_RNA")
miss <- setdiff(need, colnames(obj@meta.data))
if (length(miss)) stop("meta.data missing: ", paste(miss, collapse=", "), call.=FALSE)
md <- as.data.table(obj@meta.data[, need, drop = FALSE]); md[, cell_id := rn]
cnt <- LayerData(obj[["RNA"]], layer = "counts")
GEN <- unique(unlist(lapply(FAMILIES, `[[`, "g")))
G   <- intersect(GEN, rownames(cnt))
absent <- setdiff(GEN, G)
cat(sprintf("\nmarker panel: %d of %d genes on the 6K panel\n", length(G), length(GEN)))
if (length(absent)) cat("  absent, dropped: ", paste(absent, collapse=", "), "\n", sep="")
M <- as(cnt[G, , drop = FALSE], "dgCMatrix"); colnames(M) <- rn
tot_all <- md$nCount_RNA
rm(obj, cnt); invisible(gc())

## Defect cells are excluded from scoring. Over-fixation flattens a profile toward
## uniform, which is precisely the distortion that would blur a lineage call.
ok <- md$inSitu_clusters != DEFECT_CL & md$nCount_RNA > 0
cat(sprintf("scoring on %d of %d cells (%d dropped: cluster '%s' or zero counts)\n",
            sum(ok), nrow(md), sum(!ok), DEFECT_CL))
Mk <- M[, ok, drop = FALSE]; mdk <- md[ok]; totk <- tot_all[ok]

## Defect load per cluster, measured BEFORE the exclusion above. A cluster can be mostly
## cluster `i` and still produce a clean-looking score, because scoring only ever sees the
## surviving minority. pearson_clusters_low cluster 0 is the case that matters: 70% of it
## is defect, and it holds 86% of every cluster-`i` cell in IB1 -- it is where the
## over-fixation population landed in this clustering, not a lineage that scored badly.
## Over-fixation flattens a profile toward uniform (see the C2 diagnosis), so a lineage
## call on the remnant is not meaningful however well it scores.
defect_frac <- lapply(setNames(CLUSTERINGS, CLUSTERINGS), function(cc) {
  t <- md[, .(n = .N, n_def = sum(inSitu_clusters == DEFECT_CL)), by = c(cc)]
  setNames(t$n_def / t$n, as.character(t[[cc]]))
})
cat("\n=== fraction of each cluster that is over-fixation cluster '", DEFECT_CL, "'\n", sep="")
for (cc in CLUSTERINGS) {
  v <- sort(defect_frac[[cc]], decreasing = TRUE)
  cat(" ", cc, ": ", paste(sprintf("%s=%.0f%%", names(v), 100*v), collapse = "  "), "\n", sep="")
}
FAM_OK <- lapply(FAMILIES, function(f) { f$g <- intersect(f$g, G); f })
FAM_OK <- FAM_OK[vapply(FAM_OK, function(f) length(f$g) > 0, TRUE)]

## the family's own genes, ordered by how strongly they mark THIS cluster.
## setNames is load-bearing: Z[genes, k] on a ONE-gene family drops to an unnamed scalar,
## so names(sort(v)) was NULL and top_markers came out EMPTY. Every Neutrophil call in the
## project (the only single-gene family) had a blank evidence column because of this.
f_top <- function(Z, genes, k) {
  v <- setNames(as.numeric(Z[genes, k]), genes)
  names(sort(v, decreasing = TRUE))
}

score_one <- function(clustering) {
  cl <- as.character(mdk[[clustering]])
  KF <- factor(cl)
  IND <- sparseMatrix(i = seq_along(cl), j = as.integer(KF), x = 1,
                      dims = c(length(cl), nlevels(KF)), dimnames = list(NULL, levels(KF)))
  PB  <- as.matrix(Mk %*% IND)
  ## Denominator is TOTAL transcripts per cluster, not marker transcripts, so
  ## marker_share below means "share of the cluster's whole profile".
  cl_tot <- as.numeric(tapply(totk, KF, sum)[colnames(PB)])
  CP  <- sweep(PB, 2, pmax(cl_tot,1), "/") * 1e4
  ## z per GENE across clusters: each marker votes on which cluster it distinguishes,
  ## which removes the between-gene abundance differences that would otherwise let one
  ## highly-expressed marker decide a family on its own.
  Z <- t(scale(t(log1p(CP)))); Z[!is.finite(Z)] <- 0
  sizes <- as.numeric(table(KF)[colnames(CP)]); names(sizes) <- colnames(CP)
  FS <- t(vapply(FAM_OK, function(f) colMeans(Z[f$g, , drop=FALSE]), numeric(ncol(Z))))
  rownames(FS) <- names(FAM_OK)
  mk_share <- colSums(CP) / 1e4

  res <- rbindlist(lapply(colnames(FS), function(k) {
    v <- sort(FS[, k], decreasing = TRUE)
    frac <- unname(sizes[k]) / sum(sizes)
    uncon <- names(v)[1]                      # what scored highest, ignoring plausibility

    ## The abundance cap is a CONSTRAINT, not an annotation. A previous version merely
    ## flagged cluster 0 as "Mast cell (implausible_abundance)" at 21.5% of the node and
    ## left the wrong label in place for anything downstream that ignored the flag.
    ## Restrict to families that can credibly be this abundant and pick the best of them.
    ok_fam <- names(v)[vapply(names(v), function(f) frac <= FAM_OK[[f]]$cap, TRUE)]
    if (!length(ok_fam)) ok_fam <- names(v)      # nothing plausible; fall back, flag below
    top <- ok_fam[1]; tz <- unname(v[top])
    rest <- setdiff(ok_fam, top)
    run  <- if (length(rest)) rest[1] else NA_character_
    rz   <- if (length(rest)) unname(v[run]) else NA_real_
    mrg  <- if (length(rest)) tz - rz else NA_real_

    ## Child-vs-parent ties resolve to the parent rather than overclaiming a subset.
    note <- NA_character_
    if (!is.na(mrg) && mrg < MIN_MARGIN && !is.na(PARENT[top]) &&
        identical(unname(PARENT[top]), run)) {
      note <- paste0("resolved to parent lineage; ", top, " scored ", round(tz,2),
                     " vs ", round(rz,2), " -- subset call not supported")
      top <- run; tz <- rz
      rest <- setdiff(ok_fam, c(top, names(PARENT)[PARENT == top]))
      run <- if (length(rest)) rest[1] else NA_character_
      rz  <- if (length(rest)) unname(v[run]) else NA_real_
      mrg <- if (length(rest)) tz - rz else NA_real_
    }

    dfr <- unname(defect_frac[[clustering]][k]); if (is.na(dfr)) dfr <- 0
    st <- if (dfr > DEFECT_MAX)                             "defect_dominated"
          else if (mk_share[k] < MIN_MK_SHARE || tz < MIN_TOP_Z) "low_information"
          else if (frac > FAM_OK[[top]]$cap)               "implausible_abundance"
          else if (!is.na(mrg) && mrg < MIN_MARGIN)        "ambiguous"
          else if (tz < MIN_STRONG_Z)                      "weak_support"
          ## A family with one gene on the panel cannot be `confident`. The whole design of
          ## the family score -- z per gene, then averaged within family -- exists so that
          ## "one abundant marker cannot carry a call". With a single gene the family score
          ## IS that gene's z, so the safeguard does not apply and the call must not claim
          ## the top tier. Neutrophil is the only such family here (CSF3R is on the 6K
          ## panel, FCGR3B is not), and it was returning `confident` on two clusters.
          else if (length(FAM_OK[[top]]$g) < MIN_FAM_GENES) "weak_support"
          else "confident"
    if (st == "defect_dominated") note <- sprintf(
      "%.0f%% of this cluster is over-fixation cluster '%s'; scored on the %.0f%% remnant, best was %s (%.2f) -- %s",
      100*dfr, DEFECT_CL, 100*(1-dfr), top, tz,
      if (tz < MIN_STRONG_Z) "remnant is weak too, so no lineage call"
      else "remnant is strong, so the label stands but the cluster is defect-loaded")
    if (st == "low_information") note <- paste(na.omit(c(note,
      sprintf("best plausible family %s scores only %.2f", top, tz))), collapse = "; ")
    if (length(FAM_OK[[top]]$g) < MIN_FAM_GENES) note <- paste(na.omit(c(note, sprintf(
      "%s has only %d gene on this panel, so the call rests on a single marker and is capped below confident",
      top, length(FAM_OK[[top]]$g)))), collapse = "; ")
    if (!identical(top, uncon) && is.na(note))
      note <- sprintf("%s scored higher (%.2f) but exceeds its %.0f%% abundance cap",
                      uncon, unname(v[1]), 100*FAM_OK[[uncon]]$cap)

    data.table(clustering = clustering, cluster_id = k, cells = unname(sizes[k]),
               pct_cells = round(100*frac, 2),
               ## A defect-dominated cluster still gets its lineage label IF the clean
               ## remnant speaks clearly. high cluster 8 is 52% defect yet scores 3.97
               ## with a 3.06 margin -- the strongest call in the audit -- and blanking
               ## that to "low-information" would discard real evidence. Only a
               ## defect-dominated cluster whose remnant is ALSO weak loses its label.
               label = if (st == "low_information" ||
                           (st == "defect_dominated" && tz < MIN_STRONG_Z))
                         "low-information" else top,
               pct_defect = round(100*dfr, 1),
               best_supported = top, unconstrained_top = uncon,
               top_z = round(tz,3), runner_up = run, runner_up_z = round(rz,3),
               margin = round(mrg,3), marker_share = round(unname(mk_share[k]), 5),
               abundance_cap = FAM_OK[[top]]$cap, n_family_genes = length(FAM_OK[[top]]$g),
               status = st, note = note,
               top_markers = paste(head(f_top(Z, FAM_OK[[top]]$g, k), 5), collapse=", "))
  }))
  list(res = res, Z = Z, CP = CP, FS = FS)
}

hier <- if (file.exists(LAB_F)) {
  L <- fread(LAB_F, select = c("cluster_id","annot_path","clustering","hierarchy"))
  dcast(L, cluster_id + clustering ~ hierarchy, value.var = "annot_path")
} else { cat("\nNOTE: no HierAnnot labels for ", SAMPLE, " -- emitting marker calls only\n", sep=""); NULL }

ALL <- list()
for (cc in CLUSTERINGS) {
  cat("\n============================================================\n")
  cat(cc, "\n"); cat("============================================================\n")
  o <- score_one(cc); r <- o$res
  ## cluster_id stays CHARACTER. It used to be coerced with as.integer() here, which
  ## returns NA for every row of a letter-id clustering -- by coercion, not by error, so
  ## an inSitu_clusters audit would have completed and written a file full of NA ids.
  r[, cluster_id := as.character(cluster_id)]

  ## Scoring runs on the defect-FREE subset (`ok`, above), so when the clustering IS
  ## inSitu_clusters the defect cluster has no cells left and drops out of the result
  ## entirely -- a 10-level clustering would quietly return 9 rows. Put it back explicitly,
  ## labelled and with the circularity named: on this clustering "100% defect" is the
  ## DEFINITION of cluster `i`, not evidence about it. The over-fixation characterisation
  ## is step 07's (output/07_overfixation/), and this row must never be cited for it.
  if (cc == "inSitu_clusters" && !DEFECT_CL %in% r$cluster_id) {
    n_def <- md[get(cc) == DEFECT_CL, .N]
    r <- rbindlist(list(r, data.table(
      clustering = cc, cluster_id = DEFECT_CL, cells = n_def,
      pct_cells = round(100 * n_def / nrow(md), 2),
      label = "low-information", status = "defect_dominated",
      pct_defect = 100, best_supported = NA_character_, unconstrained_top = NA_character_,
      note = paste("the over-fixation cluster itself: excluded from scoring because scoring",
                   "drops defect cells, so 100% defect here is definitional and is NOT",
                   "independent evidence -- see step 07"))), fill = TRUE)
    cat(sprintf("\nnote: cluster '%s' (%d cells) re-inserted unscored -- it defines the defect\n",
                DEFECT_CL, n_def))
  }

  if (!is.null(hier)) {
    h <- hier[clustering == cc, .(cluster_id, hierannot_immune_core = immune_core,
                                  hierannot_tme_core = tme_core,
                                  hierannot_tonsil_tme = tonsil_tme)]
    h[, cluster_id := as.character(cluster_id)]   # fread gives integer for numeric ids
    r <- merge(r, h, by = "cluster_id", all.x = TRUE)
    ## Agreement is judged on the LINEAGE ROOT, not the exact string: HierAnnot writes
    ## "Immune > Myeloid > Macrophage" where this panel writes "Macrophage", and calling
    ## that a disagreement would drown the real errors in formatting noise.
    root <- function(x) sub(".*> *", "", x)
    r[, agrees := {
        h <- hierannot_immune_core
        label <- best_supported
        ok <- !is.na(h)
        out <- rep(NA, .N)
        out[ok] <- mapply(function(a, b)
            isTRUE(a == b) || isTRUE(grepl(a, b, fixed=TRUE)) || isTRUE(grepl(b, a, fixed=TRUE)),
            root(h[ok]), label[ok])
        out
      }]
  } else r[, agrees := NA]
  r <- r[match(cid_order(r$cluster_id), r$cluster_id)]
  print(r[, .(cluster_id, pct = pct_cells, pct_def = pct_defect, label, status, top_z, margin,
              hierannot = if (is.null(hier)) NA_character_ else hierannot_immune_core, agrees)])
  if (any(!is.na(r$note))) {
    cat("\nnotes:\n"); for (i in which(!is.na(r$note)))
      cat(sprintf("  cl %-3s %s\n", r$cluster_id[i], r$note[i]))
  }
  ALL[[cc]] <- r
  hm <- as.data.table(as.table(o$FS)); setnames(hm, c("family","cluster_id","z"))
  hm[, cluster_id := factor(cluster_id, levels = cid_order(unique(r$cluster_id)))]
  lim <- max(abs(hm$z))
  lv <- r[, setNames(sprintf("%s\n%s%s", cluster_id, label,
            ifelse(status=="confident","",paste0("\n[",status,"]"))), as.character(cluster_id))]
  p <- ggplot(hm, aes(cluster_id, family, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = DIV_LO, mid = "#f7f6f2", high = DIV_HI, midpoint = 0,
                         limits = c(-lim, lim), name = "family score (mean z of member markers)") +
    scale_x_discrete(labels = lv, expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
    labs(title = sprintf("Marker audit of %s (%s)", cc, SAMPLE),
         subtitle = paste0("Each gene is z-scored across clusters, then averaged within its ",
           "family, so one abundant marker cannot carry a call.\nColumn labels are the ",
           "assigned label and, where the call is not clean, why."),
         x = NULL, y = NULL,
         caption = sprintf("Data: %s/cluster_labels_final.csv", OUT)) +
    theme_viz(10) + theme(panel.grid = element_blank(),
                          axis.text.x = element_text(size = 7, lineheight = .95))
  ## Width floor of 11 so the subtitle is not clipped; the low clustering has only 9
  ## columns and .8*9 was too narrow to fit the caption line.
  sf(sprintf("fig_09b_family_scores_%s.png", sub("pearson_clusters_","",cc)),  # inSitu keeps its full name
     p, max(11, .85*nrow(r)), 6.5)
}

FIN <- rbindlist(ALL, fill = TRUE)
setcolorder(FIN, c("clustering","cluster_id","label","status","cells","pct_cells"))
sc("cluster_labels_final.csv", FIN)

cat("\n============================================================\n")
cat("AUDIT SUMMARY\n============================================================\n")
print(FIN[, .(clusters = .N, confident = sum(status=="confident"),
              weak_support = sum(status=="weak_support"),
              ambiguous = sum(status=="ambiguous"),
              low_information = sum(status=="low_information"),
              defect_dominated = sum(status=="defect_dominated"),
              implausible = sum(status=="implausible_abundance")), by = clustering])
cat("\n=== defect load: clusters holding the over-fixation population\n")
print(FIN[order(-pct_defect)][1:6, .(clustering, cluster_id, label, status,
      pct_cells, pct_defect)])
if (!is.null(hier)) {
  print(FIN[, .(clusters = .N, agree_with_hierannot = sum(agrees %in% TRUE),
                disagree = sum(agrees %in% FALSE)), by = clustering])
  cat("\n=== where the marker panel OVERRULES HierAnnot\n")
  print(FIN[agrees %in% FALSE, .(clustering, cluster_id, pct_cells,
        hierannot = hierannot_immune_core, markers_say = label, status, top_markers)])
}
cat("\nDownstream scripts should read these via scripts/lib/cluster_labels.R,\n",
    "not output/09_hierannot/*/hierannot_cluster_labels.csv.\n", sep="")
cat("============================================================\nOK\n")
