#pipeline v2
#this version adheres to bruker's best practices instead of standard seurat guidelines


# Claude Opus 4.8
#potential references:
# PMID 38949890 — Pourmaleki et al., Clin Cancer Res 2024 (the multiplexed HRS-neighborhood paper you just asked about)
# PMID 39896575 — genome-scale spatial mapping of the Hodgkin lymphoma microenvironment (Slide-seqV2; the HRS niche enriched for CD4+ T cells, myeloid, and FDCs)
# The Int J Mol Sci 2026 spatial-transcriptomics review didn't surface a PMID in my earlier search — it came through as an MDPI DOI (10.3390/ijms27083689). If you want it as a PMID I can look it up, or you can drop that one since the two primary papers above carry the neighborhood biology on their own.
#
# For the FL and composite side, if it's useful to have matching PMIDs on hand:
#
#   PMID 35067515 — Composite classic Hodgkin and follicular lymphoma clinicopathologic series (the source for the HRS phenotype frequencies and the "two distinct areas in 95%" point)
# PMID 11474574 — common germinal-center B-cell origin in composite HL + FL (the clonal-relationship basis for your central hypothesis)

# ====== Design Notes =========================================================
# CosMx 6K analysis pipeline — composite classic Hodgkin + follicular lymphoma
# Case study: same patient, PRIMARY (composite cHL + FL) vs RECURRENCE (FL only)
#
# Design note: this is n = 1 (one patient). Treat every result as descriptive /
# hypothesis-generating, not statistically inferential. The one well-controlled
# comparison is FL-niche vs cHL-niche WITHIN the primary (same tissue/batch);
# the primary-vs-recurrence FL comparison is confounded by time and treatment.
#
# HRS (Hodgkin/Reed-Sternberg) calls are heuristic and MUST be confirmed
# visually / with pathology — HRS cells are rare (1-5%), large, multinucleated,
# and easily lost or split by segmentation.

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(HCATonsilData)
  library(hoodscanR)
  library(dplyr)
  library(BPCells)
  library(openxlsx)
})
RDS_PATH = "/netfiles/volaric_research/DLBCL/CosMx/export1/"
# ------------------------------------------------------------------ CONFIG ----
# Edit these to match your object, then source the script section by section.
RDS_PATHS = dir(RDS_PATH, pattern = "RDS$", full.names = TRUE)

# v2
# preserve HRS cells by allowing bcell markers or TNFRSF8
# v3
# cell lineage now assigned using ratios of marker genes

OUT_DIR       <- "output/01_import_v1/"                  # results go here
ASSAY         <- "RNA"                        # CosMx assay name ("RNA" or "Nanostring")
SAMPLE_COL    <- "Run_Tissue_name"                     # colData column identifying the two samples
PRIMARY_ID    <- "IB1"                    # value in SAMPLE_COL for the composite primary
RECURR_ID     <- "RB1"                 # value in SAMPLE_COL for the FL relapse
REF_ANNO_COL  <- "annotation_20230508"        # tonsil-atlas label column (see Section 5)
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- Coarse lineage by marker-ratio argmax --------------------------------
# Cross-platform atlas label transfer (scRNA tonsil -> CosMx 6K) does NOT agree
# with canonical markers, so lineage is called here, directly from markers, using
# the same logic that makes the k:l call trustworthy (within-cell, count-based,
# ratio-driven). Panel-checked against cosmx_genes.Rds; FDCSP dropped (off-panel).
# HRS is deliberately NOT a set here: HRS shed the B program and are lineage-
# negative, so they surface as "unassigned"/"ambiguous" and are called downstream.
# LINEAGE_SETS <- list(
#   B       = c("MS4A1", "CD19", "CD79A", "CD79B", "BANK1"),
#   T       = c("CD3D", "CD3E", "CD3G", "TRAC", "CD2", "CD5"),
#   NK      = c("NKG7", "GNLY", "KLRD1", "NCAM1", "KLRB1", "NCR1"),
#   myeloid = c("CD68", "CD163", "LYZ", "ITGAX", "C1QA", "C1QB", "CD14"),
#   FDC     = c("CR2", "CR1", "FCER2", "CXCL13"),
#   endo    = c("PECAM1", "VWF", "CLDN5", "CD34"),
#   fibro   = c("COL1A1", "COL1A2", "PDGFRB", "DCN", "ACTA2"),
#   PC      = c("MZB1", "XBP1", "PRDM1", "SDC1", "DERL3", "TNFRSF17"))  # secretory, NOT Ig

LINEAGE_SETS <- list(
  B       = c("MS4A1", "CD79A", "CD79B", "BANK1"),          # drop CD19 (dropout)
  T       = c("CD3D", "CD3E", "CD3G", "TRAC"),              # drop CD2/CD5 (NK/B overlap)
  NK      = c("NKG7", "GNLY", "KLRD1", "NCR1"),             # drop NCAM1/KLRB1 (T overlap)
  myeloid = c("CD68", "CD163", "LYZ", "ITGAX", "C1QA", "C1QB", "CD14"),
  endo    = c("PECAM1", "VWF", "CLDN5", "CD34"),
  fibro   = c("COL1A1", "COL1A2", "PDGFRB", "DCN", "ACTA2"),
  PC      = c("MZB1", "XBP1", "PRDM1", "SDC1", "DERL3", "TNFRSF17"))
# FDC dropped: CR2/FCER2 are B markers; nothing FDC-exclusive survives on this panel


LINEAGE_MARGIN  <- 1      # log2 top:second set score to ASSIGN (else "ambiguous"); 1 = 2:1
LINEAGE_BG_MULT <- 2      # winning set's raw mean must beat BG_MULT x per-cell background
LINEAGE_MIN_RAW <- 0.2    # ...and an absolute per-gene floor (guards cells w/ 0 negprobes)
SCALE_CAP       <- 10     # winsorize per-gene fold-over-mean so one bright gene can't rule

RDS_PATH = RDS_PATHS[1]
all_seu = list()

.detect_gl_any = function(M, gl, req_gene_count = 1, min_count = 1){
  if(!all(gl %in% rownames(M))){
    message(setdiff(gl, rownames(M)))
    stop("above genes are missing")
  }
  det = (M[gl, , drop = FALSE] %>% as.matrix) >= min_count
  colSums(det) >= req_gene_count
}

.detect_gl_all = function(M, gl){
  .detect_gl_any(M, gl, req_gene_count = length(gl))
}

.goi_counts = function(M, goi){
  if(!goi %in% rownames(M)) stop(goi, " missing")
  bc[goi, ] %>% as.matrix %>% as.vector
}

# Coarse-lineage assignment (the k:l ratio logic, generalized to N marker sets).
# Per gene: fold-over-mean, winsorized at `cap`, so sets of differing brightness
# are comparable. Per set: mean of the scaled values (size-normalized). A cell is
# ASSIGNED its argmax set only if that set (a) clears background in RAW counts --
# else "unassigned" -- and (b) beats the runner-up by `margin` log2 -- else
# "ambiguous". The raw floor is per-cell (negative-probe rate), so spillover into
# an intrinsically-dim lineage cannot mint a call. `cts` = counts (BPCells ok),
# `neg_total` = per-cell negprobe totals, `n_neg` = number of negative probes.
assign_lineage <- function(cts, sets, neg_total, n_neg,
                           margin = 1.5, bg_mult = 1.5,
                           min_raw = .2, cap = 10, use_max_marker = TRUE) {
  sets <- lapply(sets, function(g) intersect(g, rownames(cts)))
  sets <- sets[lengths(sets) > 0]
  M    <- as.matrix(cts[unique(unlist(sets)), , drop = FALSE])   # markers x cells
  # gmean<- pmax(rowMeans(M), 1e-3)
  # browser()
  gmean <- rowSums(M) / pmax(rowSums(M > 0), 1)   # mean among detected cells
  gmean <- pmax(gmean, 1e-3)
  S    <- M / gmean;
  S[S > cap] <- cap                           # fold-over-mean, winsorized

  marker_agg_FUN = if(use_max_marker){
    colMaxs
  }else{
    colMeans
  }

  scaled <- sapply(sets, function(g) marker_agg_FUN(S[g, , drop = FALSE]))  # cells x sets
  raw    <- sapply(sets, function(g) marker_agg_FUN(M[g, , drop = FALSE]))
  bg     <- as.numeric(neg_total) / max(n_neg, 1)                # per-cell background

  ii     <- seq_len(nrow(scaled))
  top_i  <- max.col(scaled, ties.method = "random")               # argmax (vectorized)
  s2     <- scaled; s2[cbind(ii, top_i)] <- -Inf
  sec_i  <- max.col(s2, ties.method = "random")
  top_s  <- scaled[cbind(ii, top_i)]; sec_s <- s2[cbind(ii, sec_i)]
  top_raw<- raw[cbind(ii, top_i)]

  pass   <- top_raw > pmax(bg_mult * bg, min_raw)                # real signal vs background

  l2m    <- log2((top_s + 1e-6) / (sec_s + 1e-6))
  lab    <- colnames(scaled)[top_i]
  lab_sec = colnames(scaled)[sec_i]


  lab[!pass]               <- "below_background"
  lab_sec[!pass]               <- "below_background"
  all_sets_0 = rowMax(raw) == 0
  lab[all_sets_0] = "all_markers_0"
  lab_sec[all_sets_0] = "all_markers_0"
  table(lab)

  lab.raw = lab

  lab[pass & l2m < log2(margin)] <- "ambiguous"
  table(lab)
  list(label = lab, log2_margin = l2m, raw_label = lab.raw, sec_label = lab_sec)
}

RDS_PATH = RDS_PATHS[1]
for(RDS_PATH in RDS_PATHS){
  message(RDS_PATH)
  out_root = file.path(OUT_DIR, basename(RDS_PATH) %>% sub(".RDS$", "", .))
  sample_id = strsplit(basename(RDS_PATH), "[_\\.]")[[1]][2]


  FULL_RDS = paste0(out_root, "_BPCells_full.Rds")
  if(file.exists(FULL_RDS)){
    message("skip cached full rds")
    next
  }
  message("starting from initial complete seurat object (high memory)")
  #### Shrink full object  ####
  seu <- readRDS(RDS_PATH)
  seu@assays$Volaric_02032026_Normalization.RNA.1_1 = NULL
  gc()

  ### BPCells ####
  BP_DIR = paste0(out_root, "_BPCells_count_matrix")

  seu %>% object.size() %>% format(units = "GB")
  seu@assays %>% object.size() %>% format(units = "GB")
  seu@assays$RNA %>% object.size() %>% format(units = "GB")
  seu@assays$RNA@meta.features %>% object.size() %>% format(units = "GB")
  seu@assays$RNA@data[1:20, 1:10]
  seu@assays$RNA@counts %>% object.size() %>% format(units = "GB")
  seu@assays$RNA@counts[1:20, 1:10]
  # Write the counts layer to a directory
  if(!file.exists(BP_DIR)) write_matrix_dir(mat = seu@assays$RNA@counts  %>% convert_matrix_type(type = "uint32_t"), dir = BP_DIR)

  # replace in-memory counts
  counts.mat <- open_matrix_dir(dir = BP_DIR)

  stopifnot(all.equal(colnames(counts.mat), colnames(seu)))
  seu[["RNA"]] <- CreateAssay5Object(counts = counts.mat)

  seu = UpdateSeuratObject(seu)


  seu@assays$falsecode %>% object.size() %>% format(units = "GB")
  seu@assays$negprobes %>% object.size() %>% format(units = "GB")
  seu@assays$RNA %>% object.size() %>% format(units = "GB")

  #### Coarse lineage (marker-ratio argmax) ####
  # Computed on the FULL object so fold-over-mean scaling uses whole-tissue gene
  # means; the label rides into seu_b (and the step-10 merge) via subsetting.
  cts    <- LayerData(seu, layer = "counts")
  genes = rownames(seu)

  neg_mat = seu@assays$negprobes@counts
  lin <- assign_lineage(
    cts, LINEAGE_SETS,
    neg_total = colSums(neg_mat), n_neg = nrow(neg_mat),
    margin = 1.5, bg_mult = 1.5,
    min_raw = .2, cap = 10,
    use_max_marker = FALSE
  )
  label_tab = lin$label %>% table
  sec_tab = table(lin$raw_label, lin$sec_label)
  (sec_tab/1e3) %>% round
  #columns are sec, rows are first
  sec_tab.norm = apply(sec_tab, 1, function(x){
    x / sum(x)
  }) %>% t
  sec_tab.norm
  cbind(label_tab)
  neg_cat = c(
    "all_markers_0",
    "below_background",
    "ambiguous"
  )
  ct_tab = label_tab[setdiff(names(label_tab), neg_cat)]
  ct_frac = ct_tab / sum(ct_tab)
  label_tab["B"]
  ct_frac["B"]

  seu$per_cell_marker_lineage = lin$label

  b_markers <- c("MS4A1","CD79A","CD79B","BANK1")
  Bm      <- as.matrix(cts[b_markers, ])
  strongB <- colSums(Bm >= 1) >= 2 & colSums(Bm) >= 3   # >=2 B markers co-detected
  sum(strongB)/sum(ct_tab)

  seu$strongB = strongB

  seu$lineage        <- factor(lin$label)
  seu$lineage_margin <- lin$log2_margin
  message(sample_id, " coarse lineage:")
  print(sort(table(seu$lineage), decreasing = TRUE))

  old_cn = colnames(seu@meta.data)
  cn = old_cn
  seu$insitutype_deriveProfiles %>% table
  k = grepl("Volaric_02032026_Neighbor.network.expression.space.1_1_cluster_Volaric_02032026_Leiden.Clustering.1_1", cn)
  cn[k] = sub("Volaric_02032026_Neighbor.network.expression.space.1_1_cluster_Volaric_02032026_Leiden.Clustering.1_1", "neighbor_leiden", cn[k])

  cn
  k = grepl("RNA_Volaric_02032026_Cell.Typing.InSituType.1_1", cn)
  cn[k] = sub("RNA_Volaric_02032026_Cell.Typing.InSituType.1_1", "inSitu_clusters", cn[k]) %>% sub("clusters_clusters", "clusters", .)

  k = grepl("RNA_Volaric_02032026_Neighborhood.Analysis.1_1_", cn)
  cn[k] = sub("RNA_Volaric_02032026_Neighborhood.Analysis.1_1_", "inSitu_", cn[k])

  k = grepl("Volaric_02032026_Neighborhood.Analysis.1_1_assignments", cn)
  cn[k] = sub("Volaric_02032026_Neighborhood.Analysis.1_1_assignments", "inSitu_neighborhood_assignment", cn[k])

  colnames(seu@meta.data) = cn

  seu@misc$metadata_renames =  data.frame(
    orginal = old_cn[cn != old_cn],
    rename = cn[cn != old_cn]
  )

  saveRDS(seu, FULL_RDS)
}



