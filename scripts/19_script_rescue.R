## ============================================================================
## 19_script_rescue.R -- label the cells OUTSIDE the high-confidence object.
##
## The high-confidence set is 1,993,349 cells; the sample holds 4,406,940. The other
## 2,413,591 are defect, failed the symmetric QC tier, or fell to the per-cell
## low-information floor. They carry no cluster label -- and the spatial maps need them,
## because a tissue map drawn from 45% of the cells is not a tissue map.
##
## TWO COLUMNS, NEVER MERGED.
##   celltype_argmax     per-cell marker evidence. Already available for all 4.4M cells from
##                       step 16; tier-independent, unaffected by any of this.
##   celltype_projected  the cluster label, carried over by projecting each cell into the
##                       high-confidence PCA and taking a kNN vote among reference cells.
## Keeping them apart is the whole design: a figure can use either, and no statement can
## silently mix marker evidence with guilt by association.
##
## EXPRESSION SPACE ONLY. NEVER SPATIAL. A spatial kNN would manufacture precisely the halo
## artefact step 18 measured -- it would hand every cell its neighbours' identity and then
## report, as a finding, that neighbours share identity. The projection uses transcriptome
## coordinates and nothing else.
##
## THE PCA OBJECT HAS TO BE REBUILT. 19a saved only `po$reduction.data` into the Seurat
## object, but projection needs the full pcaobj -- grate, phi, the pearson-residual mean and
## sd, scale.max. The PCA is recomputed here from the same HVGs and parameters and then
## VERIFIED against the stored embeddings: if the correlation is not ~1 the recomputed space
## is not the space the clusters live in and the projection is meaningless.
##
## THE GATE, PRE-COMMITTED. HCA tonsil label transfer was disabled in this project because
## its labels anti-tracked marker evidence, and this step must not repeat that.
##   Gate A  hold 10% of high-confidence cells out of the kNN REFERENCE, transfer to them,
##           and require agreement with their own cluster label. Tests the vote.
##   Gate B  among RESCUED cells where argmax resolves, require agreement between the
##           transferred label and argmax to be comparable to the 66.2% the two callers
##           reach inside the high-confidence set. Tests projection and vote together
##           against independent evidence.
## If Gate B falls below MIN_GATE_B, `celltype_projected` ships as NA rather than as a
## label. That rule is fixed here, before the numbers are seen.
##
## GATE B IS APPLIED PER STRATUM, NOT POOLED, and the first run showed why. Pooled it read
## 0.417 and rejected everything. Split by over-fixation status it is two different
## populations:
##     non-defect   0.677   -- ABOVE the 0.662 the two callers reach inside the HC object
##     defect       0.231   -- flat at 0.26-0.29 across a 7x depth span
## Defect cells are not shallow (median nCount 623, deeper than the non-defect cells that
## succeed at 475); their reads are simply wrong, which is also why imputation is not the
## remedy. Rejecting 344k good labels because 1.6M bad ones were averaged in with them is a
## failure of the gate's resolution, not a finding. The threshold itself is unchanged.
##
## Usage: sbatch --mem=300G --cpus-per-task=4 --time=12:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_rescue.R
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(scPearsonPCA)
  library(BiocNeighbors); setDTthreads(2)
})
OUT <- "output/19_celltypes"
IMPORT <- "output/01_import_v1"
## THE NEIGHBOUR SEARCH IS THE WHOLE COST, and RANN was the wrong tool for it. A kd-tree
## degenerates towards brute force as dimension rises, and at 30 PCs with a 400,000-cell
## reference and 2.4M queries RANN::nn2 ran for six hours without finishing -- twice. The
## first failure was misattributed to the majority vote; vectorising that helped but did not
## touch the real bottleneck. HNSW is an approximate graph index built for exactly this
## regime and returns the same answer to within a hair for a vote over 15 neighbours.
K_NN      <- 15L      # neighbours voting
N_REF     <- 400000L  # reference cells sampled from the HC object for the kNN
N_PCS     <- 20L      # PCs used for the search; beyond ~20 the extra dims are noise here
MIN_VOTE  <- 0.50     # a vote below this leaves the cell unresolved
MIN_GATE_A<- 0.85
MIN_GATE_B<- 0.50     # HC baseline argmax-vs-cluster agreement is 0.662
SEED <- 1L; set.seed(SEED)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

lab <- fread(file.path(OUT,"hc_cells_labelled_high.csv.gz"), na.strings=c("NA",""))
md  <- fread("output/16_regions/cells_meta_all.csv.gz",
             select=c("cell_id","samp","is_defect","tier_matched","nCount_RNA"))
ct  <- fread("output/16_celltypes/cell_celltypes.csv.gz",
             select=c("cell_id","celltype","lineage16"))
cat(sprintf("HC labelled %d | all cells %d | rescue target %d\n",
            nrow(lab), nrow(md), nrow(md) - nrow(lab)))

## ---- 1. rebuild and VERIFY the high-confidence PCA ---------------------------
seu <- readRDS(file.path(OUT,"hc_02_clustered.Rds"))
ca  <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
hv  <- VariableFeatures(seu)
cat(sprintf("\nHC object: %d cells, %d HVGs\n", ncol(seu), length(hv)))
tc  <- Matrix::colSums(ca)
gf  <- scPearsonPCA::gene_frequency(ca)
ch  <- as(ca[hv,],"dgCMatrix")
po  <- sparse_quasipoisson_pca_seurat(ch, totalcounts=tc, grate=gf[hv], scale.max=10,
                                      do.scale=TRUE, do.center=TRUE)
emb_new <- Embeddings(po$reduction.data)
emb_old <- Embeddings(seu, "pearsonpca")
cmn <- intersect(rownames(emb_new), rownames(emb_old))
cors <- sapply(1:10, function(i) abs(cor(emb_new[cmn,i], emb_old[cmn,i])))
cat("PCA rebuild check -- |cor| with the stored embedding, PC1-10:\n"); print(round(cors,4))
if (min(cors) < 0.99)
  stop("the rebuilt PCA is not the space the clusters live in (min |cor| ",
       round(min(cors),3), "); projection would be meaningless", call. = FALSE)
rm(ch); invisible(gc())

REF <- data.table(cell_id = rownames(emb_new))
REF <- lab[, .(cell_id, celltype_cluster, hc_cl = cl)][REF, on=.(cell_id)]
rm(seu); invisible(gc())

## ---- 2. project the rescue cells --------------------------------------------
rescue_ids <- setdiff(md$cell_id, lab$cell_id)
cat(sprintf("\nrescue cells: %d\n", length(rescue_ids)))
## The projection is checkpointed: it costs ~10 minutes of BPCells reads and must not be
## repeated when the vote downstream is changed.
PROJ_F <- file.path(OUT, "rescue_projection.Rds")
proj <- if (file.exists(PROJ_F)) { cat("\nreusing ", PROJ_F, "\n", sep=""); readRDS(PROJ_F) } else
rbindlist(lapply(c("IB1","RB1"), function(S) {
  f <- file.path(IMPORT, sprintf("seuratObject_%s_BPCells_full.Rds", S))
  cat("\nreading ", f, "\n", sep="")
  o <- readRDS(f)
  keep <- intersect(colnames(o), rescue_ids)
  cat(sprintf("  %s: %d rescue cells\n", S, length(keep)))
  o <- subset(o, cells = keep)
  cq <- if (!is.null(o[["RNA"]]$counts)) o[["RNA"]]$counts else o[["RNA"]]@counts
  tq <- Matrix::colSums(cq)
  xq <- as(cq[names(po$grate), ],"dgCMatrix")
  rm(o, cq); invisible(gc())
  pr <- project_sparse_quasipoisson_pca_seurat(xq, pcaobj = po, totalcounts = tq)
  e  <- Embeddings(pr$reduction.data)
  rm(xq, pr); invisible(gc())
  data.table(cell_id = rownames(e), e)
}), fill = TRUE)
if (!file.exists(PROJ_F)) saveRDS(proj, PROJ_F)
cat(sprintf("\nprojected %d cells into the HC PCA\n", nrow(proj)))
PCS <- grep("^PC_", names(proj), value=TRUE)[1:N_PCS]

## ---- 3. kNN vote -------------------------------------------------------------
ref_ok <- REF[!is.na(celltype_cluster)]
set.seed(SEED)
hold <- sample(nrow(ref_ok), max(1L, round(0.10*nrow(ref_ok))))
ref_use <- ref_ok[-hold]
ref_use <- ref_use[sample(.N, min(.N, N_REF))]
cat(sprintf("\nkNN reference: %d cells (10%% held out for gate A)\n", nrow(ref_use)))
RE <- emb_new[ref_use$cell_id, PCS, drop=FALSE]
## Index built once and reused for both the gate and the 2.4M rescue queries.
IDX <- BiocNeighbors::buildIndex(RE, BNPARAM = BiocNeighbors::HnswParam())
cat(sprintf("HNSW index built on %d x %d\n", nrow(RE), ncol(RE)))

## VECTORISED MAJORITY VOTE, IN CHUNKS. The first version built a 2.4M x 25 CHARACTER
## matrix and ran apply() over its rows; that spent five hours and was killed for memory.
## Labels are integer-coded, votes accumulated into an n x L integer table with one pass per
## neighbour (25 passes, not 2.4M row calls), and max.col picks the winner. Chunking bounds
## the nn2 allocation.
LEV <- sort(unique(ref_use$celltype_cluster))
refcode <- match(ref_use$celltype_cluster, LEV)
vote <- function(Q, chunk = 250000L) {
  n <- nrow(Q); out_lab <- character(n); out_frac <- numeric(n)
  for (st in seq(1L, n, by = chunk)) {
    en <- min(st + chunk - 1L, n); ii <- st:en
    r  <- BiocNeighbors::queryKNN(query = Q[ii, , drop = FALSE], k = K_NN,
                                  BNINDEX = IDX, get.distance = FALSE)
    cm <- matrix(refcode[r$index], nrow = length(ii))
    cnt <- matrix(0L, length(ii), length(LEV))
    rows <- seq_len(length(ii))
    for (j in seq_len(K_NN)) {
      ij <- cbind(rows, cm[, j]); cnt[ij] <- cnt[ij] + 1L
    }
    b <- max.col(cnt, ties.method = "first")
    out_lab[ii]  <- LEV[b]
    out_frac[ii] <- cnt[cbind(rows, b)] / K_NN
    rm(r, cm, cnt); if (st %% (chunk*4L) == 1L) invisible(gc())
    cat(sprintf("  voted %d / %d\n", en, n))
  }
  cbind(out_lab, out_frac)
}
## Gate A -- held-out HC cells
QA <- emb_new[ref_ok$cell_id[hold], PCS, drop=FALSE]
vA <- vote(QA)
gateA <- mean(vA[,1] == ref_ok$celltype_cluster[hold])
cat(sprintf("\nGATE A -- held-out HC cells, transferred vs their own cluster label: %.3f\n", gateA))

## rescue cells
QB <- as.matrix(proj[, ..PCS])
rownames(QB) <- proj$cell_id
vB <- vote(QB)
res <- data.table(cell_id = proj$cell_id, celltype_projected = vB[,1],
                  vote_frac = as.numeric(vB[,2]))
res[vote_frac < MIN_VOTE, celltype_projected := NA_character_]
## The pre-gate call is KEPT even when the gate rejects it. The first run wiped it, which
## made it impossible to ask WHY the gate failed without re-running the whole projection.
## A failed gate means "do not ship this as a label", not "destroy the evidence".
res[, celltype_projected_raw := celltype_projected]
res <- ct[res, on=.(cell_id)]
res <- md[, .(cell_id, samp, is_defect, tier_matched, nCount_RNA)][res, on=.(cell_id)]

FAMS <- names(LINEAGE_SETS_COARSE <- NULL) # placeholder replaced below
source("scripts/lib/lineage.R"); FAMS <- names(LINEAGE_SETS_COARSE)
map <- c(B="B cell", T="T cell", NK="NK cell", myeloid="Myeloid", endo="Endothelial",
         fibro="Stromal/FDC", PC="Plasma cell")
gb <- res[lineage16 %in% FAMS & !is.na(celltype_projected)]
gateB <- gb[, mean(map[lineage16] == celltype_projected, na.rm=TRUE)]
cat(sprintf("GATE B -- rescued cells where argmax resolves (%s): agreement %.3f (HC baseline 0.662)\n",
            format(nrow(gb), big.mark=","), gateB))

## ---- WHY did it fail? sparsity or corruption? ------------------------------
## Imputation before projection is only defensible if the problem is SPARSITY -- real signal
## diluted by shot noise. It is harmful if the problem is CORRUPTION: over-fixation degraded
## the RNA before the reads existed, and imputing from a biased observation fabricates
## plausible values. 68% of the rescue set is the over-fixation defect, so this stratification
## decides whether imputation is even the right kind of tool.
cat("\n=== GATE B stratified: is the failure defect, or depth? ===\n")
print(gb[, .(cells=.N, agreement=round(mean(map[lineage16]==celltype_projected),3),
             med_nCount=as.numeric(median(nCount_RNA)),
             med_vote=round(median(vote_frac),2)),
         by=.(is_defect, tier_matched)][order(-cells)])
cat("\n-- within NON-DEFECT rescue cells, by depth decile --\n")
nd <- gb[is_defect == FALSE]
if (nrow(nd) > 1000) {
  nd[, dec := cut(nCount_RNA, quantile(nCount_RNA, 0:10/10), include.lowest=TRUE, labels=FALSE)]
  print(nd[, .(cells=.N, med_nCount=as.numeric(median(nCount_RNA)),
               agreement=round(mean(map[lineage16]==celltype_projected),3)), by=dec][order(dec)])
}
cat("\n-- within DEFECT rescue cells, by depth decile (if depth were the issue these rise too) --\n")
df <- gb[is_defect == TRUE]
if (nrow(df) > 1000) {
  df[, dec := cut(nCount_RNA, quantile(nCount_RNA, 0:10/10), include.lowest=TRUE, labels=FALSE)]
  print(df[, .(cells=.N, med_nCount=as.numeric(median(nCount_RNA)),
               agreement=round(mean(map[lineage16]==celltype_projected),3)), by=dec][order(dec)])
}
gateB_nd <- nd[, mean(map[lineage16]==celltype_projected)]
cat(sprintf("\nGate B on NON-DEFECT rescue cells alone: %.3f (threshold %.2f, HC baseline 0.662)\n",
            gateB_nd, MIN_GATE_B))
cat(if (gateB_nd >= MIN_GATE_B)
  "  -> the non-defect subset PASSES; the failure is the defect cells, and a restricted rescue is viable.\n" else
  "  -> the non-defect subset also fails; the failure is not merely the defect.\n")

## Gate A is global: if the vote does not work inside the space it was built from, nothing
## downstream is interpretable and no stratum can rescue it.
if (gateA < MIN_GATE_A) {
  cat("\n*** GATE A FAILED -- the vote is unreliable in its own space; nothing ships ***\n")
  res[, celltype_projected := NA_character_]
} else {
  res[, stratum := fifelse(is_defect, "defect", "non-defect")]
  gb[, stratum := fifelse(is_defect, "defect", "non-defect")]
  gbs <- gb[, .(cells = .N,
                agreement = mean(map[lineage16] == celltype_projected, na.rm = TRUE)), by = stratum]
  gbs[, passes := agreement >= MIN_GATE_B]
  cat("\n=== GATE B, per stratum (threshold ", MIN_GATE_B, ", HC baseline 0.662) ===\n", sep="")
  print(gbs[, .(stratum, cells, agreement = round(agreement,3), passes)])
  fail <- gbs[passes == FALSE, stratum]
  if (length(fail)) {
    cat("\nwithholding labels for: ", paste(fail, collapse=", "), "\n", sep="")
    res[stratum %in% fail, celltype_projected := NA_character_]
  }
  if (nrow(gbs[passes == TRUE])) cat("shipping labels for: ",
      paste(gbs[passes == TRUE, stratum], collapse=", "), "\n", sep="")
}

res[, confidence := "rescued"]
res[is.na(celltype_projected), confidence := "unresolved"]
cat("\n=== what shipped, by stratum ===\n")
print(res[, .(cells=.N, labelled=sum(!is.na(celltype_projected)),
              pct=round(100*mean(!is.na(celltype_projected)),1)),
          by=.(is_defect, tier_matched)][order(-cells)])
cat("\n=== rescued label totals ===\n")
print(res[, .N, by=.(celltype_projected)][order(-N)])
cat("\n=== by tier ===\n")
print(res[, .(cells=.N, resolved=sum(!is.na(celltype_projected))),
          by=.(samp, is_defect, tier_matched)][order(samp,-cells)])
sc("rescued_cells.csv.gz",
   res[, .(cell_id, samp, celltype_projected, celltype_projected_raw, vote_frac, confidence,
           celltype_argmax = celltype, lineage_argmax = lineage16, is_defect, tier_matched)])

## ---- 4. the single all-cell table --------------------------------------------
allc <- rbind(
  lab[, .(cell_id, samp = merge_sample, confidence = "high",
          celltype_cluster, cluster_status, is_niche, halo_suspect,
          celltype_projected = NA_character_, vote_frac = NA_real_)],
  res[, .(cell_id, samp, confidence, celltype_cluster = NA_character_,
          cluster_status = NA_character_, is_niche = NA, halo_suspect = NA,
          celltype_projected, vote_frac)], fill = TRUE)
allc <- ct[allc, on=.(cell_id)]
setnames(allc, c("celltype","lineage16"), c("celltype_argmax","lineage_argmax"))
cat(sprintf("\nall-cell table: %d rows (expected %d)\n", nrow(allc), nrow(md)))
stopifnot(nrow(allc) == nrow(md), !anyDuplicated(allc$cell_id))
cat("\n=== confidence classes ===\n"); print(allc[, .N, by=confidence][order(-N)])
sc("all_cells_celltypes.csv.gz", allc)
cat("\ndone\n")
