## ============================================================================
## 16_script_celltypes.R -- cross-sample-comparable cell types, both samples at once.
##
## WHY NOT THE EXISTING LABELS. Neither existing cell-type call can carry a cross-sample
## comparison:
##   * `pearson_clusters_low/high` ids are PER SAMPLE. IB1 cluster 3 and RB1 cluster 3 are
##     unrelated, and the step-09b marker audit that turns ids into labels has only ever
##     been run for IB1 (output/09b_marker_relabel/ holds IB1 only).
##   * `lineage` (from assign_lineage in step 01) IS comparable -- it is a marker-ratio call
##     with no clustering in it -- and it reproduces the report's headline composition
##     numbers exactly (8.5% B in IB1, 39.9% in RB1). But it has no T or myeloid SUBSETS,
##     and the step-16 checkpoint question is specifically about Tfh (the PD-1 source) and
##     macrophages/cDC (the PD-L1 source).
## So this re-runs the same caller, on both samples in ONE call so the fold-over-mean gene
## scaling is shared, with a second hierarchical stage inside T and inside myeloid.
##
## TWO STAGES, NOT ONE FLAT ARGMAX. Putting Tfh/Treg/CD8/CD4 in as siblings of B and myeloid
## would force every T cell into a subset on markers far sparser than CD3D/E/G/TRAC. That is
## exactly the failure step 09b documented -- "resolved to parent lineage; subset call not
## supported" for two of its T clusters. Instead stage 1 is step 01's coarse call verbatim
## (which also makes it a regression test), and stage 2 runs only within the parent lineage,
## where cells that cannot be resolved stay `T_unresolved` / `myeloid_unresolved` rather than
## being assigned by noise.
##
## THE SUBTYPE SETS DELIBERATELY EXCLUDE THE CHECKPOINT GENES. Tfh is normally defined with
## PDCD1 and Treg with CTLA4, but PDCD1 and CTLA4 are the READOUT in step 16e. Defining the
## cell type with them would make "Tfh express PD-1" a restatement of the definition. So
## Tfh here is CXCL13/ICOS/CXCR5/BTLA/CD40LG and Treg is FOXP3/IL2RA/IKZF2, and every
## checkpoint gene is held out. This is the same discipline step 14 applies to IGKC/IGLC.
##
## TIER-INDEPENDENT BY DESIGN. Called once over all 4.4M pre-QC cells, exactly as step 01
## does, so the per-gene means come from whole tissue rather than from whichever QC tier
## happens to be in play. Each tier then subsets this table.
##
## Usage: sbatch --mem=120G --cpus-per-task=2 --time=2:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/16_script_celltypes.R
## ============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); setDTthreads(1)
})

if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/lineage.R")

META_F <- "output/16_regions/cells_meta_all.csv.gz"
CNT_F  <- "output/16_regions/panel_counts_all.Rds"
NEG_F  <- "output/16_regions/negprobe_counts_all.Rds"
OUT    <- "output/16_celltypes"
for (f in c(META_F, CNT_F, NEG_F))
  if (!file.exists(f)) stop("missing ", f, "\n  run scripts/16_script_cell_meta.R and ",
                            "scripts/16_script_panel_counts.R first", call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

## ---- stage-2 marker sets -----------------------------------------------------
## Every checkpoint gene is HELD OUT of these on purpose (see the header).
T_SETS <- list(
  Tfh      = c("CXCL13","ICOS","CXCR5","BTLA","CD40LG"),  # NOT PDCD1 -- that is the readout
  Treg     = c("FOXP3","IL2RA","IKZF2"),                  # NOT CTLA4 -- ditto
  CD8_T    = c("CD8A","GZMK","GZMB"),
  CD4_conv = c("CD4","IL7R","CCR7"))
M_SETS <- list(
  Macrophage = c("CD68","CD163","C1QA","C1QB","MRC1","CSF1R"),
  cDC        = c("BATF3","CLEC9A","CD1C","FCER1A","LAMP3"),
  Monocyte   = c("CD14","LYZ","ITGAX"))
## The germinal-centre programme is reported as a SCORE on B cells rather than a hard call.
## FL is GC-derived, so this is the axis that matters for the FL comparison, but a 4-gene
## argmax against the B-defining markers would be circular and would invent a precision the
## panel does not support (AICDA is detected in ~1% of cells).
GC_SET  <- c("BCL6","AICDA","RGS13","LMO2")
SEED    <- 1L

cat("============================================================\n")
cat("step 16b -- cross-sample cell types, hierarchical marker call\n")
cat("============================================================\n")

md <- fread(META_F, select = c("cell_id","samp","lineage","nCount_RNA",
                               "tier_strict","tier_matched","is_defect"))
M  <- readRDS(CNT_F); NP <- readRDS(NEG_F)
if (!identical(colnames(M), colnames(NP))) stop("panel and negprobe columns disagree", call. = FALSE)
gi <- match(md$cell_id, colnames(M))
if (anyNA(gi)) stop(sum(is.na(gi)), " metadata cells absent from the count matrix", call. = FALSE)
M <- M[, gi, drop = FALSE]; NP <- NP[, gi, drop = FALSE]
NEG_TOT <- Matrix::colSums(NP); N_NEG <- nrow(NP)
cat(sprintf("cells %d | panel genes %d | negprobes %d | median background %.3f\n",
            nrow(md), nrow(M), N_NEG, median(NEG_TOT / N_NEG)))
rm(NP); invisible(gc())

## ---- stage 1: the coarse call, step 01's sets and parameters -----------------
## RUN TWICE, ON PURPOSE, because the two runs answer different questions.
##
## Step 01 calls assign_lineage INSIDE its per-sample loop (01_import_cosmx_seurat.R:190
## onwards), so the stored `lineage` column was scaled against per-sample gene means. Step
## 16 needs the opposite -- pooled means, so that "B cell" means the same thing in IB1 and
## RB1 and a composition difference between the samples is not partly an artefact of each
## sample having been standardised to itself.
##
##   lineage_ps  per-sample scaling  -> must reproduce the stored column (REGRESSION TEST)
##   lineage16   pooled scaling      -> what everything downstream uses
##
## A first version asserted that the pooled call reproduce the stored column and tripped at
## 95.99%. That assertion was simply wrong: it demanded reproduction of a per-sample call
## from a deliberately pooled one. The substantive check is the one below -- that pooling
## only moves cells across the AMBIGUOUS margin and never flips one lineage into another.
set.seed(SEED)
ps <- rbindlist(lapply(unique(md$samp), function(S) {
  ii <- md[, which(samp == S)]
  r <- do.call(assign_lineage, c(list(cts = M[, ii, drop = FALSE], sets = LINEAGE_SETS_COARSE,
                                      neg_total = NEG_TOT[ii], n_neg = N_NEG),
                                 LINEAGE_CALL_ARGS))
  data.table(row = ii, lineage_ps = r$label)
}))
md[ps$row, lineage_ps := ps$lineage_ps]
agree_ps <- md[, mean(lineage_ps == lineage)]
cat(sprintf("\n=== regression test: per-sample call vs the stored `lineage` column: %.4f%%\n",
            100*agree_ps))
sc("stage1_persample_vs_stored.csv", md[, .N, by = .(lineage, lineage_ps)][order(-N)])
if (agree_ps < 0.99) {
  print(dcast(md[, .N, by = .(lineage, lineage_ps)], lineage ~ lineage_ps,
              value.var = "N", fill = 0L))
  stop("the per-sample call does not reproduce the stored `lineage` column (",
       round(100*agree_ps,2), "%) -- the panel extraction or cell ordering is wrong",
       call. = FALSE)
}
cat("per-sample call reproduces the stored column -- extraction and ordering are correct\n")

set.seed(SEED)
s1 <- do.call(assign_lineage, c(list(cts = M, sets = LINEAGE_SETS_COARSE,
                                     neg_total = NEG_TOT, n_neg = N_NEG),
                                LINEAGE_CALL_ARGS))
md[, lineage16 := s1$label]
md[, l2m16 := s1$log2_margin]

## What did pooling change? Split the disagreement into (a) margin moves, where one call is
## a sentinel (ambiguous / below_background / all_markers_0) and the other is a lineage, and
## (b) HARD FLIPS, where both calls name a lineage and they differ. Only (b) would threaten
## the cross-sample comparison, because only (b) means a cell changed identity.
SENT <- c("ambiguous","below_background","all_markers_0")
md[, `:=`(ps_sent = lineage_ps %in% SENT, p16_sent = lineage16 %in% SENT)]
n_dis  <- md[, sum(lineage_ps != lineage16)]
n_hard <- md[, sum(lineage_ps != lineage16 & !ps_sent & !p16_sent)]
cat(sprintf("\n=== pooled vs per-sample scaling\n  agreement %.2f%% | disagreements %d (%.2f%%)\n",
            100*md[, mean(lineage_ps == lineage16)], n_dis, 100*n_dis/nrow(md)))
cat(sprintf("  of those, margin moves (one side a sentinel): %d (%.1f%%)\n",
            n_dis - n_hard, 100*(n_dis-n_hard)/max(n_dis,1)))
cat(sprintf("  HARD FLIPS (lineage -> different lineage):    %d (%.3f%% of all cells)\n",
            n_hard, 100*n_hard/nrow(md)))
sc("pooled_vs_persample.csv", md[, .N, by = .(lineage_ps, lineage16)][order(-N)])
if (n_hard / nrow(md) > 0.005) {
  print(md[lineage_ps != lineage16 & !ps_sent & !p16_sent, .N,
           by = .(lineage_ps, lineage16)][order(-N)])
  stop("pooling flipped ", n_hard, " cells between named lineages (>0.5%) -- the two ",
       "samples are too different for a shared scaling; reconsider", call. = FALSE)
}
cat("  -> pooling only moves cells across the ambiguity margin. Safe to use pooled.\n")
md[, c("ps_sent","p16_sent") := NULL]
print(md[, .(cells = .N), by = .(samp, lineage16)][order(samp, -cells)])

## ---- stage 2: within T, and within myeloid ----------------------------------
## Each stage-2 call sees only its parent's cells, so the fold-over-mean scaling is computed
## among T cells (or myeloid cells) rather than against whole tissue -- which is the point:
## CXCL13 in a T cell should be judged against other T cells.
sub_call <- function(parent, sets, prefix) {
  ii <- md[, which(lineage16 == parent)]
  cat(sprintf("\n=== stage 2: %s (%d cells) -> %s\n", parent, length(ii),
              paste(names(sets), collapse = ", ")))
  if (length(ii) < 100) { cat("  too few cells; skipped\n"); return(NULL) }
  set.seed(SEED)
  r <- do.call(assign_lineage, c(list(cts = M[, ii, drop = FALSE], sets = sets,
                                      neg_total = NEG_TOT[ii], n_neg = N_NEG),
                                 LINEAGE_CALL_ARGS))
  ## The sentinels mean the same thing they do at stage 1, but "below_background" inside a
  ## confirmed T cell is not an absent cell -- it is a T cell whose subset markers are all
  ## at background. Renamed so that distinction survives into the tables.
  lab <- r$label
  lab[lab %in% c("below_background","all_markers_0","ambiguous")] <- paste0(prefix, "_unresolved")
  data.table(row = ii, subtype = lab, sub_l2m = r$log2_margin)
}
s2 <- rbindlist(list(sub_call("T", T_SETS, "T"), sub_call("myeloid", M_SETS, "myeloid")))

md[, celltype := lineage16]
if (nrow(s2)) md[s2$row, `:=`(celltype = s2$subtype, sub_l2m = s2$sub_l2m)]

cat("\n=== final celltype, both samples\n")
ct <- md[, .N, by = .(samp, celltype)]
ct[, pct := round(100*N/sum(N), 2), by = samp]
print(dcast(ct, celltype ~ samp, value.var = c("N","pct"))[order(-N_IB1)])
sc("celltype_counts_preqc.csv", ct)

## ---- validation: do the called subtypes carry the right markers? ------------
## Self-validating rather than trusting the call: mean counts per 10k per subtype for the
## DEFINING markers and, more informatively, for markers the call did NOT use. A Tfh call
## that does not elevate held-out PDCD1 is not a Tfh call.
HOLDOUT <- c("PDCD1","CTLA4","CD274","PDCD1LG2","TNFRSF9","LAG3","TIGIT","HAVCR2","TOX",
             "CD3E","MS4A1","CD68","CD163","BATF3","LAMP3","MKI67","BCL2")
SHOW <- unique(c(unlist(T_SETS, use.names=FALSE), unlist(M_SETS, use.names=FALSE),
                 GC_SET, HOLDOUT))
SHOW <- intersect(SHOW, rownames(M))
dep <- md$nCount_RNA
prof <- rbindlist(lapply(unique(md$celltype), function(k) {
  ii <- which(md$celltype == k)
  if (length(ii) < 50) return(NULL)
  tot <- sum(dep[ii])
  data.table(celltype = k, cells = length(ii), gene = SHOW,
             cp10k = round(1e4 * as.numeric(Matrix::rowSums(M[SHOW, ii, drop=FALSE])) / tot, 2),
             pct_det = round(100 * as.numeric(Matrix::rowMeans(M[SHOW, ii, drop=FALSE] > 0)), 2))
}))
prof[, is_defining := gene %in% c(unlist(T_SETS, use.names=FALSE),
                                  unlist(M_SETS, use.names=FALSE))]
sc("celltype_marker_profile.csv", prof)
cat("\n=== held-out marker check: cp10k by celltype (the call never saw these)\n")
print(dcast(prof[gene %in% c("PDCD1","CTLA4","CD274","PDCD1LG2","CD3E","CD68","BATF3")],
            celltype ~ gene, value.var = "cp10k"))

## Rank check. The `genuinely_held_out` column is the point: a first version of this table
## validated Macrophage on CD68 and cDC on BATF3, which are their own DEFINING markers -- a
## circular check that could not fail. It is computed here rather than asserted, so the
## distinction cannot be lost again.
##
## On this panel only Tfh and Treg have real held-out validators, because every myeloid and
## DC marker present is definitional. That is a limitation of the 6K panel, not something to
## paper over: Macrophage/cDC rows below are DESCRIPTIVE separation, not independent
## validation. What does support them is the coarse `myeloid` call itself, which reproduces
## the stored `lineage` column at 100% above.
## RANKED AMONG SIBLINGS, not among all 17 cell types. Ranking a T subset against
## endothelium is trivially winnable and says nothing; the question a subtype label has to
## answer is whether it separates from the OTHER cells of its own parent lineage. A first
## version ranked across all 17 and also validated Tfh on TNFRSF9, which is an activation
## marker rather than a Tfh marker -- a validator no one claims, so its failure meant
## nothing. Both are corrected here.
DEFINING <- c(T_SETS, M_SETS)
SIBS <- list(T = c(names(T_SETS), "T_unresolved"),
             myeloid = c(names(M_SETS), "myeloid_unresolved"))
PARENT <- c(setNames(rep("T", length(T_SETS)), names(T_SETS)),
            setNames(rep("myeloid", length(M_SETS)), names(M_SETS)))
chk <- rbindlist(lapply(list(
    c(sub="Tfh",        gene="PDCD1"),
    c(sub="Tfh",        gene="CXCL13"),
    c(sub="Treg",       gene="CTLA4"),
    c(sub="CD8_T",      gene="GZMB"),
    c(sub="Macrophage", gene="CD68"),
    c(sub="cDC",        gene="BATF3")), function(z) {
  sib <- SIBS[[PARENT[[z[["sub"]]]]]]
  pp <- prof[gene == z[["gene"]] & celltype %in% sib][order(-cp10k)]
  if (!nrow(pp) || !z[["sub"]] %in% pp$celltype) return(NULL)
  data.table(subtype = z[["sub"]], gene = z[["gene"]],
             genuinely_held_out = !z[["gene"]] %in% DEFINING[[z[["sub"]]]],
             rank_among_siblings = which(pp$celltype == z[["sub"]]), of = nrow(pp),
             cp10k = pp[celltype == z[["sub"]], cp10k], top = pp$celltype[1],
             top_cp10k = pp$cp10k[1],
             ## how much better than the best sibling that is not this subtype
             margin_vs_next = round(pp[celltype == z[["sub"]], cp10k] /
                                    max(pp[celltype != z[["sub"]], cp10k]), 3))
}))
chk[, pass := rank_among_siblings == 1L]
cat("\n=== subtype rank on marker genes AMONG SIBLINGS (only held-out rows are evidence)\n")
print(chk)
sc("celltype_holdout_check.csv", chk)
real <- chk[genuinely_held_out == TRUE]
cat(sprintf("\ngenuinely held-out checks: %d, of which passing (rank <= 2): %d\n",
            nrow(real), real[, sum(pass)]))
if (nrow(real) && !all(real$pass))
  cat("\n*** WARNING: ", paste(real[pass == FALSE, paste0(subtype,"/",gene)], collapse=", "),
      " is not the top sibling on its held-out marker.\n*** That subtype does not separate",
      " from its siblings on independent evidence; label it as weakly\n*** supported wherever",
      " it appears downstream.\n", sep="")
if (nrow(chk[genuinely_held_out == FALSE]))
  cat("NOTE: ", paste(chk[genuinely_held_out == FALSE, paste0(subtype,"/",gene)], collapse=", "),
      " are DEFINING markers -- descriptive separation only, not validation.\n", sep="")

## ---- GC programme score on B cells ------------------------------------------
gs <- intersect(GC_SET, rownames(M))
md[, gc_ct := as.numeric(Matrix::colSums(M[gs, , drop = FALSE]))]
md[, gc_cp10k := 1e4 * gc_ct / nCount_RNA]
cat(sprintf("\n=== GC programme (%s) in B cells, pre-QC\n", paste(gs, collapse="+")))
print(md[lineage16 == "B", .(cells = .N, pct_gc_pos = round(100*mean(gc_ct > 0),2),
         med_gc_cp10k = round(median(gc_cp10k),2)), by = samp])

## ---- write -------------------------------------------------------------------
setcolorder(md, c("cell_id","samp","celltype","lineage16","lineage","gc_ct","gc_cp10k"))
fwrite(md[, .(cell_id, samp, celltype, lineage16, lineage_ps, lineage, l2m16, sub_l2m,
              gc_ct, gc_cp10k)],
       file.path(OUT, "cell_celltypes.csv.gz"), compress = "gzip")
cat("\nwrote ", file.path(OUT, "cell_celltypes.csv.gz"), " (", nrow(md), " rows)\n", sep="")
writeLines(format(Sys.time()), file.path(OUT, "celltypes.done"))
cat("============================================================\nOK\n")
