## ============================================================================
## 16_script_cell_meta.R -- one metadata table for ALL cells in both samples.
##
## Step 16 compares FL between IB1 and RB1 at two QC tiers (see the plan and CLAUDE.md):
##   strict   qc_passing_all           -- IB1 400,111 / RB1 1,399,189 (asymmetric: the
##                                       staining rule was never evaluated on RB1)
##   matched  fovpx ok & qcCellsPassed -- IB1 1,877,697 / RB1 1,399,189 (symmetric rules)
##
## The strict tier is already materialised in output/10_hrs/hrs_segments.csv.gz, but the
## matched tier is not, and `lineage` / `inSitu_clusters` / `x_slide_mm` exist only inside
## the step-01 objects -- the vendor flat files carry neither `lineage` (project-added by
## scripts/01_import_cosmx_seurat.R) nor slide-mm coordinates.
##
## So this runs ONCE, over the PRE-QC objects (every cell, both samples), and every later
## 16* script is pure CSV work at either tier. That is the only reason this script exists:
## to keep one heavy object read out of the per-tier loop.
##
## NOTE the pre-QC objects are the right source precisely BECAUSE they are pre-QC -- a
## tier is then a filter on this table rather than a different object, so the two tiers
## are guaranteed to be nested and comparable.
##
## Usage: sbatch --mem=120G --cpus-per-task=2 --time=2:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/16_script_cell_meta.R
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(data.table); setDTthreads(1)
})

SAMPLES <- c("IB1", "RB1")
PRE_RDS <- "output/01_import_v1/seuratObject_%s_BPCells_full.Rds"
QC_F    <- "output/02_spatial_QC/qc_spatial_cells.csv.gz"
OUT     <- "output/16_regions"
OUT_F   <- file.path(OUT, "cells_meta_all.csv.gz")

if (!file.exists("scripts/lib/cluster_labels.R")) stop("run from the project root", call. = FALSE)
if (!file.exists(QC_F)) stop("missing ", QC_F, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("step 16a-meta -- per-cell metadata, all cells, both samples\n")
cat("============================================================\n")

## ---- 1. metadata out of the pre-QC objects ----------------------------------
## Same `need` set and the same failure mode as scripts/11_script_hodgkin_niche.R:127 --
## if a column is missing, say which, rather than producing a table with a silent NA
## column that only breaks three scripts later.
NEED <- c("Run_Tissue_name", "x_slide_mm", "y_slide_mm", "nCount_RNA",
          "inSitu_clusters", "lineage", "Area.um2", "nCount_negprobes",
          "nFeature_negprobes", "qcCellsPassed")

read_one <- function(S) {
  f <- sprintf(PRE_RDS, S)
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
  cat("reading ", f, "\n", sep = "")
  o  <- readRDS(f)
  rn <- rownames(o@meta.data)
  miss <- setdiff(NEED, colnames(o@meta.data))
  if (length(miss)) stop(f, " meta.data is missing: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  m <- as.data.table(o@meta.data[, NEED, drop = FALSE])
  m[, cell_id := rn]
  rm(o); invisible(gc())
  cat(sprintf("  %s: %d cells\n", S, nrow(m)))
  m[]
}

md <- rbindlist(lapply(SAMPLES, read_one))
setnames(md, "Run_Tissue_name", "samp")
md[, samp := as.character(samp)]
if (!all(md$samp %in% SAMPLES))
  stop("unexpected Run_Tissue_name values: ",
       paste(setdiff(unique(md$samp), SAMPLES), collapse = ", "), call. = FALSE)

## FOV from the cell id (2nd-to-last token of c_<slide>_<fov>_<cell>), as
## 11_script_hodgkin_niche.R:184 does -- the objects do not agree on whether a `fov`
## metadata column is present.
md[, fov_id := as.integer(vapply(strsplit(cell_id, "_", fixed = TRUE),
                                 function(z) z[length(z) - 1L], ""))]
if (anyNA(md$fov_id)) stop("could not parse fov from every cell_id", call. = FALSE)

md[, lineage := as.character(lineage)]
md[, inSitu_clusters := as.character(inSitu_clusters)]
md[, is_defect := inSitu_clusters == "i"]

cat(sprintf("\ntotal cells: %d\n", nrow(md)))
print(md[, .(cells = .N, defect_pct = round(100 * mean(is_defect), 1),
             med_depth = as.numeric(median(nCount_RNA)),
             med_area = round(median(Area.um2), 1)), by = samp])

## ---- 2. QC flags, and the two tiers -----------------------------------------
## The staining rule is a SEPARATE column from the pixel/FOV rule, which is what makes a
## symmetric tier possible at all: dropping `qc_spatial_staining` leaves rules both
## samples were actually evaluated on. Verified against qc_spatial_attrition.csv below.
## Join on `cell`, NOT `cell_id_prefixed`. The prefixed form is "IB1_c_2_1_1" while the
## object rownames are "c_2_1_1"; the bare ids are globally unique anyway because the two
## samples are on different slides (IB1 = slide 2, RB1 = slide 1), verified below.
qc <- fread(QC_F, select = c("cell", "sample_id", "qc_spatial_fovpx",
                             "qc_spatial_staining", "qc_spatial_i_fraction",
                             "qcCellsPassed", "qc_spatial_passing", "qc_passing_all"))
setnames(qc, c("cell", "sample_id", "qcCellsPassed"),
         c("cell_id", "samp_qc", "qcCellsPassed_qc"))
if (anyDuplicated(qc$cell_id))
  stop("cell ids are not unique across samples in the QC table", call. = FALSE)
cat(sprintf("\nQC table: %d rows\n", nrow(qc)))

n0 <- nrow(md)
md <- merge(md, qc, by = "cell_id", all.x = TRUE, sort = FALSE)
if (nrow(md) != n0) stop("QC merge changed the row count: ", n0, " -> ", nrow(md),
                         " (duplicate cell ids?)", call. = FALSE)
nmiss <- md[, sum(is.na(qc_passing_all))]
if (nmiss) stop(nmiss, " cells have no QC row -- the id conventions disagree", call. = FALSE)
if (md[, any(samp != samp_qc)]) stop("sample disagrees between object and QC table", call. = FALSE)
md[, samp_qc := NULL]

## Two tiers. `strict` is qc_passing_all verbatim; `matched` drops the staining term.
md[, tier_strict  := qc_passing_all == TRUE]
md[, tier_matched := qc_spatial_fovpx == "ok" & qcCellsPassed_qc == TRUE]

## Nesting is a property we rely on downstream (a strict result must be a subset of the
## matched tissue), so assert it rather than hope for it.
if (md[, any(tier_strict & !tier_matched)])
  stop("strict tier is not nested inside matched -- the tier definitions are wrong",
       call. = FALSE)

cat("\n=== QC tiers\n")
print(md[, .(cells = .N,
             strict = sum(tier_strict), matched = sum(tier_matched),
             strict_nodefect  = sum(tier_strict  & !is_defect),
             matched_nodefect = sum(tier_matched & !is_defect)), by = samp])

## Regression check against the numbers the plan was built on. These are hard-coded on
## purpose: if step 02 is ever re-run with different parameters, this is where it surfaces.
EXP <- data.table(samp = c("IB1","RB1"), strict = c(400111L, 1399189L),
                  matched = c(1877697L, 1399189L))
got <- md[, .(strict = sum(tier_strict), matched = sum(tier_matched)), by = samp]
chk <- merge(EXP, got, by = "samp", suffixes = c("_exp","_got"))
print(chk)
if (chk[, any(strict_exp != strict_got | matched_exp != matched_got)])
  stop("tier counts do not match the expected values -- step 02 output has changed",
       call. = FALSE)
cat("tier counts match the expected values\n")

## ---- 3. write ----------------------------------------------------------------
setcolorder(md, c("cell_id","samp","fov_id","x_slide_mm","y_slide_mm","nCount_RNA",
                  "Area.um2","lineage","inSitu_clusters","is_defect",
                  "tier_strict","tier_matched"))
fwrite(md, OUT_F, compress = "gzip")
cat("\nwrote ", OUT_F, " (", nrow(md), " rows, ", ncol(md), " cols)\n", sep = "")
writeLines(format(Sys.time()), file.path(OUT, "cells_meta_all.done"))
cat("============================================================\nOK\n")
