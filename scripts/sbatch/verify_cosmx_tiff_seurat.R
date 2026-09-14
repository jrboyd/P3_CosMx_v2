## ============================================================================
## Seurat-object half of the assertion suite for scripts/lib/cosmx_tiff.R.
##
## scripts/sbatch/verify_cosmx_tiff.R covers everything reachable without a
## readRDS. This covers what is left: cosmx_cells_to_px(), cosmx_fit_mm_affine()
## and the mm->px affine path, which need a real object. Needs a compute node
## because these objects are 0.2-4.2 GB on disk.
##
##   sbatch --mem=32G --time=1:00:00 --job-name=verify-cosmx-seu \
##     scripts/sbatch/run_rscript.sh \
##     scripts/sbatch/verify_cosmx_tiff_seurat.R \
##     output/03_filter_and_process/IB1_03_HCA_cell_typed.Rds IB1
##
## Usage: Rscript scripts/sbatch/verify_cosmx_tiff_seurat.R <seurat_rds> <sample_id>
## ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("Usage: Rscript scripts/sbatch/verify_cosmx_tiff_seurat.R <seurat_rds> <sample_id>",
       call. = FALSE)
RDS       <- args[1]
SAMPLE_ID <- args[2]

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root: cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2",
       call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})
setDTthreads(1)

V1_DIR <- "/gpfs1/pi/avolaric/files_jrboyd/P3_CosMx"
FLAT   <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"

n_pass <- 0L
check <- function(label, expr) {
  ok <- isTRUE(tryCatch(expr, error = function(e) FALSE))
  if (!ok) {
    res <- tryCatch(expr, error = function(e) conditionMessage(e))
    stop("FAIL ", label, if (is.character(res)) paste0(" -- ", res) else "", call. = FALSE)
  }
  n_pass <<- n_pass + 1L
  cat("  ok  ", label, "\n", sep = "")
}

cat("reading ", RDS, "\n", sep = "")
seu <- readRDS(RDS)
md  <- seu@meta.data
cat("cells: ", format(nrow(md), big.mark = ","), "\n", sep = "")

## ---------------------------------------------------------------------------
## S1 -- what the object actually carries. Reported, not just asserted: this is
## the definitive answer to which sample column exists, which drove the
## behaviour-changing fix in cosmx_cells_to_px().
## ---------------------------------------------------------------------------
cat("\n== S1 object contract ==\n")
for (cn in c("Run_Tissue_name", "sample_id", "merge_sample")) {
  if (cn %in% colnames(md)) {
    cat("  ", cn, ": ", paste(names(table(md[[cn]])), table(md[[cn]]),
                              sep = "=", collapse = "  "), "\n", sep = "")
  } else cat("  ", cn, ": ABSENT\n", sep = "")
}
cat("  images: ", paste(names(seu@images), collapse = ", "), "\n", sep = "")
cat("  resolved sample column: ",
    .cosmx_sample_col(NULL, colnames(md)), "\n", sep = "")

check("FOV-local px columns present",
      all(c("x_FOV_px", "y_FOV_px") %in% colnames(md)))
check("an image exists for this sample", SAMPLE_ID %in% names(seu@images))
check("centroid names match meta.data rownames",
      all(seu@images[[SAMPLE_ID]]$centroids@cells %in% rownames(md)))
check("cell ids parse to a plausible fov range", {
  f <- cosmx_fov_from_id(rownames(md))
  !anyNA(f) && min(f) >= 1L
})

## ---------------------------------------------------------------------------
## S2 -- cosmx_cells_to_px() against the flat-file exact route, cell for cell.
## Independent of the library: the flat file carries CenterX/Y_local_px directly.
## ---------------------------------------------------------------------------
cat("\n== S2 cells_to_px vs the flat file ==\n")
px <- cosmx_cells_to_px(SAMPLE_ID, seu)
cat("  mapped ", format(nrow(px), big.mark = ","), " cells\n", sep = "")
check("mapped at least one cell", nrow(px) > 0L)
check("no NA coordinates", !anyNA(px$X_px) && !anyNA(px$Y_px))

flat <- fread(sprintf(FLAT, SAMPLE_ID, SAMPLE_ID),
              select = c("cell", "fov", "CenterX_local_px", "CenterY_local_px"))
## Build the key with the same parser used on the Seurat side, rather than
## assuming what the flat file's cell_ID column counts.
flat$key <- cosmx_seu_key(flat$cell)
check("flat-file keys are unique", !anyDuplicated(flat$key))
check("flat-file fov agrees with the parsed id",
      identical(as.integer(flat$fov), cosmx_fov_from_id(flat$cell)))
ref <- cosmx_local_to_px(SAMPLE_ID, flat$fov,
                         flat$CenterX_local_px, flat$CenterY_local_px)
flat$X_ref <- ref$X_px; flat$Y_ref <- ref$Y_px

px$key <- cosmx_seu_key(px$cell_id)
j <- merge(px, flat[, c("key", "X_ref", "Y_ref")], by = "key")
cat("  joined ", format(nrow(j), big.mark = ","), " / ",
    format(nrow(px), big.mark = ","), " cells to the flat file\n", sep = "")
check("every mapped cell joins the flat file", nrow(j) == nrow(px))
check("X_px identical to the flat-file route", identical(j$X_px, j$X_ref))
check("Y_px identical to the flat-file route", identical(j$Y_px, j$Y_ref))

## ---------------------------------------------------------------------------
## S3 -- old vs new on the two Seurat-facing functions.
## v1's cells_to_tiff_px() only filtered by sample when Run_Tissue_name existed,
## so an exact match is expected only in that case; otherwise the new function
## filters where v1 did not, which is the point of the fix.
## ---------------------------------------------------------------------------
cat("\n== S3 old vs new ==\n")
v1 <- new.env(parent = globalenv())
sys.source(file.path(V1_DIR, "assemble_cosmx/16_cosmx_tiff_mapping.R"), v1)

o_px <- suppressWarnings(v1$cells_to_tiff_px(seu, SAMPLE_ID))
if ("Run_Tissue_name" %in% colnames(md)) {
  check("cells_to_px identical to v1 (Run_Tissue_name present)",
        identical(o_px$cell_id, px$cell_id) &&
        identical(o_px$X_px, px$X_px) && identical(o_px$Y_px, px$Y_px))
} else {
  cat("  --  Run_Tissue_name ABSENT: v1 silently mapped all ", nrow(o_px),
      " cells with this sample's placement; new maps ", nrow(px), "\n", sep = "")
  check("new result is a subset of v1's", all(px$cell_id %in% o_px$cell_id))
  check("the fix actually excluded other samples' cells", nrow(px) < nrow(o_px))
}

o_aff <- suppressWarnings(suppressMessages(v1$fit_mm_to_px(seu, SAMPLE_ID)))
n_aff <- suppressWarnings(suppressMessages(cosmx_fit_mm_affine(SAMPLE_ID, seu)))
check("affine ax/bx identical to v1",
      identical(o_aff$ax, n_aff$ax) && identical(o_aff$bx, n_aff$bx))
check("affine ay/by identical to v1",
      identical(o_aff$ay, n_aff$ay) && identical(o_aff$by, n_aff$by))
check("affine residuals and n identical to v1",
      identical(o_aff$resid_px, n_aff$resid_px) && identical(o_aff$n, n_aff$n))

## ---------------------------------------------------------------------------
## S4 -- the fit lands where the geometry says it must, and the guards fire.
## ---------------------------------------------------------------------------
cat("\n== S4 affine sanity ==\n")
expect <- 1 / COSMX_MM_PER_PX
cat(sprintf("  ax=%.3f bx=%.1f  ay=%.3f by=%.1f  resid=(%.3f, %.3f) n=%d\n",
            n_aff$ax, n_aff$bx, n_aff$ay, n_aff$by,
            n_aff$resid_px[1], n_aff$resid_px[2], n_aff$n))
check("both slopes are negative (x and y both mirrored vs mm)",
      n_aff$ax < 0 && n_aff$ay < 0)
check("slope magnitude within 1% of 1/COSMX_MM_PER_PX",
      abs(abs(n_aff$ax) - expect) / expect < 0.01 &&
      abs(abs(n_aff$ay) - expect) / expect < 0.01)
check("residual SD is sub-pixel", max(n_aff$resid_px) < 1)

## the affine must reproduce the exact route on this sample's own centroids
co  <- seu@images[[SAMPLE_ID]]$centroids
fit <- cosmx_mm_to_px(co@coords, n_aff)
ex  <- cosmx_cells_to_px(SAMPLE_ID, seu)
m   <- match(co@cells, ex$cell_id)
ok  <- !is.na(m)
dx  <- abs(fit$X_px[ok] - ex$X_px[m[ok]])
dy  <- abs(fit$Y_px[ok] - ex$Y_px[m[ok]])
cat(sprintf("  fitted vs exact: max |dx|=%.3f  max |dy|=%.3f px  (%d cells)\n",
            max(dx), max(dy), sum(ok)))
check("fitted affine agrees with the exact route to under 2 px",
      max(dx) < 2 && max(dy) < 2)

## ---------------------------------------------------------------------------
## S5 -- the new argument behaves
## ---------------------------------------------------------------------------
cat("\n== S5 sample_col ==\n")
check("sample_col = NA maps every cell",
      nrow(suppressWarnings(cosmx_cells_to_px(SAMPLE_ID, seu, sample_col = NA))) >=
      nrow(px))
check("an unknown sample_col stops",
      inherits(tryCatch(cosmx_cells_to_px(SAMPLE_ID, seu, sample_col = "nope"),
                        error = function(e) e), "error"))
check("an unknown image stops",
      inherits(tryCatch(cosmx_fit_mm_affine(SAMPLE_ID, seu, image = "nope"),
                        error = function(e) e), "error"))

cat("\n", n_pass, " checks passed\n", sep = "")
