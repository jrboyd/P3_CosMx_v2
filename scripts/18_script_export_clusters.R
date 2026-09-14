## ============================================================================
## 18_script_export_clusters.R -- export per-cell cluster assignments to CSV.
##
## WHY THIS EXISTS. The step-09b marker audit labels CLUSTERS, and the labels live in
## output/09b_marker_relabel/. But the per-cell cluster ASSIGNMENT lives only inside the
## Seurat object, so every script that wants "which cluster is this cell in" has had to
## load a 1.3 GB .Rds. Step 18 needs it for ~1M cells joined against transcripts, so it
## gets exported once, here, and read as a CSV from then on.
##
## inSitu_clusters is exported too even though output/16_regions/cells_meta_all.csv.gz
## already carries it -- as a join check. If the two disagree the export is wrong.
##
## Usage: sbatch --mem=64G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/18_script_export_clusters.R IB1
## ============================================================================
suppressPackageStartupMessages({ library(data.table); setDTthreads(1) })
args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"
IN_RDS <- sprintf("output/03_filter_and_process/%s_03_HCA_cell_typed.Rds", SAMPLE)
OUT    <- "output/18_leakage"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(IN_RDS)) stop("missing ", IN_RDS, call. = FALSE)

cat("reading ", IN_RDS, "\n", sep = "")
obj <- readRDS(IN_RDS)
md  <- obj@meta.data
cat("cells: ", nrow(md), "\n", sep = "")
cat("\n=== all meta.data columns ===\n"); print(colnames(md))

WANT <- c("pearson_clusters_low", "pearson_clusters_high", "inSitu_clusters")
have <- intersect(WANT, colnames(md))
miss <- setdiff(WANT, colnames(md))
if (length(miss)) cat("\nNOT PRESENT (exported as NA): ", paste(miss, collapse = ", "), "\n", sep = "")
if (!length(have)) stop("none of the wanted cluster columns are in this object", call. = FALSE)

d <- data.table(cell_id = rownames(md))
for (k in WANT) d[[k]] <- if (k %in% have) as.character(md[[k]]) else NA_character_

for (k in have) {
  cat("\n=== ", k, " ===\n", sep = "")
  print(d[, .N, by = k][order(-N)])
}

f <- file.path(OUT, sprintf("%s_clusters.csv.gz", SAMPLE))
fwrite(d, f)
cat("\nwrote ", f, " (", nrow(d), " cells)\n", sep = "")
