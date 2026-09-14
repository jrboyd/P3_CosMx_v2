#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------------------
## step 25a -- FOV-LEVEL UNCERTAINTY ON COMPOSITION.  Closes gap G8.
##
## THE GAP. `output/19_celltypes/composition_by_sample.csv` is a point estimate per label:
## "IB1 is 23.6% B cells". At n = 1 patient and two sections that number has no interval on
## it, and a bar chart of point estimates invites the reader to compare 23.6% against 71.4%
## as though the difference were the only thing uncertain. The honest presentation attaches
## the spread ACROSS FIELDS OF VIEW, which is the project's unit of analysis everywhere else
## ("region is the unit of every p-value").
##
## WHAT THIS IS AND IS NOT. The FOV interquartile range describes how much composition VARIES
## ACROSS THE SECTION. It is a descriptive spread, not a confidence interval on a population
## mean -- there is one patient, so no such population exists and none is implied. A wide IQR
## means the tissue is heterogeneous, not that the point estimate is unreliable.
##
## CELLS. `confidence == "high"` only, per step 19's statement rule: rescued labels are
## validated by extrapolation from their stratum and should not carry a composition claim.
## FOV is parsed from the cell id, where it is the second-to-last token -- sanctioned by
## CLAUDE.md. The SAMPLE is taken from the `samp` COLUMN and never from the id prefix, because
## `c_2_` is IB1 and `c_1_` is RB1 and guessing it the obvious way round silently swaps every
## row.
##
## THE HODGKIN CHEMOKINE NICHE IS KEPT BUT FLAGGED. It is a REGION label, not a lineage
## (CLAUDE.md), and step 22's snapshots showed what happens when that is forgotten. Dropping
## it would renormalise every other share and quietly inflate them, so it is retained in the
## table with `is_region = TRUE` and the figure must render it distinctly and say so.
##
## Usage: sbatch --mem=32G --time=0:30:00 scripts/sbatch/run_rscript.sh \
##          scripts/25_script_composition_fov.R
## ---------------------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table)})
setDTthreads(1); options(width = 200)

OUT <- "output/25_figures"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sc <- function(f, x) { fwrite(x, file.path(OUT, f)); cat("wrote", file.path(OUT, f), "\n") }
MIN_CELLS_FOV <- 50L      # below this a per-FOV share is noise, not heterogeneity
NICHE <- "Hodgkin chemokine niche"

cat("========================================================================\n")
cat("step 25a -- FOV-level composition spread (G8)\n")
cat("========================================================================\n")

CE <- fread("output/19_celltypes/all_cells_celltypes.csv.gz",
            select = c("cell_id","samp","confidence","celltype_cluster"))
cat(sprintf("\nall labelled cells: %s\n", format(nrow(CE), big.mark = ",")))
CE <- CE[confidence == "high" & !is.na(celltype_cluster) & celltype_cluster != ""]
## FOV is the second-to-last token of c_<slide>_<fov>_<cell>.
CE[, fov := as.integer(tstrsplit(cell_id, "_", keep = 3L)[[1]])]
stopifnot(!anyNA(CE$fov))
cat(sprintf("high-confidence, labelled: %s in %d FOVs\n",
            format(nrow(CE), big.mark = ","), uniqueN(CE[, .(samp, fov)])))

## per FOV: share of each label among that FOV's labelled cells
tot <- CE[, .(n_fov_cells = .N), by = .(samp, fov)][n_fov_cells >= MIN_CELLS_FOV]
cat(sprintf("FOVs with >= %d labelled cells: %d of %d\n",
            MIN_CELLS_FOV, nrow(tot), uniqueN(CE[, .(samp, fov)])))
PF <- CE[tot, on = .(samp, fov)][, .(n = .N), by = .(samp, fov, celltype_cluster, n_fov_cells)]
## a label absent from an FOV is a share of 0, not a missing row
grid <- CJ(samp = unique(PF$samp), fov = unique(PF$fov),
           celltype_cluster = unique(PF$celltype_cluster), unique = TRUE)
grid <- grid[tot[, .(samp, fov, n_fov_cells)], on = .(samp, fov), nomatch = 0L]
PF <- PF[grid, on = .(samp, fov, celltype_cluster, n_fov_cells)][is.na(n), n := 0L]
PF[, share := n / n_fov_cells]
sc("composition_by_fov.csv.gz", PF)

S <- PF[, .(fovs = .N, cells = sum(n),
            share_med = median(share), share_q1 = quantile(share, .25),
            share_q3 = quantile(share, .75),
            share_p10 = quantile(share, .10), share_p90 = quantile(share, .90),
            share_pooled = sum(n)/sum(n_fov_cells)),
        by = .(samp, celltype_cluster)]
S[, is_region := celltype_cluster == NICHE]
S <- S[order(samp, -share_pooled)]
cat("\n=== composition with FOV-level spread\n")
print(S[, .(samp, celltype_cluster, fovs, cells,
            pooled_pct = round(100*share_pooled, 2),
            med_pct = round(100*share_med, 2),
            IQR_pct = sprintf("%.2f-%.2f", 100*share_q1, 100*share_q3),
            is_region)])
sc("composition_fov_summary.csv", S)

cat("\n=== how much does the pooled estimate differ from the FOV median?\n")
cat("    (a large gap means a few dense FOVs carry the pooled number)\n")
print(S[, .(samp, celltype_cluster, pooled_pct = round(100*share_pooled,2),
            med_pct = round(100*share_med,2),
            diff_pp = round(100*(share_pooled - share_med), 2))][order(-abs(diff_pp))][1:10])
writeLines(format(Sys.time()), file.path(OUT, "composition_fov.done"))
cat("\nOK\n")
