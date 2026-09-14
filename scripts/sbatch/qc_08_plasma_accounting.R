## Diagnostic for step 08: where do plasma cells actually end up?
## Cross-tabs the step-07 tier / step-08 group assignment against the independent
## marker-based `lineage` call from step 01, for the IB1 r40um pool.
## Answers "how much of each compartment is plasma, and what is the rest".
##   sbatch --mem=32G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##     scripts/sbatch/qc_08_plasma_accounting.R
suppressPackageStartupMessages({library(Seurat); library(data.table); setDTthreads(1)})

S <- "IB1"; R <- 40
obj <- readRDS(sprintf("output/01_import_v1/seuratObject_%s_BPCells_full.Rds", S))
md  <- as.data.table(obj@meta.data[, c("lineage","strongB","nCount_RNA"), drop = FALSE])
md[, cell_id := rownames(obj@meta.data)]
rm(obj); invisible(gc())

lc  <- fread(sprintf("output/07_lightchain/%s/lightchain_cells_r%dum.csv", S, R))
nbt <- fread(sprintf("output/07_lightchain/%s/lightchain_neighborhoods_r%dum.csv", S, R))[
         , .(nb, class, artifact_suspect)]
x <- merge(lc, nbt, by = "nb", all.x = TRUE)
x[, grp := fifelse(tier != "B_clean", tier,
          fifelse(is.na(class),                          "B_clean: diffuse",
          fifelse(artifact_suspect %in% TRUE,            "B_clean: dropped",
          fifelse(class == "kappa-restricted (FL-like)", "B_clean: FL-like",
                                                         "B_clean: clustered"))))]
d <- merge(x[, .(cell_id, grp, inSitu_clusters)], md, by = "cell_id")
d <- d[inSitu_clusters != "i"]          # the primary filter step 08 uses

ct <- d[, .N, by = .(grp, lineage)]
ct[, pct := round(100 * N / sum(N), 1), by = grp]
w <- dcast(ct, grp ~ lineage, value.var = "pct", fill = 0)
w <- merge(d[, .(cells = .N, pct_strongB = round(100*mean(strongB %in% TRUE), 1),
                 med_nCount = as.numeric(median(nCount_RNA))), by = grp], w, by = "grp")
cat("=== lineage composition of each compartment (% of its own cells), cluster i excluded\n")
print(w[order(-B)])

cat("\n=== the same for ALL IB1 cells outside the light-chain pool, as a baseline\n")
bg <- md[!cell_id %in% x$cell_id]
print(bg[, .(cells = .N, pct = round(100*.N/nrow(bg), 1)), by = lineage][order(-cells)])

fwrite(w, "output/08_ib1_clustered_vs_diffuse/plasma_accounting_lineage.csv")
cat("\nwrote output/08_ib1_clustered_vs_diffuse/plasma_accounting_lineage.csv\n")
cat("OK\n")
