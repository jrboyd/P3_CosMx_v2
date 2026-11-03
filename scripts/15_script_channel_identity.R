## ============================================================================
## Which marker is morphology channel "G"?
##
## The vendor names four of the five morphology channels (PanCK, Membrane, CD45,
## DAPI) and leaves the fifth as "G" -- in the OME-TIFF (Channel B/G/Y/R/U) and in
## the per-cell metadata (Mean.G / Max.G) alike. The question raised is whether it
## is CD68.
##
## Test: for every channel, rank all panel genes by how well their per-cell
## expression (log1p cp10k) tracks the channel's per-cell mean intensity
## (log). An antibody stain should be topped by its own gene's RNA. The four
## NAMED channels are the positive controls: CD45 must surface PTPRC and PanCK
## the keratins, or the method is not working and G's ranking means nothing.
##
## A random ~15% of cells per sample is read (spread over every FOV, unlike a
## head -n), InSituType cluster `i` (over-fixation defect, CLAUDE.md) is dropped,
## and genes detected in < 1% of cells are excluded so sparse genes cannot top a
## list on noise.
##
##   sbatch --mem=48G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##       scripts/15_script_channel_identity.R
## ============================================================================
suppressPackageStartupMessages(library(data.table))

FLAT     <- "input/CosMx/export1/flatFiles"
OUT_DIR  <- "output/15_imagej/channel_identity"
FRAC     <- 0.15
MIN_DET  <- 0.01
MIN_NCT  <- 50
TOP_N    <- 25
CHANNELS <- c("PanCK", "G", "Membrane", "CD45", "DAPI")
## Genes whose rank is reported for every channel regardless of where they fall:
## the myeloid set CD68 belongs to, plus each named channel's expected gene and
## the markers commonly offered in the CosMx optional slot.
WATCH <- c("CD68", "CD163", "CD14", "CSF1R", "LYZ", "C1QA", "C1QB", "MRC1",
           "PTPRC", "KRT8", "KRT18", "KRT19", "EPCAM", "B2M", "ATP1B3",
           "CD3E", "CD3D", "MS4A1", "CD19", "TNFRSF8", "CD8A", "CD4")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

res <- rbindlist(lapply(c("IB1", "RB1"), function(s) {
  message("== ", s)
  meta <- fread(file.path(FLAT, s, paste0(s, "_metadata_file.csv.gz")),
                select = c("fov", "cell_ID", "nCount_RNA",
                           "RNA_Volaric_02032026_Cell.Typing.InSituType.1_1_clusters",
                           paste0("Mean.", CHANNELS)))
  setnames(meta, "RNA_Volaric_02032026_Cell.Typing.InSituType.1_1_clusters", "insitu")

  ## Random subsample spread across all FOVs; awk's rand() seeded for reproducibility.
  cmd <- sprintf("zcat %s | awk 'BEGIN{srand(1)} NR==1 || rand() < %s'",
                 file.path(FLAT, s, paste0(s, "_exprMat_file.csv.gz")), FRAC)
  ex <- fread(cmd = cmd)
  message("  read ", nrow(ex), " cells x ", ncol(ex) - 2, " genes")

  ex <- merge(meta, ex, by = c("fov", "cell_ID"))
  n0 <- nrow(ex)
  ex <- ex[insitu != "i" & nCount_RNA >= MIN_NCT &
           Reduce(`&`, lapply(paste0("Mean.", CHANNELS), function(m) get(m) > 0))]
  message("  kept ", nrow(ex), " of ", n0, " (dropped cluster i, nCount < ", MIN_NCT,
          ", zero intensity)")

  genes <- setdiff(names(ex), c(names(meta)))
  genes <- genes[!grepl("^(Negative|SystemControl|NegPrb|FalseCode)", genes)]
  X <- as.matrix(ex[, ..genes])
  det <- colMeans(X > 0)
  X <- X[, det >= MIN_DET, drop = FALSE]
  message("  ", ncol(X), " genes detected in >= ", MIN_DET * 100, "% of cells")

  X <- log1p(X / rowSums(X) * 1e4)                 # rowSums over kept genes only
  Y <- log(as.matrix(ex[, paste0("Mean.", CHANNELS), with = FALSE]))
  R <- cor(X, Y)                                   # genes x channels, Pearson
  colnames(R) <- CHANNELS

  out <- rbindlist(lapply(CHANNELS, function(ch) {
    r <- R[, ch]
    data.table(sample = s, channel = ch, gene = names(r), r = unname(r),
               rank = frank(-r, ties.method = "min"), n_cells = nrow(ex),
               detect = unname(det[names(r)]))
  }))
  ## Channels' own mutual correlation, so a G that just tracks DAPI or autofluor
  ## is visible as such.
  cc <- as.data.table(cor(Y), keep.rownames = "channel")
  setnames(cc, c("channel", CHANNELS))
  fwrite(cc[, sample := s], file.path(OUT_DIR, paste0("channel_cor_", s, ".csv")))
  rm(X, ex); gc()
  out
}))

fwrite(res, file.path(OUT_DIR, "gene_channel_cor.csv.gz"))
top <- res[rank <= TOP_N][order(sample, channel, rank)]
fwrite(top, file.path(OUT_DIR, "top_genes_by_channel.csv"))
watch <- dcast(res[gene %in% WATCH], sample + gene ~ channel, value.var = "rank")
fwrite(watch, file.path(OUT_DIR, "watch_gene_ranks.csv"))

for (s in unique(res$sample)) for (ch in CHANNELS) {
  t <- top[sample == s & channel == ch][rank <= 12]
  cat(sprintf("%s %-8s %s\n", s, ch,
              paste(sprintf("%s(%.2f)", t$gene, t$r), collapse = " ")))
}
cat("\nwatch-gene ranks (out of ~", uniqueN(res$gene), " genes):\n", sep = "")
print(watch[order(sample, G)])
