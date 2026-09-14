# Is IB1's weak B compartment (7.5% of NBC_MBC calls are marker-B vs RB1's 35.3%)
# because the tissue genuinely HAS fewer B cells, or because label transfer is
# over-calling NBC_MBC on signal-free cells?
#
# Precision alone cannot tell these apart: the transfer must label every cell, so
# if true B cells are rare, ~30% of cells still get called NBC_MBC and precision
# collapses. Marker-defined B prevalence + recall separates them.
suppressPackageStartupMessages(library(Seurat))
md <- readRDS("output/04_merge/merged_01_merged.Rds")@meta.data
for (s in c("IB1", "RB1")) {
  m <- md[md$merge_sample == s, ]
  mb <- m$lineage == "B"; sb <- m$strongB; nb <- m$celltype_coarse == "NBC_MBC"
  inf <- !m$lineage %in% c("all_markers_0", "below_background", "ambiguous")
  cat("\n##########", s, "(", nrow(m), "cells )##########\n")
  cat(sprintf("marker-B prevalence          : %6.2f%% (n=%d)\n", 100*mean(mb), sum(mb)))
  cat(sprintf("  ... among informative cells: %6.2f%% (n=%d informative)\n",
              100*mean(mb[inf]), sum(inf)))
  cat(sprintf("strongB prevalence           : %6.2f%% (n=%d)\n", 100*mean(sb), sum(sb)))
  cat(sprintf("NBC_MBC called               : %6.2f%% (n=%d)\n", 100*mean(nb), sum(nb)))
  cat(sprintf("RECALL  (marker-B called NBC_MBC)     : %6.2f%%\n", 100*mean(nb[mb])))
  cat(sprintf("RECALL  (strongB   called NBC_MBC)    : %6.2f%%\n", 100*mean(nb[sb])))
  cat(sprintf("PRECISION (NBC_MBC that are marker-B) : %6.2f%%\n", 100*mean(mb[nb])))
  cat("where marker-B cells actually got called (%):\n")
  print(round(100*sort(table(m$celltype_coarse[mb]), decreasing=TRUE)/sum(mb), 1))
  cat("where strongB cells got called (%):\n")
  print(round(100*sort(table(m$celltype_coarse[sb]), decreasing=TRUE)/sum(sb), 1))
}
cat("\nREAD: similar marker-B prevalence + high recall + low precision => transfer over-calls.\n")
cat("      much lower marker-B prevalence in IB1 => genuinely fewer B cells.\n\nOK\n")
