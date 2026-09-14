## ============================================================================
## cluster_labels.R -- the canonical cluster-label lookup for this project.
##
## Downstream scripts must NOT read output/09_hierannot/*/hierannot_cluster_labels.csv
## directly. HierAnnot's annot_path is wrong for a substantial minority of clusters on
## this data (see scripts/09b_script_marker_relabel.R for the audit and the evidence),
## and the errors are not subtle -- stroma called macrophage, endothelium called
## neutrophil, pDC called plasma cell. Read through here instead, so every script gets
## the marker-verified label and the same one.
##
##   source("scripts/lib/cluster_labels.R")
##   L   <- cluster_labels("IB1", "pearson_clusters_low")   # data.table, one row/cluster
##   vec <- cluster_label_vec("IB1", "pearson_clusters_low") # named chr, for plot labels
## ============================================================================

suppressPackageStartupMessages(library(data.table))

CLUSTER_LABEL_DIR <- "output/09b_marker_relabel"

## Full table: cluster_id, label, status, hierannot_path, agrees, evidence, and the
## scoring columns behind the call.
cluster_labels <- function(sample, clustering = NULL, path = CLUSTER_LABEL_DIR) {
  f <- file.path(path, sample, "cluster_labels_final.csv")
  if (!file.exists(f))
    stop("no marker-verified labels for ", sample, " at ", f,
         "\n  run: Rscript scripts/09b_script_marker_relabel.R ", sample, call. = FALSE)
  d <- fread(f)
  if (!is.null(clustering)) {
    want <- clustering                 # avoid the column/arg name collision in i
    d <- d[d$clustering == want]
    if (!nrow(d)) stop("no rows for clustering '", clustering, "' in ", f, call. = FALSE)
  }
  d[]
}

## Named character vector keyed by cluster_id, for axis/facet labels. `with_status`
## appends a marker on anything not confidently resolved, so a figure cannot quietly
## present a shaky call as settled.
cluster_label_vec <- function(sample, clustering, with_status = TRUE,
                              path = CLUSTER_LABEL_DIR) {
  d <- cluster_labels(sample, clustering, path)
  lab <- d$label
  if (with_status) {
    ## Newline, not " (...)": these land in facet strips and axis text, where a long
    ## single line is silently clipped -- "GC B cell (weak_support)" rendered as
    ## "B cell (weak_suppor" in a 9-panel heatmap.
    lab <- ifelse(d$status == "confident", lab, paste0(lab, "\n(", d$status, ")"))
  }
  setNames(lab, as.character(d$cluster_id))
}
