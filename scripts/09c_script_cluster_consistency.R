## ============================================================================
## 09c_script_cluster_consistency.R -- do the per-sample clusterings and the merged
## clustering name the same cells the same way?
##
## WHY THIS EXISTS. Step 16 decided that cross-sample cell type must come from the argmax
## marker caller rather than from clusters, and recorded the reason in a header comment:
## "mapped through sample_pearson_clusters_low, which is the only honest route". That was an
## assumption. This script tests it, and bounds any future cluster-level cross-sample claim.
##
## THREE CLUSTERINGS, AND ONLY ONE IS DIRECTLY COMPARABLE.
##   pearson_clusters_low/high -- derived INDEPENDENTLY per object, so IB1's cluster 3 and
##     RB1's cluster 3 are unrelated by construction. Comparison must go through the merged
##     object's `sample_pearson_clusters_*` columns, which record which per-sample cluster
##     each cell came from.
##   inSitu_clusters -- the VENDOR InSituType call, assigned once and carried into every
##     object, so ids "a".."j" mean the same thing everywhere. This is the only clustering
##     where "does IB1's cluster c mean what RB1's cluster c means?" is directly answerable,
##     and it is answered by comparing the three independent marker audits of it.
##
## SAMPLE COMES FROM `merge_sample`, NEVER FROM THE CELL-ID PREFIX. The prefix mapping is
## counter-intuitive on this data -- IB1 is `c_2_` and RB1 is `c_1_` -- and guessing it the
## other way round silently swaps every row of the output. Verified against
## output/16_regions/cells_meta_all.csv.gz, which carries both columns.
##
## Usage: Rscript scripts/09c_script_cluster_consistency.R
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); setDTthreads(1)
})
if (!file.exists("scripts/lib/cluster_labels.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/cluster_labels.R")

IN_RDS <- "output/04_merge/merged_processed.Rds"
OUT    <- "output/09c_cluster_consistency"
FIG    <- file.path(OUT, "figures")
MAP_F  <- file.path(OUT, "merged_cluster_map.csv.gz")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n, d) { fwrite(d, file.path(OUT, n)); cat("wrote ", file.path(OUT, n), "\n", sep="") }

SURF <- "#fcfcfb"; INK <- "#0b0b0b"; GRID <- "#e1e0d9"
DIV_LO <- "#2c7bb6"; DIV_HI <- "#d7191c"          # blue -> red, red is HIGH (CLAUDE.md)
theme_viz <- function(base = 10) theme_minimal(base_size = base) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid = element_line(colour = GRID, linewidth = .25),
        plot.title = element_text(face = "bold", colour = INK),
        plot.caption = element_text(colour = "#898781", size = base - 2.5))
sf <- function(n, p, w, h) { ggsave(file.path(FIG, n), p, width = w, height = h, dpi = 200); cat("wrote ", file.path(FIG, n), "\n", sep="") }

cat("============================================================\n")
cat("step 09c -- cluster label consistency, per-sample vs merged\n")
cat("============================================================\n")

## ---- 1. the per-cell mapping, cached so the Rds is read once -----------------------
NEED <- c("merge_sample", "inSitu_clusters", "pearson_clusters_low", "pearson_clusters_high",
          "sample_pearson_clusters_low", "sample_pearson_clusters_high")
if (file.exists(MAP_F)) {
  cat("reading cached ", MAP_F, "\n", sep="")
  MP <- fread(MAP_F, colClasses = "character")
} else {
  if (!file.exists(IN_RDS)) stop("missing ", IN_RDS, call. = FALSE)
  cat("reading ", IN_RDS, "\n", sep="")
  obj <- readRDS(IN_RDS)
  miss <- setdiff(NEED, colnames(obj@meta.data))
  if (length(miss)) stop("merged meta.data missing: ", paste(miss, collapse=", "), call.=FALSE)
  MP <- as.data.table(lapply(obj@meta.data[, NEED, drop = FALSE], as.character))
  MP[, cell_id := rownames(obj@meta.data)]
  rm(obj); invisible(gc())
  fwrite(MP, MAP_F); cat("wrote ", MAP_F, "\n", sep="")
}
cat(sprintf("cells %d | samples %s\n", nrow(MP),
            paste(sprintf("%s=%d", names(table(MP$merge_sample)), table(MP$merge_sample)),
                  collapse=" ")))

## The prefix check that this script refuses to rely on, asserted rather than assumed.
pfx <- MP[, .N, by = .(merge_sample, p = sub("^(c_[0-9]+)_.*$", "\\1", cell_id))][order(-N)]
cat("\nsample <-> cell-id prefix (from the object, not guessed):\n"); print(pfx)

## ---- 2. audited labels for all three objects --------------------------------------
LAB <- rbindlist(lapply(c("IB1","RB1","merged"), function(S) {
  d <- tryCatch(cluster_labels(S), error = function(e) NULL)
  if (is.null(d)) { cat("NOTE: no audit for ", S, "\n", sep=""); return(NULL) }
  d[, .(obj = S, clustering, cluster_id = as.character(cluster_id), label, status,
        cells, pct_cells, top_z, pct_defect, best_supported, n_family_genes)]
}), fill = TRUE)
cat("\naudits loaded:\n"); print(LAB[, .N, by = .(obj, clustering)])

lab_of <- function(S, cc) {
  d <- LAB[obj == S & clustering == cc]
  setNames(sprintf("%s%s", d$label, ifelse(d$status == "confident", "", paste0(" [", d$status, "]"))),
           d$cluster_id)
}

## ---- 3. pearson: merged cluster x per-sample audited label ------------------------
## `sample_pearson_clusters_*` is "IB1_3" / "RB1_0"; split it to recover the per-sample id
## and look the label up in THAT sample's audit.
CONF <- list(); AGREE <- list()
for (res in c("low","high")) {
  mcol <- paste0("pearson_clusters_", res)
  scol <- paste0("sample_pearson_clusters_", res)
  d <- MP[, .(merged_cl = get(mcol), spl = get(scol), samp = merge_sample)]
  d[, samp_cl := sub("^[A-Za-z0-9]+_", "", spl)]
  ## the sample recorded in the id must match merge_sample, or the join is meaningless
  d[, samp_from_id := sub("_.*$", "", spl)]
  bad <- d[samp_from_id != samp, .N]
  if (bad) stop("sample_", res, ": ", bad, " cells where the id's sample disagrees with merge_sample",
                call. = FALSE)
  d[, samp_label := ifelse(samp == "IB1", lab_of("IB1", mcol)[samp_cl],
                                          lab_of("RB1", mcol)[samp_cl])]
  d[, merged_label := lab_of("merged", mcol)[merged_cl]]

  ## confusion: one row per merged cluster x sample x per-sample label
  ## `resolution = res` cannot go in `by`: data.table requires every by-item to be as long
  ## as the table, and a loop scalar is length 1. Add it as a column afterwards.
  cf <- d[, .N, by = .(merged_cl, merged_label, samp, samp_cl, samp_label)]
  cf[, resolution := res]
  cf[, merged_cl_n := sum(N), by = merged_cl]
  cf[, share_of_merged_cl := round(N / merged_cl_n, 4)]
  setorder(cf, merged_cl, -N)
  CONF[[res]] <- cf

  ## per-merged-cluster: the dominant per-sample label WITHIN each sample, and whether the
  ## two samples' dominant labels name the same cell type. A merged cluster drawn almost
  ## entirely from one sample is the batch-effect mode; one whose two pluralities disagree
  ## is the mislabelling mode. They are different problems, so both are reported.
  dom <- cf[, .SD[which.max(N)], by = .(merged_cl, samp)]
  ag <- dcast(dom, resolution + merged_cl + merged_label ~ samp,
              value.var = c("samp_label","N"))
  setnames(ag, c("samp_label_IB1","samp_label_RB1","N_IB1","N_RB1"),
           c("IB1_dominant","RB1_dominant","IB1_n","RB1_n"), skip_absent = TRUE)
  for (cc in c("IB1_n","RB1_n")) if (cc %in% names(ag)) ag[is.na(get(cc)), (cc) := 0L]
  tot <- d[, .N, by = .(merged_cl, samp)]
  ag <- merge(ag, dcast(tot, merged_cl ~ samp, value.var = "N", fill = 0L),
              by = "merged_cl", all.x = TRUE, suffixes = c("","_tot"))
  setnames(ag, c("IB1","RB1"), c("IB1_cells","RB1_cells"), skip_absent = TRUE)
  ag[, one_sample_frac := round(pmax(IB1_cells, RB1_cells) / pmax(IB1_cells + RB1_cells, 1), 3)]
  ## strip the "[status]" suffix before comparing: a status difference is not a type disagreement
  base <- function(x) trimws(sub("\\[.*\\]$", "", x))
  ag[, types_agree := base(IB1_dominant) == base(RB1_dominant)]
  ag[, mode := fifelse(one_sample_frac >= 0.90, "single-sample cluster",
               fifelse(is.na(types_agree), "one sample only",
               fifelse(types_agree, "agrees", "DISAGREES")))]
  setorder(ag, -one_sample_frac)
  AGREE[[res]] <- ag

  cat("\n=== ", mcol, ": merged cluster vs each sample's dominant audited label\n", sep="")
  print(ag[, .(merged_cl, merged_label, IB1_dominant, RB1_dominant,
               IB1_cells, RB1_cells, one_sample_frac, mode)])
}
sc("confusion_pearson.csv", rbindlist(CONF, fill = TRUE))
sc("agreement_pearson.csv", rbindlist(AGREE, fill = TRUE))

## ---- 4. inSitu: three independent audits of the SAME vendor ids -------------------
IS <- LAB[clustering == "inSitu_clusters"]
if (nrow(IS)) {
  w <- dcast(IS, cluster_id ~ obj, value.var = c("label","status","pct_cells"))
  base <- function(x) trimws(sub("\\[.*\\]$", "", x))
  w[, all_three_agree := base(label_IB1) == base(label_RB1) &
                         base(label_IB1) == base(label_merged)]
  w[, samples_agree := base(label_IB1) == base(label_RB1)]
  setorder(w, cluster_id)
  cat("\n=== inSitu_clusters: the only clustering whose ids mean the same thing everywhere\n")
  print(w[, .(cluster_id, label_IB1, label_RB1, label_merged, samples_agree, all_three_agree)])
  cat(sprintf("\nvendor clusters where BOTH SAMPLES agree: %d of %d\n",
              sum(w$samples_agree, na.rm = TRUE), nrow(w)))
  sc("agreement_insitu.csv", w)

  ## Why each disagreement happened -- the abundance cap and single-gene families are the
  ## two mechanisms, and both are visible in the audit's own columns.
  dis <- IS[cluster_id %in% w[samples_agree == FALSE, cluster_id]]
  if (nrow(dis)) {
    cat("\n=== the disagreeing vendor clusters, with the audit's own reasons\n")
    print(dis[order(cluster_id, obj), .(cluster_id, obj, label, status, top_z,
                                        best_supported, n_family_genes, pct_cells)])
  }

  hm <- melt(w, id.vars = "cluster_id",
             measure.vars = patterns("^pct_cells_"), variable.name = "obj",
             value.name = "pct")[, obj := sub("pct_cells_", "", obj)]
  lv <- w[, setNames(sprintf("%s\n%s", cluster_id, label_merged), cluster_id)]
  p <- ggplot(hm, aes(cluster_id, obj, fill = as.numeric(pct))) + geom_tile() +
    geom_text(aes(label = sprintf("%.1f%%", as.numeric(pct))), size = 2.6, colour = INK) +
    scale_fill_gradient(low = SURF, high = DIV_HI, name = "% of object's scored cells") +
    scale_x_discrete(labels = lv, expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
    labs(title = "Vendor InSituType clusters: prevalence across the three audits",
         subtitle = paste("Same ids in every object, so prevalence is comparable.",
                          "Column labels are the merged audit's call.\nDisagreements between",
                          "the three audits are listed in agreement_insitu.csv."),
         x = NULL, y = NULL) +
    theme_viz(10) + theme(panel.grid = element_blank())
  sf("fig_09c_insitu_prevalence.png", p, 11, 3.4)
}

## ---- 5. the headline -------------------------------------------------------------
cat("\n============================================================\nSUMMARY\n")
cat("============================================================\n")
for (res in c("low","high")) {
  a <- AGREE[[res]]
  cat(sprintf("pearson_clusters_%-5s %2d merged clusters | %2d single-sample (>=90%% one sample) | %2d agree | %2d DISAGREE\n",
              res, nrow(a), sum(a$mode == "single-sample cluster"),
              sum(a$mode == "agrees"), sum(a$mode == "DISAGREES")))
}
if (nrow(IS)) cat(sprintf("inSitu_clusters       %2d vendor clusters | %2d agree between samples | %2d disagree\n",
                          nrow(w), sum(w$samples_agree, na.rm=TRUE),
                          sum(!w$samples_agree, na.rm=TRUE)))
cat("\nRead agreement_pearson.csv `mode` first: 'single-sample cluster' is a batch effect,\n",
    "'DISAGREES' is a labelling conflict. They are different failures.\n", sep="")
cat("============================================================\nOK\n")
