## Read the CSVs a reader sends back, unblind them against the key, summarise.
##
##   Rscript --vanilla pathology_review_app/ingest/ingest_reviews.R [dir_or_csv ...]
##   (default: every CSV in pathology_review_app/returns/; small data, fine inline)
##
## Every export holds EVERYTHING that reader had saved at export time, so a later export
## supersedes an earlier one; records are de-duplicated by record_id, then the latest
## revision per (reader, bundle, item) is the answer. All revisions are kept in the long
## table, so a change of mind is visible.
##
## Writes pathology_review_app/output/reviews/:
##   records_all.csv      every revision, unblinded
##   answers_latest.csv   one row per reader x item, unblinded
##   summary_by_group.csv per reader x sample x kind x stratum: n, adequate rate, yes-rate
##                        per checklist item (among adequate), mean score, share scored >= 4
##   repeat_agreement.csv per reader x repeated cell: original vs repeat answers
##
## Reading order matters for interpretation: the comparison that means something is
## candidate vs matched CONTROL within a sample, and IB1 vs RB1 at matched brightness.
## Never pool across samples without saying so (CLAUDE.md fact 5: the slides differ).
suppressMessages(library(data.table))
ROOT <- "pathology_review_app"
args <- commandArgs(TRUE)
if (!length(args)) args <- file.path(ROOT, "returns")
files <- unlist(lapply(args, function(a) if (dir.exists(a))
  list.files(a, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE) else a))
if (!length(files)) stop("no export CSVs found in: ", paste(args, collapse = ", "))
OUT <- file.path(ROOT, "output", "reviews"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rec <- rbindlist(lapply(files, function(f) fread(f, colClasses = "character", encoding = "UTF-8")[, src_file := basename(f)]),
                 fill = TRUE)
rec <- unique(rec, by = "record_id")
rec[, `:=`(revision = as.integer(revision), hrs_confidence = as.numeric(hrs_confidence),
           view_seconds = as.numeric(view_seconds))]
cat(sprintf("%d export file(s), %d unique records, readers: %s, bundles: %s\n", length(files), nrow(rec),
            paste(unique(rec$reader), collapse = ", "), paste(unique(rec$bundle_id), collapse = ", ")))

key <- rbindlist(lapply(unique(rec$bundle_id), function(b) {
  kf <- file.path(ROOT, "keys", paste0(b, "_key.csv"))
  if (!file.exists(kf)) stop("no key for bundle ", b, " (expected ", kf, ")")
  fread(kf, colClasses = c(item_id = "character"))[, bundle_id := b]
}))
key[, seq := NULL]
rec <- merge(rec, key, by = c("bundle_id", "item_id"), all.x = TRUE)
stopifnot(!anyNA(rec$cell_id))
setorder(rec, reader, bundle_id, saved_at_utc)
fwrite(rec, file.path(OUT, "records_all.csv"))

lat <- rec[, .SD[.N], by = .(reader, bundle_id, item_id)]
fwrite(lat, file.path(OUT, "answers_latest.csv"))

## ---- summary: originals only (repeats are for agreement, not counted twice) --------
chk <- intersect(c("single_cell", "large_cell", "multinuc", "nucleolus", "pale_chrom", "cd45_dim", "rosette"),
                 names(lat))
o <- lat[repeat_of == "" | is.na(repeat_of)]
summ <- o[, c(list(n = .N, adequate = mean(adequate == "yes")),
              setNames(lapply(chk, function(k) { a <- get(k)[adequate == "yes"]; if (length(a)) mean(a == "yes") else NA_real_ }),
                       paste0("yes_", chk)),
              list(mean_score = mean(hrs_confidence, na.rm = TRUE),
                   score_ge4 = mean(hrs_confidence >= 4, na.rm = TRUE),
                   median_seconds = median(view_seconds, na.rm = TRUE))),
          by = .(reader, samp, kind, stratum)]
setorder(summ, reader, samp, kind, stratum)
fwrite(summ, file.path(OUT, "summary_by_group.csv"))

## ---- intra-reader agreement on planned repeats -------------------------------------
rp <- lat[repeat_of != "" & !is.na(repeat_of)]
if (nrow(rp)) {
  ag <- merge(rp, lat, by.x = c("reader", "bundle_id", "repeat_of"), by.y = c("reader", "bundle_id", "item_id"),
              suffixes = c("_rep", "_orig"))
  cols <- c("adequate", chk, "hrs_confidence")
  rpt <- ag[, c(list(cell_id = cell_id_rep), setNames(lapply(cols, function(k)
    paste(get(paste0(k, "_orig")), get(paste0(k, "_rep")), sep = " -> ")), cols)), by = .(reader, repeat_of)]
  fwrite(rpt, file.path(OUT, "repeat_agreement.csv"))
  both <- ag[adequate_orig == "yes" & adequate_rep == "yes"]
  cat(sprintf("repeats: %d answered twice; score exact agreement %.0f%%, within 1 point %.0f%% (n = %d both adequate)\n",
              nrow(ag), 100 * mean(both$hrs_confidence_orig == both$hrs_confidence_rep),
              100 * mean(abs(both$hrs_confidence_orig - both$hrs_confidence_rep) <= 1), nrow(both)))
}
print(summ[, .(reader, samp, kind, stratum, n, adequate = round(adequate, 2),
               mean_score = round(mean_score, 2), score_ge4 = round(score_ge4, 2))])
cat("wrote", OUT, "\n")
