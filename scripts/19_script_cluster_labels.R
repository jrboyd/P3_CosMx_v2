## ============================================================================
## 19_script_cluster_labels.R -- name the high-confidence clusters. Guilt by association.
##
## A cluster is named by the argmax calls of the cells inside it that DID resolve. Cells that
## resolved nothing inherit that name. That is the "guilt by association" concept, made
## literal: marker evidence names the cluster, co-clustering spreads the name.
##
## LABELS COME FROM ABSOLUTE COMPOSITION, NOT ACROSS-CLUSTER z. Step 09b scored clusters by
## z across clusters within a sample (09b:199, `Z <- t(scale(t(log1p(CP))))`), which
## collapses when one lineage dominates -- RB1 is ~60% B and its largest cluster, 836,909
## cells, audits `low-information` at 3% defect while its own table records best_supported =
## B cell with a 0.75 margin. The high-confidence object is 62.6% RB1 and B is 50.5% of its
## resolved cells, so a relative score would fail here for the same reason. Modal argmax
## share is absolute and cannot collapse.
##
## THE B LABEL IS "B cell", NOT "B cell, clonal". The aggregate evidence for one clone is
## strong -- kappa (fold 1.31, CV 0.073) and IgM (1.45, 0.073) are the least variable
## measurements in the object, and 178/178 RB1 neighbourhoods are kappa-restricted -- but
## that is a statement about the COMPARTMENT, not a per-cell one. Baking "clonal" into the
## name asserts it of every one of ~1.08M cells, including IB1's reactive B cells, which no
## per-cell measurement here supports. The clonality evidence lives in the cluster table
## and the substructure outputs, where it can be read with its denominator attached.
##
## B GETS A SECOND STAGE, AND IT COLLAPSES RATHER THAN SPLITS. 20 of 27 clusters have modal
## call B. The module analysis (19_script_b_substructure.R) shows why, and it is not that
## there are 20 B cell types: normalised to B content, kappa (fold 1.31, CV 0.073) and IgM
## (fold 1.45, CV 0.073) are the LEAST variable measurements in the whole object. Every B
## cluster is IgM+, kappa-restricted and BCL2-high -- one clone, matching the 178/178
## kappa-restricted neighbourhoods already on record. Eleven of the twenty sit within 25% of
## the median on every module. So the second stage MERGES the monomorphic mass into a single
## clonal-B label and pulls out only the clusters that genuinely depart:
##   naive/mantle   naive_mantle ratio > NAIVE_MULT x the monomorphic-core median
##   otherwise      clonal B
## Germinal-centre substructure is deliberately NOT assigned: dark/light zone vary only
## 2.9-5.2 fold and the variation traces to clusters with low B content, i.e. small
## numerators over small denominators, not biology.
##
## A CLUSTER TOO MIXED TO NAME IS RECORDED AS MIXED, NOT NAMED ANYWAY. That is the whole
## lesson of 09c: a cluster acquires a lineage label because it happened to co-cluster, and
## the label then travels. `mixed` clusters keep their module profile in the table so the
## evidence is visible, but they carry no lineage name and no cell-type statement may use
## them.
##
## Usage: sbatch --mem=32G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/19_script_cluster_labels.R [high|low]
## ============================================================================
suppressPackageStartupMessages({ library(data.table); setDTthreads(2) })
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/lineage.R")

args <- commandArgs(trailingOnly = TRUE)
RES  <- if (length(args) >= 1) args[1] else "high"
CL   <- paste0("hc_clusters_", RES)
OUT  <- "output/19_celltypes"
HC_F <- file.path(OUT, "hc_cells.csv.gz")
CT_F <- "output/16_celltypes/cell_celltypes.csv.gz"
MOD_F<- file.path(OUT, sprintf("b_substructure_modules_%s.csv", RES))
for (f in c(HC_F, CT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

JUNK_FRAC  <- 0.15   # a cluster resolving below this is junk at any size
## TWO ROUTES TO CONFIDENT, because purity alone is the wrong yardstick for a minor lineage.
## The object is ~50% B among resolved cells, so a B cluster can reach purity 0.82 while an
## NK cluster structurally cannot: NK is only ~4.8% of resolved cells, and a cluster at 33.5%
## NK -- seven times its baseline -- was being called `mixed` and losing the label entirely.
## Enrichment is capped at 1/baseline, so it can never certify the DOMINANT family (B tops
## out near 2x); purity can never certify a minor one. Each route covers the other's blind
## spot, so a cluster qualifies on either.
CONF_PUR   <- 0.50; CONF_MRG <- 0.20   # route 1: dominant families
WEAK_PUR   <- 0.35; WEAK_MRG <- 0.10
CONF_ENR   <- 3.0;  CONF_ENR_PUR <- 0.25   # route 2: minor families, enrichment over baseline
WEAK_ENR   <- 2.0;  WEAK_ENR_PUR <- 0.18
NAIVE_MULT <- 2.0    # naive/mantle ratio above this multiple of the core median
FDC_MULT   <- 3.0    # CR2/FCER2 ratio above this -> FDC, not B (tested BEFORE naive)
FDC_MAX_PUR<- 0.60   # ...but NEVER override a cluster whose argmax is strongly B
FDC_BCORE  <- 0.80   # ...and require B_core to be depressed against the core median
PDC_MULT   <- 5.0    # pDC module above this multiple of the object median -> pDC
HODG_MULT  <- 5.0    # Hodgkin-chemokine module above this -> a NICHE, not a cell type
CORE_TOL   <- 0.25   # "monomorphic core" = within this of the median on every module
FAMS <- names(LINEAGE_SETS_COARSE)
NICE <- c(B="B cell", T="T cell", NK="NK cell", myeloid="Myeloid",
          endo="Endothelial", fibro="Stromal/FDC", PC="Plasma cell")

cat("============================================================\n")
cat("step 19d -- cluster labels by argmax composition (", CL, ")\n", sep="")
cat("============================================================\n")

hc <- fread(HC_F); ct <- fread(CT_F, select = c("cell_id","celltype","lineage16"))
d  <- ct[hc, on = .(cell_id)]
setnames(d, CL, "cl"); d[, cl := as.character(cl)]
d[, resolved := lineage16 %in% FAMS]

comp <- d[resolved == TRUE, .N, by = .(cl, lineage16)]
comp[, share := N/sum(N), by = cl]
lab <- comp[order(-share), .(top1 = lineage16[1], purity = round(share[1],3),
                             top2 = lineage16[2], margin = round(share[1]-share[2],3)), by = cl]
tot <- d[, .(cells = .N, n_resolved = sum(resolved),
             frac_resolved = round(mean(resolved),3),
             pct_IB1 = round(100*mean(merge_sample == "IB1"),1)), by = cl]
r <- lab[tot, on = .(cl)]

## baseline share of each family among ALL resolved cells in the object
base <- d[resolved == TRUE, .N, by = lineage16][, .(lineage16, base = N/sum(N))]
cat("\n=== family baselines among resolved cells ===\n"); print(base[order(-base)])
r <- base[r, on = c(lineage16 = "top1")]
setnames(r, "lineage16", "top1")
r[, enrichment := round(purity / base, 2)]
r[, status := fifelse(frac_resolved < JUNK_FRAC, "junk",
              fifelse((purity >= CONF_PUR & margin >= CONF_MRG) |
                      (enrichment >= CONF_ENR & purity >= CONF_ENR_PUR), "confident",
              fifelse((purity >= WEAK_PUR & margin >= WEAK_MRG) |
                      (enrichment >= WEAK_ENR & purity >= WEAK_ENR_PUR), "weak_support", "mixed")))]
r[, n_fam_genes := lengths(LINEAGE_SETS_COARSE)[top1]]
## MIN_FAM_GENES, carried from 09b: a family with fewer than two panel genes can never be
## confident. All seven coarse families clear it here (4-7 genes each), but the guard stays
## because the Neutrophil family -- one gene on this panel, CSF3R, with FCGR3B absent --
## is the standing reason it exists.
r[status == "confident" & n_fam_genes < 2L, status := "weak_support"]

## ---- B second stage ---------------------------------------------------------
r[, label := NA_character_]
r[status %in% c("confident","weak_support"), label := NICE[top1]]
b_note <- NULL
if (file.exists(MOD_F)) {
  m <- fread(MOD_F); m[, cl := as.character(cl)]
  w <- dcast(m[top == "B"], cl ~ module, value.var = "cp10k")
  mods <- setdiff(names(w), c("cl","B_core"))
  for (k in mods) w[[k]] <- w[[k]] / w[["B_core"]]
  med <- sapply(mods, function(k) median(w[[k]]))
  core <- w[, Reduce(`&`, lapply(mods, function(k) abs(get(k)/med[k]-1) < CORE_TOL))]
  cat(sprintf("\nmonomorphic B core: %d of %d B-dominant clusters\n", sum(core), nrow(w)))
  ref <- sapply(mods, function(k) median(w[[k]][core]))
  w[, naive_ratio := naive_mantle / ref[["naive_mantle"]]]
  w[, fdc_ratio := if ("fdc_shared" %in% names(w)) fdc_shared / ref[["fdc_shared"]] else 0]
  w[, bcore_ratio := B_core / median(w$B_core[core])]
  w <- r[, .(cl, .pur = purity)][w, on = .(cl)]
  ## FDC is tested FIRST. CR2/FCER2 are shared between naive B and follicular dendritic
  ## cells, so a cluster high on them is only naive B if its B_core holds up; an FDC has
  ## CR2 high and MS4A1/CD79A depressed. Getting this order wrong is what put an FDC
  ## cluster under a B label in the first pass.
  ## The FDC override needs all three conditions. CR2/FCER2 alone is not enough: the ratio
  ## is taken over B_core, so ANY cluster with a thin B signal inflates it automatically.
  ## On the first pass that relabelled a 102,783-cell cluster with argmax B purity 0.806 as
  ## FDC. A real FDC cluster has CR2 high AND B_core depressed AND argmax not strongly B --
  ## per-cell marker evidence is the more trustworthy signal and must not be overridden by
  ## a module ratio.
  w[, b_sub := fifelse(fdc_ratio > FDC_MULT & bcore_ratio < FDC_BCORE & .pur < FDC_MAX_PUR,
                       "Stromal/FDC",
               fifelse(naive_ratio > NAIVE_MULT, "B cell, naive/mantle",
                       "B cell"))]
  cat("\n=== B second stage: FDC tested before naive ===\n")
  print(w[order(-fdc_ratio), .(cl, fdc_ratio = round(fdc_ratio,2),
                               bcore_ratio = round(bcore_ratio,2), argmax_B_purity = .pur,
                               naive_ratio = round(naive_ratio,2), b_sub)])
  r <- w[, .(cl, naive_ratio = round(naive_ratio,2), fdc_ratio = round(fdc_ratio,2),
             bcore_ratio = round(bcore_ratio,2), b_sub)][r, on = .(cl)]
  r[label == "B cell" & !is.na(b_sub), label := b_sub]
  b_note <- sprintf("B split on naive/mantle > %.1fx the monomorphic-core median; %d clusters merged into one B label",
                    NAIVE_MULT, r[label == "B cell", .N])
} else cat("\nNOTE: no ", MOD_F, " -- B left as a single label\n", sep="")

## ---- supplementary rules: populations the coarse marker sets cannot name ------
## Applied ONLY where the argmax composition left the cluster unnamed, so a module score can
## never override per-cell marker evidence that already succeeded. A labelled cluster that
## also scores high is reported as a warning rather than silently relabelled.
r[, is_niche := FALSE]
if (exists("m") && all(c("pDC","hodgkin_niche") %in% m$module)) {
  sup <- dcast(m[module %in% c("pDC","hodgkin_niche")], cl ~ module, value.var = "cp10k")
  sup[, `:=`(pdc_ratio = pDC / median(pDC), hodg_ratio = hodgkin_niche / median(hodgkin_niche))]
  r <- sup[, .(cl, pdc_ratio = round(pdc_ratio,2), hodg_ratio = round(hodg_ratio,2))][r, on=.(cl)]
  cat("\n=== supplementary module scores (ratio to the object median) ===\n")
  print(r[pdc_ratio > 2 | hodg_ratio > 2,
          .(cl, cells, pct_IB1, top1, purity, status, label, pdc_ratio, hodg_ratio)][order(-hodg_ratio)])
  clash <- r[!is.na(label) & (pdc_ratio > PDC_MULT | hodg_ratio > HODG_MULT)]
  if (nrow(clash)) { cat("\nWARNING -- already-labelled clusters also scoring high; NOT relabelled:\n"); print(clash[, .(cl, label, pdc_ratio, hodg_ratio)]) }
  r[is.na(label) & pdc_ratio  > PDC_MULT,  `:=`(label = "pDC", status = "module_call")]
  r[is.na(label) & hodg_ratio > HODG_MULT,
    `:=`(label = "Hodgkin chemokine niche", status = "module_call", is_niche = TRUE)]
}

setorder(r, -cells)
cat("\n=== cluster labels ===\n")
print(r[, .(cl, cells, pct_IB1, frac_resolved, top1, purity, enrichment, margin, status, label)])

cat("\n=== label totals (high-confidence cells) ===\n")
d2 <- r[, .(cl, label, status)][d, on = .(cl)]
print(d2[, .(cells = .N, pct = round(100*.N/nrow(d2),2)),
         by = .(label = fifelse(is.na(label), paste0("<", status, ">"), label))][order(-cells)])
## halo_suspect travels WITH the data, not only in a log: a caveat that lives in a console
## transcript is a caveat that will be lost the first time someone reads the csv.
HS <- sprintf("%s/halo_score_%s.csv", OUT, RES)
if (file.exists(HS)) {
  hs <- fread(HS)[, .(label = label, halo_fold = fold, halo_suspect)]
  r  <- hs[r, on = .(label)]
  r[is.na(halo_suspect), halo_suspect := FALSE]
  cat("\n=== halo QC joined from 19e ===\n")
  print(r[!is.na(label), .(cl, label, halo_fold, halo_suspect)][order(-halo_fold)])
} else { r[, `:=`(halo_fold = NA_real_, halo_suspect = NA)] ; cat("\nNOTE: no ", HS, " yet -- run 19e then re-run this\n", sep="") }
sc(sprintf("hc_cluster_labels_%s.csv", RES), r)
d2 <- r[, .(cl, label, status, is_niche, halo_suspect)][d, on = .(cl)]
sc(sprintf("hc_cells_labelled_%s.csv.gz", RES),
   d2[, .(cell_id, merge_sample, cl, celltype_cluster = label, cluster_status = status,
          is_niche, halo_suspect,
          celltype_argmax = celltype, lineage_argmax = lineage16, U1, U2)])
if (!is.null(b_note)) cat("\n", b_note, "\n", sep="")
cat("\n============================================================\n")
cat("cell-type STATEMENTS use cluster_status in {confident, weak_support, module_call}\n")
cat("AND is_niche == FALSE AND halo_suspect == FALSE.\n")
cat("  * `Hodgkin chemokine niche` is a REGION, not a cell type: per-cell HRS calling is a\n")
cat("    closed avenue on this panel (TNFRSF8 maxes at 6 counts, 78% singletons).\n")
cat("  * halo_suspect comes from 19e and currently fires on Plasma cell at 12.7x -- worse\n")
cat("    than the 8.4x step 18 measured for the InSituType PC call. Use celltype_argmax\n")
cat("    for plasma cells, which step 18 showed is 2.7x less neighbourhood-contaminated.\n")
cat("============================================================\n")
cat("done\n")
