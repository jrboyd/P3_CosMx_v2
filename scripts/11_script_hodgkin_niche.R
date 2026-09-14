## ============================================================================
## 11_script_hodgkin_niche.R -- where is the Hodgkin niche in IB1?
##
## WHY THIS AND NOT A PER-CELL HRS CALL. Step 10 established that per-cell HRS
## identification fails on this panel: TNFRSF8 (CD30) tops out at 6 counts with 78% of
## positives a single transcript, its detection is predicted by sequencing depth at AUC
## 0.716, and the rendered top-scoring cells are ordinary lymphocytes indistinguishable
## from RB1 controls.
##
## But the step-10 follow-up (output/10_hrs/chl_program_per_gene.csv) found real
## Hodgkin-specific signal that the per-cell score had buried. On a common
## depth-standardised counts-per-10k scale, IB1 vs RB1:
##       CCL17  4.37x  (8.56 vs 1.96 -- the highest-expressed gene of the whole set)
##       CCL22  2.71x
##       BATF3  2.60x
##       CD40   1.60x
## while the other 17 scored genes sat at 0.83-1.19. CCL17/TARC and CCL22/MDC are the
## canonical HRS chemokines (serum TARC is a clinical cHL biomarker). Two choices in
## step 10 hid this: equal-weighting 21 genes, and z-scoring WITHIN sample, which makes
## a between-sample comparison impossible by construction.
##
## Those chemokines are SECRETED, so they cannot localise an HRS cell -- but they can
## localise a REGION. That is the question this step asks, and it is the one the data
## supports. The unit of analysis is a spatial bin, not a cell.
##
## RB1 -- the FL recurrence with no Hodgkin component -- is the null throughout: it goes
## through the identical pipeline and its bin-score distribution defines the threshold.
##
## Usage: Rscript scripts/11_script_hodgkin_niche.R
##        scripts/11_submit_hodgkin_niche.sh
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(RANN)
  library(ggplot2); library(data.table); setDTthreads(1)
})

args    <- commandArgs(trailingOnly = TRUE)
QC_MODE <- if (length(args) >= 1 && args[1] %in% c("postqc","preqc")) args[1] else "postqc"
PREQC   <- QC_MODE == "preqc"
IN_RDS  <- "output/04_merge/merged_processed.Rds"
PRE_RDS <- "output/01_import_v1/seuratObject_%s_BPCells_full.Rds"
FOCI_F  <- "output/07_lightchain/IB1/lightchain_restricted_foci_persistence.csv"
## The step-02 output directory was renamed 02_filter_QC -> 02_spatial_QC partway
## through this project. Resolve whichever exists rather than hardcoding either, so this
## script works before, during and after that rename.
.qc_dir <- local({
  cand <- c("output/02_spatial_QC", "output/02_filter_QC")
  hit  <- cand[dir.exists(cand)]
  if (!length(hit)) stop("no step-02 QC directory found; looked for: ",
                         paste(cand, collapse=", "), call. = FALSE)
  hit[1]
})
FOVQC_F <- file.path(.qc_dir, "qc_spatial_fov.csv")
OUT     <- if (PREQC) "output/11_hodgkin_niche_preqc" else "output/11_hodgkin_niche"
FIG     <- file.path(OUT, "figures")

SAMPLES   <- c("IB1", "RB1")
SAMPLE_COL<- "merge_sample"
DEFECT_CL <- "i"
BIN_UM    <- c(50, 100, 200)   # multi-scale, like step 07's radius scan
FOCUS_BIN <- 100               # primary
## MIN_CELLS raised from 20. At 20, bins within 25 um of an FOV seam had median 41
## cells against 85 in the interior, and hotspots were 1.6x over-represented there --
## a small-denominator variance effect, not a real chemokine difference (edge cells
## carry 10.56 cp10k vs 10.13 interior, only 1.04x). 40 removes most of the degraded
## bins; the sweep below reports what each choice costs.
MIN_CELLS <- 40L
MIN_SWEEP <- c(20L, 30L, 40L, 60L, 80L)
## SEAM_UM tightened 50 -> 25. At 50 um the flag caught 42-49% of all bins, which is
## mostly geometry: FOVs are ~440 um and bins are 100 um, so about half the tissue is
## near an edge by construction and the flag carried almost no information. The damage
## is concentrated below 25 um (median 41 cells there vs 69 at 25-50 and 85 interior).
SEAM_UM   <- 25      # a bin this close to an FOV edge is flagged, not silently dropped
N_PERM    <- 200               # spatial-coherence permutations
NULL_Q    <- 0.99              # RB1 quantile: defines the tissue-wide FIELD
NICHE_Q   <- 0.99              # IB1's OWN quantile: defines the focal HOTSPOTS

## The niche signal: the two canonical secreted HRS chemokines, and nothing else.
## Deliberately NOT the 21-gene programme -- that is exactly what diluted them.
CHEMO   <- c("CCL17", "CCL22")
## Secondary cHL-associated genes, reported but not part of the niche score.
SECOND  <- c("BATF3", "CD40", "CD274", "PDCD1LG2")
CD30    <- "TNFRSF8"
## Niche composition. cHL niches are T-rich with Tregs and tumour-associated
## macrophages; these are read out, never used to define the niche.
NICHE_G <- c("CD3E", "CD4", "FOXP3", "PDCD1", "CD68", "CD163")
GENES   <- unique(c(CHEMO, SECOND, CD30, NICHE_G))

if (!file.exists(IN_RDS)) stop("missing ", IN_RDS, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#2a78d6"; DIV_LO<-"#eb6834"; ACC<-"#1baf7a"
SAMP_COLS <- c(IB1 = DIV_HI, RB1 = DIV_LO)
theme_viz <- function(base=11, xgrid=FALSE) {
  th <- theme_minimal(base_size=base) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.major.y=element_line(colour=GRID,linewidth=.25),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.line.x=element_line(colour=AXIS,linewidth=.3), axis.ticks=element_blank(),
          axis.text=element_text(colour=MUTED,size=base-2),
          axis.title=element_text(colour=INK2,size=base-1),
          plot.title=element_text(colour=INK,face="bold",size=base+2),
          plot.subtitle=element_text(colour=INK2,size=base-1),
          plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
          strip.text=element_text(colour=INK2,face="bold",size=base-1),
          legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left",
          legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF)
                       cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

cat("============================================================\n")
cat("step 11 -- where is the Hodgkin niche in IB1?\n")
cat("============================================================\n")

## ---- 1. cells and marker counts ---------------------------------------------
## Post-QC reads one merged object (sample column `merge_sample`); pre-QC reads the two
## step-01 per-sample objects, whose sample column is `Run_Tissue_name`. Extract genes
## and metadata per object and drop it before loading the next, so peak memory is one
## object rather than two.
read_one <- function(path, samp_col) {
  cat("reading ", path, "\n", sep = "")
  o  <- readRDS(path); rn <- rownames(o@meta.data)
  need <- c(samp_col, "x_slide_mm", "y_slide_mm", "nCount_RNA", "inSitu_clusters",
            "lineage", "Area.um2")
  miss <- setdiff(need, colnames(o@meta.data))
  if (length(miss)) stop(path, " meta.data is missing: ", paste(miss, collapse=", "),
                         call. = FALSE)
  m <- as.data.table(o@meta.data[, need, drop = FALSE]); m[, cell_id := rn]
  setnames(m, samp_col, "samp_chr")
  cn <- LayerData(o[["RNA"]], layer = "counts")
  g  <- intersect(GENES, rownames(cn))
  if (!all(CHEMO %in% g)) stop("niche chemokines absent: ",
                               paste(setdiff(CHEMO, g), collapse=", "), call. = FALSE)
  mm <- as(cn[g, , drop = FALSE], "dgCMatrix"); colnames(mm) <- rn
  rm(o, cn); invisible(gc())
  list(md = m, M = mm, G = g)
}

if (PREQC) {
  parts <- lapply(SAMPLES, function(S) {
    f <- sprintf(PRE_RDS, S)
    if (!file.exists(f)) stop("missing ", f, call. = FALSE)
    read_one(f, "Run_Tissue_name")
  })
  G  <- Reduce(intersect, lapply(parts, `[[`, "G"))
  md <- rbindlist(lapply(parts, `[[`, "md"))
  M  <- do.call(cbind, lapply(parts, function(z) z$M[G, , drop = FALSE]))
  rm(parts); invisible(gc())
} else {
  one <- read_one(IN_RDS, SAMPLE_COL)
  md <- one$md; M <- one$M; G <- one$G; rm(one); invisible(gc())
}
if (length(G) < length(GENES))
  message("absent and dropped: ", paste(setdiff(GENES, G), collapse = ", "))
setnames(md, "samp_chr", SAMPLE_COL)

md <- md[get(SAMPLE_COL) %in% SAMPLES][!is.na(x_slide_mm) & !is.na(y_slide_mm) & nCount_RNA > 0]
md[, samp := factor(get(SAMPLE_COL), levels = SAMPLES)]
md <- md[inSitu_clusters != DEFECT_CL]
cat(sprintf("cells after dropping cluster \"%s\": %d\n", DEFECT_CL, nrow(md)))
print(md[, .(cells = .N, med_depth = as.numeric(median(nCount_RNA))), by = samp])

gi <- match(md$cell_id, colnames(M)); stopifnot(!anyNA(gi))
M  <- M[, gi, drop = FALSE]
for (g in G) set(md, j = paste0("g_", g), value = as.vector(M[g, ]))
md[, chemo_ct := rowSums(.SD), .SDcols = paste0("g_", CHEMO)]
md[, lin := as.character(lineage)]
## FOV parsed from the cell id (2nd-to-last token of c_<slide>_<fov>_<cell>), which
## works identically in both modes -- the pre-QC and merged objects do not agree on
## whether a `fov` metadata column is present.
md[, fov_id := as.integer(vapply(strsplit(cell_id, "_", fixed = TRUE),
                                 function(z) z[length(z)-1L], ""))]
stopifnot(!anyNA(md$fov_id))

## ---- 2. bin, at several scales ----------------------------------------------
## Square bins in slide mm. A bin score is a SHARE of the bin's transcriptome, so it
## is controlled for both cellularity and sequencing depth by construction -- the two
## things that wrecked the per-cell score in step 10.
bin_at <- function(B, min_cells = MIN_CELLS) {
  d <- copy(md)
  d[, `:=`(bx = floor(x_slide_mm * 1000 / B), by_ = floor(y_slide_mm * 1000 / B))]
  agg <- d[, c(.(cells = .N, nCount = sum(nCount_RNA),
                 chemo = sum(chemo_ct), cd30 = sum(g_TNFRSF8),
                 x_mm = mean(x_slide_mm), y_mm = mean(y_slide_mm),
                 med_area = median(Area.um2),
                 f_T = mean(lin == "T"), f_myeloid = mean(lin == "myeloid"),
                 f_B = mean(lin == "B"), f_none = mean(lin %in% c("all_markers_0",
                                                                  "below_background")),
                 fov = as.numeric(names(sort(table(fov_id), decreasing = TRUE))[1])),
               lapply(.SD, sum)),
           by = .(samp, bx, by_), .SDcols = paste0("g_", G)]
  agg <- agg[cells >= min_cells]
  agg[, chemo_cp10k := 1e4 * chemo / nCount]
  agg[, cd30_cp10k  := 1e4 * cd30  / nCount]
  for (g in setdiff(G, CHEMO))
    agg[, (paste0("cp_", g)) := 1e4 * get(paste0("g_", g)) / nCount]
  agg[, bin_um := B]
  agg[]
}
## Per-FOV mm bounding boxes, empirical from the cells, so no assumption about tiling.
## 2nd/98th percentiles guard against a few mis-assigned cells.
FOVBOX <- md[, .(x1 = quantile(x_slide_mm, .02), x2 = quantile(x_slide_mm, .98),
                 y1 = quantile(y_slide_mm, .02), y2 = quantile(y_slide_mm, .98), n = .N),
             by = .(samp, fov_id)][n >= 50]
cat(sprintf("\nFOV boxes: %s | median extent %.3f x %.3f mm\n",
            paste(FOVBOX[, .N, by=samp][, paste0(samp,"=",N)], collapse=" "),
            median(FOVBOX$x2-FOVBOX$x1), median(FOVBOX$y2-FOVBOX$y1)))
add_seam <- function(b) {
  out <- rbindlist(lapply(unique(b$samp), function(S) {
    bb <- b[samp == S]; f <- FOVBOX[samp == S]
    if (!nrow(f)) { bb[, seam_um := NA_real_]; return(bb[]) }
    nn <- RANN::nn2(as.matrix(f[, .((x1+x2)/2, (y1+y2)/2)]),
                    as.matrix(bb[, .(x_mm, y_mm)]), k = 1)$nn.idx[,1]
    ff <- f[nn]
    bb[, seam_um := 1000 * pmin(abs(x_mm-ff$x1), abs(x_mm-ff$x2),
                                abs(y_mm-ff$y1), abs(y_mm-ff$y2))][]
  }))
  out[, seam_flag := !is.na(seam_um) & seam_um <= SEAM_UM][]
}

## MIN_CELLS sweep, so the choice is visible rather than asserted
swm <- rbindlist(lapply(MIN_SWEEP, function(mc) {
  b <- add_seam(rbindlist(lapply(BIN_UM, function(B) bin_at(B, mc)), fill = TRUE))
  b <- b[bin_um == FOCUS_BIN]
  ib <- b[samp == "IB1"]
  if (!nrow(ib)) return(NULL)
  cutq <- quantile(ib$chemo_cp10k, NICHE_Q)
  data.table(min_cells = mc, IB1_bins = nrow(ib),
             med_cells = as.numeric(median(ib$cells)),
             pct_seam_flagged = round(100*mean(ib$seam_flag), 1),
             hotspots = sum(ib$chemo_cp10k >= cutq),
             hotspots_seam_free = sum(ib$chemo_cp10k >= cutq & !ib$seam_flag),
             RB1_above = sum(b[samp=="RB1", chemo_cp10k] >= cutq))
}))
cat("\n=== MIN_CELLS sweep at the primary bin size\n"); print(swm)
sc("niche_min_cells_sweep.csv", swm)

BINS <- add_seam(rbindlist(lapply(BIN_UM, function(B) bin_at(B, MIN_CELLS)), fill = TRUE))
cat("\n=== bins per scale (>= ", MIN_CELLS, " cells)\n", sep = "")
print(BINS[, .(bins = .N, med_cells = as.numeric(median(cells)),
               med_chemo = round(median(chemo_cp10k), 2),
               q99_chemo = round(quantile(chemo_cp10k, .99), 2)), by = .(bin_um, samp)][
        order(bin_um, samp)])

## ---- 3. RB1 defines the null -------------------------------------------------
## Every RB1 bin above the threshold is a false positive by construction, so the
## enrichment below is a like-for-like specificity measure, not a p-value.
thr_tab <- rbindlist(lapply(BIN_UM, function(B) {
  b   <- BINS[bin_um == B]
  cut <- quantile(b[samp == "RB1", chemo_cp10k], NULL_Q)
  ib  <- b[samp == "IB1"]; rb <- b[samp == "RB1"]
  data.table(bin_um = B, cut = round(cut, 2),
             IB1_bins = nrow(ib), IB1_above = sum(ib$chemo_cp10k >= cut),
             RB1_bins = nrow(rb), RB1_above = sum(rb$chemo_cp10k >= cut),
             pct_IB1 = round(100*mean(ib$chemo_cp10k >= cut), 2),
             pct_RB1 = round(100*mean(rb$chemo_cp10k >= cut), 2),
             enrichment = round(mean(ib$chemo_cp10k >= cut) /
                                max(mean(rb$chemo_cp10k >= cut), 1e-9), 2))
}))
cat(sprintf("\n=== IB1 bins above RB1's %.0fth percentile (RB1 is %g%% by construction)\n",
            100*NULL_Q, 100*(1-NULL_Q)))
print(thr_tab); sc("niche_threshold_by_scale.csv", thr_tab)

B  <- BINS[bin_um == FOCUS_BIN]
## TWO DIFFERENT QUESTIONS, TWO DIFFERENT THRESHOLDS.
##
## RB1's q99 answers "is IB1 elevated at all" -- and it is, so pervasively that 29% of
## IB1 bins clear it. That is a tissue-wide FIELD, not a niche, and thresholding there
## to find a niche was a mistake: IB1's whole distribution is shifted up, so an
## RB1-relative cut just re-reports the global shift and paints most of the slide.
## Confirmed by scripts/sbatch/qc_11_niche_composition_confound.R: the elevation holds
## within EVERY lineage (2.5x myeloid to 10.8x lineage-negative) and strengthens under
## composition standardisation, so it is a real paracrine field rather than RB1 being
## a differently-composed tissue.
##
## HOTSPOTS are the actual answer to "where", and they need IB1's OWN distribution:
## the top NICHE_Q of IB1 bins, with RB1 kept as the specificity read-out rather than
## the threshold. Measured at 100 um: q0.95 -> 177 bins, q0.99 -> 36 bins, and ZERO RB1
## bins clear any of those cuts.
FIELD_CUT <- quantile(B[samp == "RB1", chemo_cp10k], NULL_Q)
B[, field := chemo_cp10k >= FIELD_CUT]
CUT <- quantile(B[samp == "IB1", chemo_cp10k], NICHE_Q)
B[, niche := samp == "IB1" & chemo_cp10k >= CUT]
## Seam-proximal hotspots are FLAGGED, not dropped: a real region can sit beside a
## seam, and silently excluding it would trade one bias for another. Both counts are
## reported and both lists are written.
B[, niche_seam_free := niche & !seam_flag]
NI <- B[samp == "IB1" & niche == TRUE]
cat(sprintf("\nprimary scale %d um\n  FIELD (>= RB1 q%.0f = %.1f): %d of %d IB1 bins (%.0f%%)\n",
            FOCUS_BIN, 100*NULL_Q, FIELD_CUT, B[samp=="IB1", sum(field)],
            B[samp == "IB1", .N], 100*B[samp=="IB1", mean(field)]))
cat(sprintf("  HOTSPOT (>= IB1 q%.1f = %.1f): %d bins (%d seam-free) | RB1 bins clearing the same cut: %d%s\n",
            100*NICHE_Q, CUT, nrow(NI), B[samp=="IB1", sum(niche_seam_free)],
            B[samp == "RB1", sum(chemo_cp10k >= CUT)],
            if (PREQC) "  [NOT COMPARABLE pre-QC -- see header]" else ""))
cat(sprintf("  MIN_CELLS %d | %.1f%% of IB1 bins are within %g um of an FOV seam\n",
            MIN_CELLS, 100*B[samp=="IB1", mean(seam_flag)], SEAM_UM))

## ---- 4. is it a FOV artifact? -------------------------------------------------
## CosMx has per-FOV batch effects, and a FOV is ~512 um -- larger than these bins --
## so a single bad FOV would produce a block of adjacent high bins that looks exactly
## like a niche. This is the first thing that has to be ruled out.
FOV_R2 <- NA_real_
if (B[samp == "IB1", uniqueN(fov)] > 1) {
  ib <- B[samp == "IB1" & !is.na(fov)]
  fit <- lm(chemo_cp10k ~ factor(fov), data = ib)
  r2  <- summary(fit)$r.squared; FOV_R2 <- round(r2, 4)
  nf  <- ib[niche == TRUE, uniqueN(fov)]
  top <- ib[niche == TRUE, .N, by = fov][order(-N)]
  cat(sprintf("\n=== FOV artifact check\n  FOV identity explains %.1f%% of bin-score variance\n",
              100*r2))
  cat(sprintf("  %d niche bins spread over %d FOVs (%.1f bins/FOV; a single-FOV artifact would concentrate)\n",
              nrow(ib[niche == TRUE]), nf, nrow(ib[niche == TRUE])/max(nf,1)))
  print(head(top, 8))
  fov_ok <- r2 < 0.5 && nf >= 5
  ## and against the step-02 spatial QC: only 145 of 758 IB1 FOVs pass
  if (file.exists(FOVQC_F)) {
    q <- fread(FOVQC_F)
    if (all(c("fov","fov_flag") %in% names(q))) {
      q <- unique(q[, .(fov = as.numeric(fov), fov_flag)])
      ib2 <- merge(ib, q, by = "fov", all.x = TRUE)
      cat("\n  niche bins by step-02 FOV QC flag:\n")
      print(ib2[, .(bins = .N, niche = sum(niche),
                    pct = round(100*mean(niche), 2)), by = fov_flag][order(-niche)])
    }
  }
  sc("niche_fov_check.csv", data.table(fov_r2 = round(r2, 4), niche_bins = nrow(ib[niche==TRUE]),
                                       n_fovs = nf, passes = fov_ok))
} else fov_ok <- NA

## ---- 5. is it spatially coherent? --------------------------------------------
## A real niche is contiguous. Permute the bin scores across IB1 bin POSITIONS and
## compare the mean nearest-neighbour distance among top bins -- if the observed set
## is tighter than permuted, the signal has spatial structure rather than being
## scattered high-variance bins.
ib <- B[samp == "IB1"]
mnn <- function(sel) {
  if (sum(sel) < 3) return(NA_real_)
  xy <- as.matrix(ib[sel, .(x_mm, y_mm)])
  nn <- RANN::nn2(xy, k = 2)
  mean(nn$nn.dists[, 2]) * 1000
}
obs_mnn <- mnn(ib$niche)
set.seed(1); k <- sum(ib$niche)
perm <- vapply(seq_len(N_PERM), function(i) {
  s <- rep(FALSE, nrow(ib)); s[sample.int(nrow(ib), k)] <- TRUE; mnn(s)
}, numeric(1))
p_coh <- (1 + sum(perm <= obs_mnn)) / (N_PERM + 1)
cat(sprintf("\n=== spatial coherence\n  observed mean NN distance among %d niche bins: %.0f um\n",
            k, obs_mnn))
cat(sprintf("  permuted: %.0f um (sd %.0f) -- p = %.4f %s\n", mean(perm), sd(perm), p_coh,
            ifelse(p_coh < 0.05, "(clustered)", "(NOT clustered)")))
coh_ok <- is.finite(p_coh) && p_coh < 0.05
sc("niche_spatial_coherence.csv",
   data.table(niche_bins = k, obs_mnn_um = round(obs_mnn, 1),
              perm_mnn_um = round(mean(perm), 1), perm_sd = round(sd(perm), 1), p = p_coh))

## ---- 6. what is in the niche? -------------------------------------------------
## .SDcols is explicit: a bare lapply(.SD, mean) would hit the character/factor columns
## (samp, fov) and return NA with a warning for each.
SC <- c("chemo_cp10k","cd30_cp10k","f_T","f_myeloid","f_B","f_none","med_area","cells")
.grp <- function(d, g) d[, c(.(group = g, bins = .N), lapply(.SD, mean)), .SDcols = SC]
stat <- rbind(.grp(ib[niche == TRUE],  "IB1 niche"),
              .grp(ib[niche == FALSE], "IB1 background"),
              .grp(B[samp == "RB1"],   "RB1 all"))
stat <- stat[, .(group, bins, chemo_cp10k = round(chemo_cp10k, 2),
                 cd30_cp10k = round(cd30_cp10k, 3),
                 f_T = round(100*f_T, 1), f_myeloid = round(100*f_myeloid, 1),
                 f_B = round(100*f_B, 1), f_none = round(100*f_none, 1),
                 med_area = round(med_area, 1), cells = round(cells, 1))]
cat("\n=== what the niche bins contain\n"); print(stat)
sc("niche_composition.csv", stat)

## secondary markers, niche vs background
sec <- rbindlist(lapply(intersect(c(SECOND, NICHE_G), G), function(g) {
  cn <- paste0("cp_", g)
  data.table(gene = g,
             niche = round(ib[niche == TRUE, mean(get(cn))], 3),
             background = round(ib[niche == FALSE, mean(get(cn))], 3),
             RB1 = round(B[samp == "RB1", mean(get(cn))], 3))
}))
sec[, ratio_niche_bg := round(niche / pmax(background, 1e-9), 2)]
setorder(sec, -ratio_niche_bg)
cat("\n=== secondary and composition markers, niche vs IB1 background (cp10k)\n")
print(sec); sc("niche_secondary_markers.csv", sec)

## ---- 7. is the Hodgkin niche distinct from the FL foci? -----------------------
## IB1 is a COMPOSITE node. If both diseases are real and this niche is Hodgkin, it
## should not simply sit on top of step 07's kappa-restricted FL foci.
if (file.exists(FOCI_F)) {
  fo <- fread(FOCI_F)
  if (all(c("x","y") %in% names(fo)) && nrow(fo)) {
    fxy <- as.matrix(fo[, .(x, y)])
    d_ni <- RANN::nn2(fxy, as.matrix(ib[niche == TRUE, .(x_mm, y_mm)]), k = 1)$nn.dists[,1]*1000
    d_bg <- RANN::nn2(fxy, as.matrix(ib[niche == FALSE, .(x_mm, y_mm)]), k = 1)$nn.dists[,1]*1000
    cat(sprintf("\n=== distance to the nearest step-07 FL focus (%d foci)\n", nrow(fo)))
    cat(sprintf("  niche bins      median %.0f um\n  background bins median %.0f um\n",
                median(d_ni), median(d_bg)))
    pw <- suppressWarnings(tryCatch(wilcox.test(d_ni, d_bg)$p.value,
                                    error = function(e) NA_real_))
    cat(sprintf("  Wilcoxon p = %s -- %s\n", signif(pw, 3),
                ifelse(is.finite(pw) && pw < 0.05 && median(d_ni) > median(d_bg),
                       "niche sits FURTHER from FL foci than background (spatially distinct)",
                       ifelse(is.finite(pw) && pw < 0.05,
                              "niche sits CLOSER to FL foci -- check for confounding",
                              "no difference"))))
    sc("niche_vs_fl_foci.csv",
       data.table(n_foci = nrow(fo), med_dist_niche_um = round(median(d_ni), 1),
                  med_dist_background_um = round(median(d_bg), 1), p = pw))
  }
} else message("missing ", FOCI_F, " -- FL-focus comparison skipped")

## ---- 8. outputs and figures ---------------------------------------------------
sc("niche_bins_all_scales.csv.gz", BINS)
setorder(NI, -chemo_cp10k)
sc("niche_bins_called.csv", NI)
cat("\n=== top 10 niche bins\n")
print(head(NI[, .(bx, by_, x_mm = round(x_mm,2), y_mm = round(y_mm,2), cells,
                  chemo_cp10k = round(chemo_cp10k,1), cd30_cp10k = round(cd30_cp10k,2),
                  f_T = round(100*f_T,1), fov)], 10))

p1 <- ggplot(BINS[bin_um == FOCUS_BIN], aes(chemo_cp10k + .1, fill = samp)) +
  geom_density(alpha = .45, colour = NA) +
  geom_vline(xintercept = CUT, colour = INK, linetype = "22", linewidth = .5) +
  scale_x_log10() + scale_fill_manual(values = SAMP_COLS, name = NULL) +
  labs(title = "CCL17+CCL22 per bin: IB1 against the RB1 null",
       subtitle = sprintf(paste0("%d um bins, >= %d cells. Dashed = RB1's %.0fth percentile. ",
                 "RB1 is the FL recurrence with no\nHodgkin component, so every RB1 bin ",
                 "past the line is a false positive."), FOCUS_BIN, MIN_CELLS, 100*NULL_Q),
       x = "CCL17+CCL22, counts per 10k of the bin's transcriptome (log, +0.1)",
       y = "density", caption = sprintf("Data: %s/niche_bins_all_scales.csv.gz", OUT)) +
  theme_viz()
sf("fig_11_niche_score_vs_null.png", p1, 8.5, 5)

## One colour scale, not two: the background is a fixed pale grey and only the niche
## bins are on the gradient. Two scales would need ggnewscale, which nothing else in
## this project uses, and the second scale adds no information here.
p2 <- ggplot(ib, aes(x_mm, y_mm)) +
  geom_point(data = ib[niche == FALSE], colour = "#e4e3dc", size = 1.1, stroke = 0) +
  geom_point(data = ib[niche == TRUE], aes(colour = chemo_cp10k), size = 2.2, stroke = 0) +
  scale_colour_gradient(low = "#7ba3d0", high = "#0b3d7a",
                        name = "CCL17+CCL22\ncp10k") +
  coord_equal() +
  labs(title = "Where the Hodgkin chemokine signal concentrates in IB1",
       subtitle = sprintf(paste0("%d um bins. Blue = the top %.0f%% of IB1 bins (hotspots); ",
                    "ZERO RB1 bins clear the same cut.\nThe pale background is itself ",
                    "elevated over RB1 -- the field is tissue-wide, the hotspots are focal."),
                          FOCUS_BIN, 100*(1-NICHE_Q)),
       x = "x (mm)", y = "y (mm)",
       caption = sprintf("Data: %s/niche_bins_called.csv", OUT)) +
  theme_viz() + theme(panel.grid.major.y = element_blank())
sf("fig_11_niche_map_IB1.png", p2, 9, 9)

p3 <- ggplot(thr_tab, aes(factor(bin_um), enrichment)) +
  geom_col(fill = DIV_HI, width = .6) +
  geom_hline(yintercept = 1, colour = INK, linetype = "22", linewidth = .4) +
  geom_text(aes(label = sprintf("%.1fx\n%d bins", enrichment, IB1_above)),
            vjust = -0.25, size = 3, colour = INK2) +
  scale_y_continuous(expand = expansion(mult = c(0, .22))) +
  labs(title = "Specificity holds across bin scale",
       subtitle = paste0("Fold-enrichment of IB1 bins over RB1 bins above the same ",
                         "threshold. 1.0 = no better than\nthe FL-only control."),
       x = "bin size (um)", y = "IB1:RB1 enrichment",
       caption = sprintf("Data: %s/niche_threshold_by_scale.csv", OUT)) + theme_viz()
sf("fig_11_scale_robustness.png", p3, 7, 5)

sl <- melt(sec, id.vars = "gene", measure.vars = c("niche","background","RB1"),
           variable.name = "group", value.name = "cp10k")
p4 <- ggplot(sl, aes(reorder(gene, cp10k), cp10k, fill = group)) +
  geom_col(position = position_dodge(.8), width = .74) + coord_flip() +
  scale_fill_manual(values = c(niche = DIV_HI, background = MUTED, RB1 = DIV_LO), name = NULL) +
  labs(title = "What else is elevated in the niche bins",
       subtitle = paste0("None of these define the niche -- it is called on CCL17+CCL22 ",
                         "alone -- so they are\nindependent corroboration."),
       x = NULL, y = "counts per 10k",
       caption = sprintf("Data: %s/niche_secondary_markers.csv", OUT)) + theme_viz(xgrid = TRUE)
sf("fig_11_niche_secondary.png", p4, 8, 6)

## ---- verdict -------------------------------------------------------------------
.one <- function(v) if (length(v) == 1) as.numeric(v) else NA_real_
enr <- .one(thr_tab[bin_um == FOCUS_BIN, enrichment])
SUM <- data.table(
  metric = c("bin_um","niche_bins_IB1","IB1_bins_total","RB1_bins_above","enrichment",
             "fov_r2","spatial_coherence_p","niche_chemo_cp10k","background_chemo_cp10k",
             "niche_pct_T","background_pct_T"),
  value = c(FOCUS_BIN, nrow(NI), B[samp=="IB1", .N],
            .one(thr_tab[bin_um==FOCUS_BIN, RB1_above]), enr,
            FOV_R2,
            p_coh, .one(stat[group=="IB1 niche", chemo_cp10k]),
            .one(stat[group=="IB1 background", chemo_cp10k]),
            .one(stat[group=="IB1 niche", f_T]), .one(stat[group=="IB1 background", f_T])))
sc("summary.csv", SUM)
cat("\n============================================================\nSUMMARY\n")
print(SUM)
cat(sprintf("\n[%s] enrichment over the RB1 null: %.2fx at %d um\n",
            ifelse(is.finite(enr) && enr >= 3, "ok", "NO"), enr, FOCUS_BIN))
cat(sprintf("[%s] spatially coherent (p = %.4f)\n",
            ifelse(isTRUE(coh_ok), "ok", "NO"), p_coh))
cat(sprintf("[%s] not a FOV artifact\n", ifelse(isTRUE(fov_ok), "ok",
            ifelse(is.na(fov_ok), "--", "NO"))))
cat("\nThis locates a REGION, not a cell. CCL17/CCL22 are secreted, so an elevated bin\n")
cat("says Hodgkin-associated chemokine signalling is happening there -- it does not\n")
cat("identify which cell is HRS, and step 10 showed this panel cannot.\n")
cat("============================================================\nOK\n")
