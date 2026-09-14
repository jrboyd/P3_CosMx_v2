## Is the IB1 CCL17/CCL22 elevation a SELECTION EFFECT of the spatial QC filter?
##
## THE CONCERN. Step 11 runs on post-QC cells, which is 15.8% of IB1 (400,111 of
## 2,540,021) but 74.9% of RB1. The dominant filter is IB1-only: `qc_spatial_staining`
## marks 2,013,225 IB1 cells "dim_membrane" against 526,796 "ok", while every RB1 cell
## is "not_evaluated" -- the test was built for IB1's over-fixation defect and never
## run on RB1. And 613 of 758 IB1 FOVs have ZERO passing cells, so most of the imaged
## tissue contributes nothing.
##
## If the surviving 16% is better-preserved tissue that also detects more transcripts,
## the chemokine elevation could be a filter artifact rather than biology.
##
## THIS TEST IS WITHIN IB1 ONLY, so it does not depend on the RB1 null at all:
##   1. CCL17+CCL22 per cell, staining "ok" vs "dim_membrane", depth-standardised
##   2. the same split by FOV-edge status, since edge cells are the seam question
##   3. bins recomputed on ALL pre-QC cells: how much coverage is gained, do the 36
##      hotspots survive, and do new ones appear in previously-dark FOVs
##
## Deliberately NOT attempted: recomputing the IB1-vs-RB1 contrast pre-QC. RB1 has no
## staining-failed population at all, so pre-QC IB1 would gain 2.0M cells with no RB1
## counterpart and the comparison would be meaningless.
##
##   sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 \
##     scripts/sbatch/run_rscript.sh scripts/sbatch/qc_11_preqc_sensitivity.R
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); setDTthreads(1)
})

SAMPLE <- "IB1"
IN_RDS <- sprintf("output/01_import_v1/seuratObject_%s_BPCells_full.Rds", SAMPLE)
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
QC_CSV <- file.path(.qc_dir, "qc_spatial_cells.csv.gz")
BIN_F  <- "output/11_hodgkin_niche/niche_bins_all_scales.csv.gz"
OUT    <- "output/11_hodgkin_niche"; FIG <- file.path(OUT, "figures")
CHEMO  <- c("CCL17","CCL22"); EXTRA <- c("TNFRSF8","CD3E","MS4A1")
BIN_UM <- 100; MIN_CELLS <- 20L; DEFECT_CL <- "i"; NICHE_Q <- 0.99

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
GCOLS <- c(ok = "#2a78d6", dim_membrane = "#eb6834")
theme_viz <- function(base=11) theme_minimal(base_size=base) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour=GRID,linewidth=.25),
        plot.title=element_text(colour=INK,face="bold",size=base+2),
        plot.subtitle=element_text(colour=INK2,size=base-1),
        plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
        legend.position="top", legend.justification="left")
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF)
                       cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

cat("============================================================\n")
cat("step 11 QC -- is the chemokine signal a spatial-QC selection effect?\n")
cat("============================================================\n")

obj <- readRDS(IN_RDS); rn <- rownames(obj@meta.data)
NEED <- c("Run_Tissue_name","x_slide_mm","y_slide_mm","nCount_RNA","inSitu_clusters")
md <- as.data.table(obj@meta.data[, NEED, drop = FALSE]); md[, cell_id := rn]
cnt <- LayerData(obj[["RNA"]], layer = "counts")
G <- intersect(c(CHEMO, EXTRA), rownames(cnt))
M <- as(cnt[G, , drop = FALSE], "dgCMatrix"); colnames(M) <- rn
rm(obj, cnt); invisible(gc())

md <- md[Run_Tissue_name == SAMPLE][!is.na(x_slide_mm) & !is.na(y_slide_mm) & nCount_RNA > 0]
qc <- fread(QC_CSV, select = c("cell","sample_id","qc_spatial_staining",
                               "qc_spatial_fovpx","qc_passing_all"))
md <- merge(md, qc[sample_id == SAMPLE, .(cell_id = cell, qc_spatial_staining,
                                          qc_spatial_fovpx, qc_passing_all)],
            by = "cell_id", all.x = TRUE, sort = FALSE)
gi <- match(md$cell_id, colnames(M)); stopifnot(!anyNA(gi))
md[, chemo := as.vector(Matrix::colSums(M[intersect(CHEMO, G), gi, drop = FALSE]))]
md[, chemo_cp10k := 1e4 * chemo / nCount_RNA]
cat(sprintf("\n%s pre-QC cells with coords/counts: %d\n", SAMPLE, nrow(md)))
print(md[, .N, by = qc_spatial_staining])

## ---- 1. staining pass vs fail, depth-standardised ---------------------------
md[, dec := cut(nCount_RNA, unique(quantile(nCount_RNA, 0:10/10)),
                include.lowest = TRUE, labels = FALSE)]
w <- md[, .(w = as.numeric(.N)), by = dec][, w := w/sum(w)]
std_by <- function(d, col) {
  s <- merge(d[, .(v = mean(chemo_cp10k), n = .N), by = c("dec", col)], w, by = "dec")
  s[, .(cp10k_std = round(sum(v*w)/sum(w), 2), cells = sum(n)), by = c(col)]
}
cat("\n=== CCL17+CCL22 per cell by staining flag, depth-standardised\n")
cat("    A selection effect -> 'ok' cells far higher than 'dim_membrane'\n")
s1 <- std_by(md, "qc_spatial_staining")
s1 <- merge(s1, md[, .(raw_mean = round(mean(chemo_cp10k), 2),
                       med_depth = as.numeric(median(nCount_RNA))),
                  by = qc_spatial_staining], by = "qc_spatial_staining")
print(s1); sc("preqc_staining_split.csv", s1)
r <- s1[qc_spatial_staining == "ok", cp10k_std] /
     s1[qc_spatial_staining == "dim_membrane", cp10k_std]
cat(sprintf("ok : dim_membrane ratio = %.2fx\n", r))

## ---- 2. FOV-edge cells, the seam question ------------------------------------
md[, edge := fifelse(qc_spatial_fovpx == "ok", "interior", "edge_or_excluded")]
cat("\n=== the same split by FOV-edge status\n")
s2 <- std_by(md, "edge")
s2 <- merge(s2, md[, .(med_depth = as.numeric(median(nCount_RNA))), by = edge], by = "edge")
print(s2); sc("preqc_edge_split.csv", s2)

## ---- 3. bins on ALL pre-QC cells --------------------------------------------
mk_bins <- function(d, tag) {
  d <- d[inSitu_clusters != DEFECT_CL]
  b <- d[, .(cells = .N, nCount = sum(nCount_RNA), chemo = sum(chemo),
             fov_n = uniqueN(cell_id)),
         by = .(bx = floor(x_slide_mm*1000/BIN_UM), by_ = floor(y_slide_mm*1000/BIN_UM))]
  b <- b[cells >= MIN_CELLS][, chemo_cp10k := 1e4*chemo/nCount][, set := tag][]
}
pre  <- mk_bins(md, "pre-QC (all cells)")
post <- mk_bins(md[qc_passing_all %in% TRUE], "post-QC (current)")
cat("\n=== bin coverage, pre-QC vs the current post-QC analysis\n")
print(rbind(pre, post)[, .(bins = .N, med_cells = as.numeric(median(cells)),
                           med_cp10k = round(median(chemo_cp10k), 1),
                           q99_cp10k = round(quantile(chemo_cp10k, .99), 1),
                           max_cp10k = round(max(chemo_cp10k))), by = set])
sc("preqc_bin_coverage.csv", rbind(pre, post))

## do the current hotspots survive, and where are the new ones?
hot_post <- post[chemo_cp10k >= quantile(post$chemo_cp10k, NICHE_Q)]
hot_pre  <- pre[ chemo_cp10k >= quantile(pre$chemo_cp10k,  NICHE_Q)]
k <- function(d) paste(d$bx, d$by_)
cat(sprintf("\nhotspots post-QC %d | pre-QC %d | overlapping bins %d\n",
            nrow(hot_post), nrow(hot_pre), sum(k(hot_pre) %in% k(hot_post))))
cat(sprintf("post-QC hotspot bins that still exist pre-QC at all: %d of %d\n",
            sum(k(hot_post) %in% k(pre)), nrow(hot_post)))
newbins <- pre[!k(pre) %in% k(post)]
cat(sprintf("bins that exist ONLY pre-QC: %d (median %.1f cp10k) | of the pre-QC hotspots, %d are in that new territory\n",
            nrow(newbins), median(newbins$chemo_cp10k),
            sum(k(hot_pre) %in% k(newbins))))
sc("preqc_hotspot_overlap.csv",
   data.table(hot_post = nrow(hot_post), hot_pre = nrow(hot_pre),
              overlap = sum(k(hot_pre) %in% k(hot_post)),
              post_hot_present_pre = sum(k(hot_post) %in% k(pre)),
              bins_only_pre = nrow(newbins),
              pre_hot_in_new_territory = sum(k(hot_pre) %in% k(newbins))))

p1 <- ggplot(md[!is.na(qc_spatial_staining)], aes(chemo_cp10k + .1,
                                                  fill = qc_spatial_staining)) +
  geom_density(alpha = .45, colour = NA) + scale_x_log10() +
  scale_fill_manual(values = GCOLS, name = NULL) +
  labs(title = "Is the chemokine signal an artifact of the spatial QC filter?",
       subtitle = paste0("Per cell, IB1 only. The staining filter removes 2.0M of 2.5M ",
                         "IB1 cells and was never run on RB1.\nIf 'ok' sat far above ",
                         "'dim_membrane', the step-11 result would be a selection effect."),
       x = "CCL17+CCL22 cp10k (log, +0.1)", y = "density",
       caption = sprintf("Data: %s/preqc_staining_split.csv", OUT)) + theme_viz()
sf("fig_11_preqc_staining_split.png", p1, 8.5, 5)

cat("\n============================================================\n")
if (is.finite(r) && r > 1.5)
  cat(sprintf("VERDICT: 'ok' cells carry %.2fx the chemokine of 'dim_membrane'. The QC\n         filter IS selecting chemokine-rich tissue -- the step-11 magnitude is\n         inflated by selection and must be re-stated on a common footing.\n", r))
if (is.finite(r) && r <= 1.5 && r >= 0.67)
  cat(sprintf("VERDICT: 'ok' vs 'dim_membrane' ratio is %.2fx -- the filter is not\n         selecting the signal. The step-11 estimate is not a selection artifact,\n         and expanding coverage pre-QC is defensible.\n", r))
if (is.finite(r) && r < 0.67)
  cat(sprintf("VERDICT: 'dim_membrane' cells carry MORE chemokine (%.2fx). The current\n         analysis is UNDER-estimating, and hotspots are being missed in the 613\n         FOVs with no passing cells.\n", r))
cat("============================================================\nOK\n")
