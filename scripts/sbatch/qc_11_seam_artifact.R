## Are the step-11 CCL17/CCL22 hotspots an FOV-SEAM artifact?
##
## Raised from the rendered views: the red 100 um bin sits on or beside the FOV seam
## cross in the panels inspected, and hotspot (77,61) had 24 cells against 50-75 in its
## neighbours, with the adjacent bin dropped for <20.
##
## WHY THE EXISTING CHECK MISSED IT. scripts/11_script_hodgkin_niche.R tests whether
## hotspots concentrate in PARTICULAR FOVs (they do not: 25 FOVs, 1.4 bins each, FOV
## identity explains 31% of variance). A seam artifact is present in EVERY FOV, so it
## produces no FOV-identity signal at all. That test cannot see this failure mode.
##
## MECHANISM IT WOULD IMPLY. At a seam, cells are clipped by the FOV edge -> fewer
## transcripts per cell -> a small nCount denominator; and transcripts in the overlap
## can be assigned on both sides -> an inflated numerator. Either inflates a
## share-of-transcriptome score exactly where the hotspots appear to sit.
##
## TEST. Per-FOV mm bounding boxes are derived EMPIRICALLY from the cells themselves
## (no assumption about tiling), then each bin's distance to the nearest FOV edge is
## compared between hotspots and background, in both samples.
##
##   sbatch --mem=32G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##     scripts/sbatch/qc_11_seam_artifact.R
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); setDTthreads(1)
})
source("scripts/lib/cosmx_tiff.R")

SEG_F <- "output/10_hrs/hrs_segments.csv.gz"
BIN_F <- "output/11_hodgkin_niche/niche_bins_all_scales.csv.gz"
OUT   <- "output/11_hodgkin_niche"; FIG <- file.path(OUT, "figures")
BINUM <- 100; NICHE_Q <- 0.99

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
COLS <- c(hotspot = "#2a78d6", background = "#a9a7a1")
theme_viz <- function(base=11) theme_minimal(base_size=base) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour=GRID,linewidth=.25),
        plot.title=element_text(colour=INK,face="bold",size=base+2),
        plot.subtitle=element_text(colour=INK2,size=base-1),
        plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
        strip.text=element_text(colour=INK2,face="bold"),
        legend.position="top", legend.justification="left")
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF)
                       cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

cat("============================================================\n")
cat("step 11 QC -- are the hotspots an FOV-seam artifact?\n")
cat("============================================================\n")

SEG <- fread(SEG_F, select = c("cell_id","samp","nCount_RNA","x_slide_mm","y_slide_mm"))
SEG[, fov := cosmx_fov_from_id(cell_id)]
BIN <- fread(BIN_F)[bin_um == BINUM]

## ---- 1. empirical mm bounding box per FOV -----------------------------------
## Derived from the cells, so no assumption about how FOVs tile or which way the axes
## run. The 2nd/98th percentiles guard against a few stray mis-assigned cells.
FB <- SEG[, .(x1 = quantile(x_slide_mm, .02), x2 = quantile(x_slide_mm, .98),
              y1 = quantile(y_slide_mm, .02), y2 = quantile(y_slide_mm, .98),
              n = .N), by = .(samp, fov)][n >= 50]
cat(sprintf("\nFOVs with >=50 cells: %s\n",
            paste(FB[, .N, by=samp][, paste0(samp,"=",N)], collapse=" ")))
cat(sprintf("median FOV mm extent: %.3f x %.3f (expected ~%.3f)\n",
            median(FB$x2-FB$x1), median(FB$y2-FB$y1), 4256*COSMX_MM_PER_PX))

## ---- 2. distance from each bin centre to the nearest FOV edge ---------------
seam_dist <- function(S) {
  b  <- BIN[samp == S]; f <- FB[samp == S]
  ## nearest FOV by centre, then distance to that FOV's four edges
  fc <- as.matrix(f[, .((x1+x2)/2, (y1+y2)/2)])
  nn <- RANN::nn2(fc, as.matrix(b[, .(x_mm, y_mm)]), k = 1)$nn.idx[,1]
  ff <- f[nn]
  d  <- pmin(abs(b$x_mm - ff$x1), abs(b$x_mm - ff$x2),
             abs(b$y_mm - ff$y1), abs(b$y_mm - ff$y2)) * 1000
  b[, seam_um := d][]
}
B <- rbindlist(lapply(unique(BIN$samp), seam_dist))
B[, hot := samp == "IB1" & chemo_cp10k >= quantile(BIN[samp=="IB1", chemo_cp10k], NICHE_Q)]
B[, grp := fifelse(hot, "hotspot", "background")]

cat("\n=== distance to the nearest FOV edge, hotspot vs background\n")
tab <- B[, .(bins = .N, med_seam_um = round(median(seam_um), 1),
             pct_within_50um = round(100*mean(seam_um <= 50), 1),
             med_cells = as.numeric(median(cells)),
             med_chemo = round(median(chemo_cp10k), 1)), by = .(samp, grp)][order(samp, grp)]
print(tab); sc("niche_seam_check.csv", tab)
pw <- suppressWarnings(tryCatch(
  wilcox.test(B[grp=="hotspot", seam_um], B[samp=="IB1" & grp=="background", seam_um])$p.value,
  error = function(e) NA_real_))
cat(sprintf("Wilcoxon hotspot vs IB1 background: p = %s\n", signif(pw, 3)))

## ---- 3. does the score itself track seam proximity? --------------------------
ib <- B[samp == "IB1"]
ib[, seam_bin := cut(seam_um, c(-1, 25, 50, 100, 150, 200, Inf),
                     labels = c("<25","25-50","50-100","100-150","150-200",">200"))]
prof <- ib[, .(bins = .N, med_chemo = round(median(chemo_cp10k), 1),
               mean_chemo = round(mean(chemo_cp10k), 1),
               med_cells = as.numeric(median(cells)),
               med_nCount = as.numeric(median(nCount)),
               pct_hot = round(100*mean(hot), 2)), by = seam_bin][order(seam_bin)]
cat("\n=== IB1 bins profiled by distance to the nearest FOV edge\n")
cat("    A seam artifact -> chemokine rises and cells/nCount fall as distance -> 0\n")
print(prof); sc("niche_seam_profile.csv", prof)

## ---- 4. the decisive version: does the result survive dropping seam bins? ---
## If hotspots are seam artifacts, excluding seam-proximal bins should collapse the
## IB1-vs-RB1 contrast. If the contrast holds on interior bins only, it is real.
for (cutoff in c(0, 25, 50, 75)) {
  bb <- B[seam_um > cutoff]
  ibb <- bb[samp == "IB1"]; rbb <- bb[samp == "RB1"]
  if (!nrow(ibb) || !nrow(rbb)) next
  cut99 <- quantile(ibb$chemo_cp10k, NICHE_Q)
  cat(sprintf("  drop bins within %3d um of a seam -> IB1 %5d bins (med %5.1f), RB1 %5d (med %4.1f) | IB1 q99 %6.1f | RB1 bins above it: %d\n",
              cutoff, nrow(ibb), median(ibb$chemo_cp10k), nrow(rbb),
              median(rbb$chemo_cp10k), cut99, sum(rbb$chemo_cp10k >= cut99)))
}

p1 <- ggplot(ib, aes(seam_um, chemo_cp10k + .1)) +
  geom_point(aes(colour = grp), size = .5, stroke = 0, alpha = .5) +
  scale_colour_manual(values = COLS, name = NULL,
                      guide = guide_legend(override.aes = list(size = 3))) +
  scale_y_log10() +
  labs(title = "Do the hotspots sit on FOV seams?",
       subtitle = paste0("Each point is a 100 um IB1 bin. A seam artifact would push the ",
                         "blue points to the left edge.\nThe existing FOV check could not ",
                         "see this: seams occur in every FOV, so they leave no\n",
                         "FOV-identity signal."),
       x = "distance from bin centre to the nearest FOV edge (um)",
       y = "CCL17+CCL22 cp10k (log, +0.1)",
       caption = sprintf("Data: %s/niche_seam_check.csv", OUT)) + theme_viz()
sf("fig_11_seam_check.png", p1, 8.5, 5.5)

cat("\n============================================================\n")
h <- tab[samp=="IB1" & grp=="hotspot"]; g <- tab[samp=="IB1" & grp=="background"]
if (nrow(h) && nrow(g)) {
  if (h$med_seam_um < 0.6 * g$med_seam_um)
    cat("VERDICT: hotspots ARE seam-proximal. Treat the call as contaminated and re-run\n",
        "         step 11 excluding seam-adjacent bins.\n", sep = "")
  else
    cat("VERDICT: hotspots are NOT closer to seams than background. Check section 4 --\n",
        "         if the IB1:RB1 contrast also survives dropping seam bins, the call stands.\n", sep="")
}
cat("============================================================\nOK\n")
