## ============================================================================
## 19_script_cluster_health.R -- THE GATE. No cell type is called until this passes.
##
## Dropping `is_defect` removes InSituType cluster `i`. It does NOT remove low-information
## tissue, and the two are not the same set: IB1's pearson_clusters_low cluster 0 audits as
## defect_dominated at 70.1% defect, which means roughly THIRTY PERCENT of that mass is
## non-defect cells that survive this filter and can re-form a low-signal cluster in the new
## embedding. If they do and we label it, we repeat 09c's laundering failure by a different
## route -- unresolvable tissue acquiring a lineage name because it happened to co-cluster.
##
## So this script asks one question: did the clustering come back clean? It answers it from
## EVIDENCE ABOUT THE CELLS, never from cluster membership, so the check is not circular.
##
##   frac_resolved   share of the cluster's cells that the step-16 argmax caller gives a real
##                   lineage. The object-wide baseline is printed; a cluster far below it is
##                   tissue that carries no marker evidence in ANY cell.
##   depth           median nCount_RNA and nFeature_RNA against the object median.
##   negprobe rate   POOLED sum(negprobes)/sum(counts) per cluster -- the honest
##                   signal-to-noise measure. It must be pooled, not a median of per-cell
##                   ratios: nCount_negprobes is nonzero in only 21.8% of cells, so the
##                   per-cell median is 0 and a median-of-ratios is 0 for every cluster and
##                   for the object. The first version of this gate used the median, which
##                   made the term `negrate > baseline` read `0 > 0` -- FALSE always -- and
##                   since it sat inside an AND, the gate COULD NOT FAIL. It reported PASS
##                   vacuously. Pooled object-wide ratio is 0.000497.
##   sample skew     the 09c cluster-1 signature: one sample supplying almost all of a
##                   cluster its own per-sample audit declines.
##
## THE GATE IS TWO CONDITIONS, NOT FOUR. Requiring size AND resolution AND depth AND
## negprobe to all fail makes failure nearly unreachable -- a four-way AND is not a strict
## test, it is a broken one. The gate proper is: is a LARGE cluster resolution-poor? Depth,
## negprobe rate and cell area are reported alongside as supporting evidence for reading
## that result, not as extra hurdles a problem cluster has to clear.
##
##   low_signal   pct > MIN_PCT AND frac_resolved < FRAC_MULT x baseline   <- the gate
##   junk         frac_resolved < JUNK_FRAC at ANY size                    <- drop these
##   unlabelable  argmax purity < MIN_PURITY or margin < MIN_MARGIN        <- a DIFFERENT
##                problem: the cluster has signal but is a mixture, so 19d cannot name it
##
## If the gate fails, do NOT label and do NOT loosen these numbers. Add per-cell
## low-information criteria (negprobe ratio, nFeature floor) -- defined on the cell's own
## evidence, never on its cluster -- and re-run 19a. Never reach for the FOV-level
## `qc_spatial_i_fraction` flag: it is far too aggressive and cuts IB1 to 196,316.
##
## Usage: sbatch --mem=64G --time=2:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/19_script_cluster_health.R [low|high]
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); setDTthreads(2)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/lineage.R")

args <- commandArgs(trailingOnly = TRUE)
RES  <- if (length(args) >= 1) args[1] else "low"
CL   <- paste0("hc_clusters_", RES)
OUT  <- "output/19_celltypes"; FIG <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
HC_F <- file.path(OUT, "hc_cells.csv.gz")
CT_F <- "output/16_celltypes/cell_celltypes.csv.gz"
MD_F <- "output/16_regions/cells_meta_all.csv.gz"
for (f in c(HC_F, CT_F, MD_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}

MIN_PCT    <- 5.0   # a failing cluster must be at least this share of the object
FRAC_MULT  <- 0.5   # ...and resolve at less than this multiple of the object baseline
JUNK_FRAC  <- 0.15  # a cluster resolving below this is junk at any size
MIN_PURITY <- 0.35  # 19d cannot name a cluster whose modal argmax call is weaker than this
MIN_MARGIN <- 0.10  # ...or whose lead over the runner-up is smaller than this
FAMS      <- names(LINEAGE_SETS_COARSE)

cat("============================================================\n")
cat("step 19b -- THE GATE: did the clustering come back clean? (", CL, ")\n", sep="")
cat("============================================================\n")

hc <- fread(HC_F)
if (!CL %in% names(hc)) stop("no column ", CL, " in ", HC_F, call. = FALSE)
ct <- fread(CT_F, select = c("cell_id","celltype","lineage16"))
md <- fread(MD_F, select = c("cell_id","nCount_RNA","nCount_negprobes","Area.um2",
                             "inSitu_clusters","is_defect"))
d <- md[ct[hc, on = .(cell_id)], on = .(cell_id)]
cat(sprintf("cells %d | clusters %d\n", nrow(d), uniqueN(d[[CL]])))

## nFeature_RNA lives in the vendor metadata, pulled out by 18_script_centroids.sh.
nf <- rbindlist(lapply(c("IB1","RB1"), function(S) {
  f <- sprintf("output/18_leakage/%s_centroids.csv.gz", S)
  if (!file.exists(f)) { cat("NOTE: no ", f, " -- nFeature_RNA omitted\n", sep=""); return(NULL) }
  fread(f, select = c("cell_id","nFeature_RNA"))
}), fill = TRUE)
if (nrow(nf)) d <- nf[d, on = .(cell_id)]
if (!"nFeature_RNA" %in% names(d)) d[, nFeature_RNA := NA_real_]

d[, resolved := lineage16 %in% FAMS]
d[, negrate  := nCount_negprobes / pmax(nCount_RNA, 1)]
BASE_FRAC <- d[, mean(resolved)]
BASE_NC   <- d[, median(nCount_RNA)]
BASE_NEG  <- d[, sum(nCount_negprobes)/sum(nCount_RNA)]
cat(sprintf("\nobject baselines: frac_resolved %.3f | median nCount %.0f | pooled negprobe rate %.6f\n",
            BASE_FRAC, BASE_NC, BASE_NEG))
if (d[, sum(is_defect)] > 0) stop("defect cells are present in the high-confidence object", call.=FALSE)

setnames(d, CL, "cl")
r <- d[, .(cells = .N,
           pct = 100*.N/nrow(d),
           pct_IB1 = round(100*mean(merge_sample == "IB1"), 1),
           frac_resolved = round(mean(resolved), 3),
           med_nCount = as.numeric(median(nCount_RNA)),
           med_nFeature = as.numeric(median(nFeature_RNA, na.rm = TRUE)),
           negrate = round(sum(nCount_negprobes)/sum(nCount_RNA), 6),
           med_area = round(median(Area.um2), 1)), by = cl][order(-cells)]
r[, pct := round(pct, 2)]
r[, low_signal := pct > MIN_PCT & frac_resolved < FRAC_MULT*BASE_FRAC]
r[, junk := frac_resolved < JUNK_FRAC]
r[, skew := pct_IB1 > 90 | pct_IB1 < 10]

## The argmax composition of each cluster -- what 19d will label on, previewed here so a
## cluster that resolves but resolves to a MIXTURE is visible before it is named.
comp <- d[resolved == TRUE, .N, by = .(cl, lineage16)]
comp[, share := N/sum(N), by = cl]
top <- comp[order(-share), .(top = lineage16[1], purity = round(share[1],3),
                             runner = lineage16[2], margin = round(share[1]-share[2],3)), by = cl]
r <- top[r, on = .(cl)]
r[, unlabelable := purity < MIN_PURITY | margin < MIN_MARGIN]
setcolorder(r, c("cl","cells","pct","pct_IB1","frac_resolved","med_nCount","med_nFeature",
                 "negrate","med_area","top","purity","runner","margin",
                 "low_signal","junk","unlabelable","skew"))
cat("\n=== per-cluster health and argmax composition ===\n"); print(r)
sc(sprintf("hc_cluster_health_%s.csv", RES), r)

## ---- the verdict ------------------------------------------------------------
bad <- r[low_signal == TRUE]
cat("\n============================================================\n")
if (nrow(bad) == 0) {
  cat("GATE PASSES for ", CL, ": no cluster over ", MIN_PCT, "% of the object resolves below\n", sep="")
  cat(sprintf("%.0f%% of the object baseline (%.3f). No large low-signal mass survived the\n",
              100*FRAC_MULT, BASE_FRAC))
  cat("defect filter.\n")
} else {
  cat("GATE FAILS for ", CL, ". ", nrow(bad), " cluster(s) are surviving low-signal tissue:\n\n", sep="")
  print(bad[, .(cl, cells, pct, frac_resolved, med_nCount, negrate)])
  cat("\nDo not label and do not loosen these numbers. Add per-cell low-information\n")
  cat("criteria (negprobe ratio, nFeature floor) defined on the cell's own evidence,\n")
  cat("never on its cluster, and re-run 19a.\n")
}
if (r[, any(junk)]) {
  cat("\nJUNK clusters (resolve below ", JUNK_FRAC, " at any size) -- drop before labelling:\n", sep="")
  print(r[junk == TRUE, .(cl, cells, pct, frac_resolved, med_nCount, med_nFeature)])
}
if (r[, any(unlabelable & !junk)]) {
  cat("\nUNLABELABLE clusters -- they have signal but are MIXTURES, so 19d cannot name them.\n")
  cat("This is a resolution problem, not a QC problem: try hc_clusters_high.\n")
  print(r[unlabelable == TRUE & junk == FALSE,
          .(cl, cells, pct, frac_resolved, top, purity, runner, margin)])
}
if (r[, any(skew)]) {
  cat("\nNOTE -- sample-skewed clusters (>90% one sample), the 09c cluster-1 signature:\n")
  print(r[skew == TRUE, .(cl, cells, pct, pct_IB1, frac_resolved)])
}
cat("============================================================\n")

## ---- figures ----------------------------------------------------------------
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
th <- theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.major=element_line(colour=GRID,linewidth=.25), panel.grid.minor=element_blank(),
        axis.text=element_text(colour=MUTED), axis.title=element_text(colour=INK2),
        plot.title=element_text(colour=INK,face="bold"), plot.subtitle=element_text(colour=INK2,size=9),
        strip.text=element_text(colour=INK2,face="bold"))
## Junk clusters set the axis if allowed to: 159 cells sit at U1 ~ -18 against a data range
## of -4.4..4.2. Range from the non-junk clusters, drawn with coord_cartesian so nothing is
## dropped -- only zoomed.
jc <- r[junk == TRUE, cl]
rl <- d[!cl %in% jc]
XL <- range(rl$U1); YL <- range(rl$U2)
pad <- 0.04 * c(diff(XL), diff(YL))
XL <- XL + c(-pad[1], pad[1]); YL <- YL + c(-pad[2], pad[2])
## coord_fixed, not coord_cartesian: a UMAP is a distance embedding, so unequal x and y
## scaling misrepresents it and makes panels of different widths incomparable. Limits are
## passed in ascending order -- ggplot SORTS the limits given to coord_fixed, so a reversed
## pair silently does nothing (see CLAUDE.md).
CO <- coord_fixed(ratio = 1, xlim = XL, ylim = YL, expand = FALSE)
set.seed(1L); sub <- d[sample(.N, min(.N, 300000L))]
## `cl` is an integer, so mapping it to colour directly gives a CONTINUOUS blue gradient and
## eleven clusters render as eleven shades of the same blue -- unreadable, and it makes a
## perfectly separated embedding look like one undifferentiated mass. Categorical always.
sub[, clf := factor(cl)]
lab <- sub[, .(U1=median(U1), U2=median(U2)), by=clf]
pal <- grDevices::hcl.colors(uniqueN(sub$clf), "Dark 3")
p1 <- ggplot(sub, aes(U1,U2,colour=clf)) + geom_point(size=.05, alpha=.4, show.legend=FALSE) +
  scale_colour_manual(values=pal) +
  geom_text(data=lab, aes(label=clf), colour=INK, size=3.4, fontface="bold") +
  labs(title=sprintf("high-confidence object -- %s", CL),
       subtitle=sprintf("%d cells, defect dropped before clustering; %d shown", nrow(d), nrow(sub))) + CO + th
p2 <- ggplot(sub, aes(U1,U2,colour=resolved)) + geom_point(size=.05, alpha=.35) +
  scale_colour_manual(values=c(`TRUE`="#184f95", `FALSE`="#d03b3b"),
                      labels=c(`TRUE`="argmax resolves", `FALSE`="no marker call"), name=NULL) +
  guides(colour=guide_legend(override.aes=list(size=2,alpha=1))) +
  labs(title="where marker evidence exists",
       subtitle=sprintf("object baseline %.1f%% resolved -- a cluster red throughout is the gate's concern",
                        100*BASE_FRAC)) + CO + th
LP <- c("B","T","NK","myeloid","endo","fibro","PC")
sub[, lin := factor(fifelse(lineage16 %in% LP, lineage16, "no call"), levels=c(LP,"no call"))]
p3 <- ggplot(sub[order(lin=="no call", decreasing=TRUE)], aes(U1,U2,colour=lin)) +
  geom_point(size=.05, alpha=.4) +
  scale_colour_manual(values=c(B="#184f95", T="#d03b3b", NK="#7b3f9e", myeloid="#e08214",
                               endo="#1b7837", fibro="#8c6d31", PC="#c51b7d",
                               `no call`="#dcdcd6"), name=NULL) +
  guides(colour=guide_legend(override.aes=list(size=2,alpha=1))) +
  labs(title="argmax lineage on the same embedding",
       subtitle="does the manifold carry lineage structure the clustering can find?") + CO + th
sf(sprintf("hc_gate_%s.png", RES), p1 + p2 + p3 + plot_layout(nrow = 1), 20, 6)
cat("\ndone\n")
