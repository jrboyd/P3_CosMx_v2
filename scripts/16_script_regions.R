## ============================================================================
## 16_script_regions.R -- the canonical per-cell region table for BOTH samples.
##
## Step 16 asks whether FL looks different when it sits next to a Hodgkin field. IB1 is
## the composite node (cHL + FL); RB1 is the FL recurrence, same patient, same clone, no
## Hodgkin. Answering it requires REGIONS, because IB1 is 8.5% B cells post-QC and RB1 is
## 39.9%: compare whole samples and every T-cell or stromal difference is composition.
##
## This script is the only new primitive in step 16. Everything downstream is reuse.
##
## TWO AXES, DELIBERATELY KEPT SEPARATE.
##
##   `region` -- DISEASE regions, sample-specific, from the readouts that already exist:
##               IB1 cHL_hotspot (step 11 chemokine bins) / FL_focus (step 07 light-chain
##               neighbourhoods) / buffer / background.
##
##   `arch`   -- ARCHITECTURE, computed IDENTICALLY in both samples from per-bin B-cell
##               density: B_dense / B_mid / B_sparse. This is the axis that makes the
##               cross-sample comparison like-for-like, and it is the ONLY way to
##               regionalise RB1 at all -- RB1's light-chain readout is degenerate
##               (178/178 neighbourhoods kappa-restricted, 1,276,791 of 1,286,223 cells
##               inside one), so clonality cannot partition it. See CLAUDE.md.
##
## WHY WITHIN-SAMPLE QUANTILES FOR `arch`. IB1's per-bin f_B has q90 = 0.19 and RB1's has
## q10 = 0.25 -- the distributions barely overlap. An absolute cut puts essentially all of
## RB1 in one class and all of IB1 in the other, which is not a regional analysis, it is a
## sample label. Quantiles match on RANK within each tissue; the residual difference in
## absolute composition is then handled by the matching in 16c and by composition
## standardisation, not swept up here. The cut values are printed and swept.
##
## CLUSTER `i` IS KEPT FOR GEOMETRY AND EXCLUDED FROM STATISTICS. Step 14 learned this the
## hard way (its header, lines 122-129): dropping the over-fixation cluster before drawing
## the footprints cut 154 neighbourhoods to 68, including 11 of the 17 kappa-restricted
## ones -- the class the whole question is about. Step 07 defined those neighbourhoods WITH
## those cells, so removing them first redefines the object of study. Bins, however, are
## scored on defect-free cells only, exactly as step 11 does, so the hotspot calls here
## reproduce step 11's.
##
## Usage: Rscript scripts/16_script_regions.R [tier]        tier = strict | matched
##        scripts/16_submit_regions.sh
## ============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(RANN); library(ggplot2); library(patchwork)
  library(data.table); setDTthreads(1)
})

if (!file.exists("scripts/lib/cluster_labels.R")) stop("run from the project root", call. = FALSE)

args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "strict"
if (!TIER %in% c("strict","matched","preqc"))
  stop("tier must be 'strict', 'matched' or 'preqc', got '", TIER, "'", call. = FALSE)
## `preqc` applies NO spatial QC at all -- every segmented cell, 4,406,940 of them. It exists
## for one purpose: MAPPING. Step 11 established that pre-QC is defensible for POSITION and not
## for MAGNITUDE (the spatial filter discards cells carrying MORE chemokine, not less, so a
## pre-QC cp10k is not comparable to a post-QC one), and RB1's staining flag is `not_evaluated`
## so the cross-sample control cannot follow either. Every number this tier produces is
## therefore positional; the thresholds it derives are its own and must not be compared across
## tiers. IB1 is ~59% over-fixation cluster `i` pre-QC, which is dropped per cell as always.

SAMPLES <- c("IB1","RB1")
META_F  <- "output/16_regions/cells_meta_all.csv.gz"
CNT_F   <- "output/16_regions/panel_counts_all.Rds"
LC_F    <- "output/07_lightchain/IB1/lightchain_cells_r40um.csv"
NB_F    <- "output/07_lightchain/IB1/lightchain_neighborhoods_r40um.csv"
OUT     <- file.path("output/16_regions", TIER)
FIG     <- file.path(OUT, "figures")
for (f in c(META_F, CNT_F, LC_F, NB_F))
  if (!file.exists(f)) stop("missing ", f,
    "\n  run scripts/16_script_cell_meta.R and scripts/16_script_panel_counts.R first",
    call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

## ---- parameters --------------------------------------------------------------
CHEMO     <- c("CCL17","CCL22")   # the cHL niche score, step 11's choice, unchanged
BIN_UM    <- c(50, 100, 200)
FOCUS_BIN <- 100L
MIN_CELLS <- 40L                  # step 11's value; the seam analysis behind it is there
MIN_SWEEP <- c(20L, 40L, 80L)
SEAM_UM   <- 25
NULL_Q    <- 0.99                 # RB1 q99 -> the tissue-wide FIELD
NICHE_Q   <- 0.99                 # IB1's own q99 -> the focal HOTSPOTS
ARCH_HI   <- 0.75                 # within-sample fb_a quantile for B_dense
ARCH_LO   <- 0.25                 # ...and for B_sparse
ARCH_MIN_ASSIGNED <- 0.25         # bins below this within-sample f_assigned quantile are
                                  # flagged low-information (flagged, never dropped)
ARCH_SWEEP<- list(c(.90,.10), c(.75,.25), c(.66,.33))
R_IN      <- 40                   # um from an FL member: the focus footprint (step 07's r)
R_BUF     <- 100                  # 40-100 um is a boundary zone, excluded from both sides
BIN_FL_FRAC <- 0.50               # a BIN counts as FL_focus when this share of its cells
                                  # sit inside the footprint; swept below
N_PERM    <- 200
SEED      <- 1L

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#2a78d6"; DIV_LO<-"#eb6834"; ACC<-"#1baf7a"
SAMP_COLS <- c(IB1 = DIV_HI, RB1 = DIV_LO)
REG_COLS  <- c(cHL_hotspot = "#0b3d7a", FL_focus = "#eb6834", buffer = "#d9d7cf",
               background = "#f0efe9", FL_follicular = "#eb6834",
               interfollicular = "#7ba3d0", follicular_mid = "#d9d7cf")
ARCH_COLS <- c(B_dense = "#eb6834", B_mid = "#d9d7cf", B_sparse = "#2a78d6")
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
cat("step 16a -- per-cell regions, both samples, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

## ---- 1. cells at this tier ---------------------------------------------------
md <- fread(META_F)
cat(sprintf("metadata: %d cells, both samples\n", nrow(md)))
if (TIER == "preqc") {
  cat("tier 'preqc': no spatial QC applied -- positional only, see the header\n")
} else {
  tcol <- paste0("tier_", TIER)
  if (!tcol %in% names(md)) stop("no column ", tcol, " in ", META_F, call. = FALSE)
  md <- md[get(tcol) == TRUE]
}
md <- md[!is.na(x_slide_mm) & !is.na(y_slide_mm) & nCount_RNA > 0]
md[, samp := factor(samp, levels = SAMPLES)]
cat(sprintf("\ntier '%s': %d cells\n", TIER, nrow(md)))
print(md[, .(cells = .N, defect = sum(is_defect),
             defect_pct = round(100*mean(is_defect),1),
             med_depth = as.numeric(median(nCount_RNA)),
             med_area = round(median(Area.um2),1), fovs = uniqueN(fov_id)), by = samp])

## ---- 2. chemokine counts -----------------------------------------------------
M <- readRDS(CNT_F)
gi <- match(md$cell_id, colnames(M))
if (anyNA(gi)) stop(sum(is.na(gi)), " cells have no column in ", CNT_F, call. = FALSE)
miss <- setdiff(CHEMO, rownames(M))
if (length(miss)) stop("niche chemokines absent from the panel: ",
                       paste(miss, collapse=", "), call. = FALSE)
md[, chemo_ct := as.vector(Matrix::colSums(M[CHEMO, gi, drop = FALSE]))]
rm(M); invisible(gc())
cat(sprintf("\nCCL17+CCL22: detected in %.2f%% of cells, mean %.3f counts\n",
            100*md[, mean(chemo_ct > 0)], md[, mean(chemo_ct)]))

## ---- 3. bins -----------------------------------------------------------------
## Scored on DEFECT-FREE cells only, exactly as step 11 does, so the hotspot calls here
## are the same construction. A bin's score is a share of the bin's transcriptome, which
## controls for cellularity and depth by construction -- the two things that wrecked the
## per-cell score in step 10.
mdc <- md[is_defect == FALSE]
cat(sprintf("\nbinning on %d defect-free cells (%d dropped)\n", nrow(mdc), md[, sum(is_defect)]))

bin_at <- function(B, min_cells = MIN_CELLS) {
  d <- mdc[, .(samp, x_slide_mm, y_slide_mm, nCount_RNA, chemo_ct, Area.um2,
               lineage, fov_id)]
  d[, `:=`(bx = floor(x_slide_mm * 1000 / B), by_ = floor(y_slide_mm * 1000 / B))]
  agg <- d[, .(cells = .N, nCount = sum(nCount_RNA), chemo = sum(chemo_ct),
               x_mm = mean(x_slide_mm), y_mm = mean(y_slide_mm),
               ## as.numeric() is load-bearing: median() of an INTEGER vector returns
               ## an integer for odd-length groups and a double for even-length ones, and
               ## data.table rejects the type flip between groups.
               med_area = as.numeric(median(Area.um2)),
               med_depth = as.numeric(median(nCount_RNA)),
               f_T = mean(lineage == "T"), f_B = mean(lineage == "B"),
               f_myeloid = mean(lineage == "myeloid"), f_PC = mean(lineage == "PC"),
               f_none = mean(lineage %in% c("all_markers_0","below_background")),
               f_amb = mean(lineage == "ambiguous"),
               fov = as.numeric(names(sort(table(fov_id), decreasing = TRUE))[1])),
           by = .(samp, bx, by_)]
  agg <- agg[cells >= min_cells]
  agg[, chemo_cp10k := 1e4 * chemo / nCount]
  ## ASSIGNED-FRACTION composition. The first version of this script called architecture on
  ## raw f_B and its own diagnostic caught the flaw: B_sparse came out lineage-NEGATIVE
  ## dominated in both samples (median f_none 0.21 in RB1, 0.36 in IB1 matched) while f_T
  ## barely moved, so "B-poor" was tracking how much of the bin could be typed at all
  ## rather than what it was made of. Renormalising to the ASSIGNABLE cells removes that
  ## quality confound: fb_a is the B share of cells that got a lineage call, so a bin full
  ## of unassignable cells is no longer automatically "interfollicular".
  agg[, f_assigned := 1 - f_none - f_amb]
  agg[, fb_a := fifelse(f_assigned > 0, f_B / f_assigned, NA_real_)]
  agg[, ft_a := fifelse(f_assigned > 0, f_T / f_assigned, NA_real_)]
  agg[, fm_a := fifelse(f_assigned > 0, f_myeloid / f_assigned, NA_real_)]
  agg[, bin_um := B]
  agg[]
}

## Per-FOV mm bounding boxes, empirical from the cells (no tiling assumption), 2nd/98th
## percentiles against a few mis-assigned cells -- 11_script_hodgkin_niche.R:213.
FOVBOX <- mdc[, .(x1 = quantile(x_slide_mm,.02), x2 = quantile(x_slide_mm,.98),
                  y1 = quantile(y_slide_mm,.02), y2 = quantile(y_slide_mm,.98), n = .N),
              by = .(samp, fov_id)][n >= 50]
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

B <- add_seam(rbindlist(lapply(BIN_UM, function(b) bin_at(b)), fill = TRUE))
cat("\n=== bins per sample per scale\n")
print(dcast(B[, .N, by = .(samp, bin_um)], samp ~ bin_um, value.var = "N"))

## ---- 4. the cHL axis: field and hotspots ------------------------------------
Bf <- B[bin_um == FOCUS_BIN]
if (!nrow(Bf[samp == "RB1"])) stop("no RB1 bins at the focus scale", call. = FALSE)
FIELD_CUT <- quantile(Bf[samp == "RB1", chemo_cp10k], NULL_Q)
CUT       <- quantile(Bf[samp == "IB1", chemo_cp10k], NICHE_Q)
B[, field := chemo_cp10k >= FIELD_CUT]
B[, niche := samp == "IB1" & chemo_cp10k >= CUT]
B[, niche_seam_free := niche & !seam_flag]
## RB1's own top-1% chemokine bins: the negative-control analog of IB1's hotspots. Kept as
## a FLAG, not a region level -- it closes the "hotspots are just the top 1% of anything"
## objection without pretending RB1 has a Hodgkin niche.
RB_CUT <- quantile(Bf[samp == "RB1", chemo_cp10k], NICHE_Q)
B[, chemo_top := samp == "RB1" & chemo_cp10k >= RB_CUT]

cat(sprintf("\n=== cHL chemokine axis at %d um\n", FOCUS_BIN))
cat(sprintf("  FIELD cut (RB1 q%.0f) = %.1f cp10k -> %d of %d IB1 bins (%.0f%%)\n",
            100*NULL_Q, FIELD_CUT, Bf[samp=="IB1", sum(chemo_cp10k >= FIELD_CUT)],
            Bf[samp=="IB1", .N], 100*Bf[samp=="IB1", mean(chemo_cp10k >= FIELD_CUT)]))
cat(sprintf("  HOTSPOT cut (IB1 q%.0f) = %.1f cp10k -> %d IB1 bins, %d seam-free; %d RB1 bins clear it\n",
            100*NICHE_Q, CUT, B[bin_um==FOCUS_BIN & niche, .N],
            B[bin_um==FOCUS_BIN & niche_seam_free, .N],
            Bf[samp=="RB1", sum(chemo_cp10k >= CUT)]))
cat(sprintf("  RB1 own-q%.0f cut = %.1f cp10k -> %d RB1 bins flagged chemo_top\n",
            100*NICHE_Q, RB_CUT, B[bin_um==FOCUS_BIN & chemo_top, .N]))

## ---- 5. the architecture axis -----------------------------------------------
## Within-sample quantiles of per-bin B-cell fraction. Assigned per scale so the
## sensitivity analysis at 50/200 um works without re-deriving cuts.
arch_call <- function(b, hi = ARCH_HI, lo = ARCH_LO) {
  b <- copy(b)
  b[, `:=`(fb_hi = quantile(fb_a, hi, na.rm = TRUE),
           fb_lo = quantile(fb_a, lo, na.rm = TRUE),
           fa_lo = quantile(f_assigned, ARCH_MIN_ASSIGNED, na.rm = TRUE)),
    by = .(samp, bin_um)]
  b[, arch := fifelse(is.na(fb_a), NA_character_,
              fifelse(fb_a >= fb_hi, "B_dense",
              fifelse(fb_a <= fb_lo, "B_sparse", "B_mid")))]
  b[, arch := factor(arch, levels = c("B_dense","B_mid","B_sparse"))]
  ## Quality flag, kept SEPARATE from the architecture call rather than folded into it.
  ## Flagging beats dropping: a low-information bin can still be genuine tissue, and
  ## silently excluding it trades one bias for another (the same argument step 11 makes
  ## about seam-proximal hotspots). Downstream contrasts report both ways.
  b[, bin_low_info := !is.na(f_assigned) & f_assigned <= fa_lo]
  b[]
}
B <- arch_call(B)
cat(sprintf("\n=== architecture axis, within-sample f_B quantiles (hi q%.2f / lo q%.2f)\n",
            ARCH_HI, ARCH_LO))
at <- B[bin_um == FOCUS_BIN, .(bins = .N, cut_hi = round(fb_hi[1],3), cut_lo = round(fb_lo[1],3),
        med_fb_a = round(median(fb_a, na.rm=TRUE),3),
        med_ft_a = round(median(ft_a, na.rm=TRUE),3),
        med_fm_a = round(median(fm_a, na.rm=TRUE),3),
        med_f_assigned = round(median(f_assigned),3),
        low_info = sum(bin_low_info), med_cells = as.numeric(median(cells)),
        med_depth = as.numeric(median(med_depth))), by = .(samp, arch)][order(samp, arch)]
print(at)
sc("arch_summary.csv", at)
## f_T inside B_sparse is the check that "interfollicular" means T-zone rather than merely
## acellular. If f_none dominates B_sparse instead, say so -- it is a damage class, not
## an architectural one, and the plan's RB1 `interfollicular` level would be mislabelled.
## Does B_sparse actually mean T-zone? Compare the T share among ASSIGNED cells against
## B_dense. If ft_a does not rise as fb_a falls, "interfollicular" is the wrong name for
## this class and the label is reported as unsupported rather than used.
for (S in SAMPLES) {
  hi <- B[bin_um == FOCUS_BIN & samp == S & arch == "B_dense"]
  lo <- B[bin_um == FOCUS_BIN & samp == S & arch == "B_sparse"]
  if (!nrow(hi) || !nrow(lo)) next
  rt <- median(lo$ft_a, na.rm=TRUE) / median(hi$ft_a, na.rm=TRUE)
  rm_ <- median(lo$fm_a, na.rm=TRUE) / median(hi$fm_a, na.rm=TRUE)
  cat(sprintf(paste0("  %s B_sparse vs B_dense: T share among assigned %.3f -> %.3f (%.2fx), ",
    "myeloid %.3f -> %.3f (%.2fx)\n    f_assigned %.3f -> %.3f -> %s\n"), S,
    median(hi$ft_a,na.rm=TRUE), median(lo$ft_a,na.rm=TRUE), rt,
    median(hi$fm_a,na.rm=TRUE), median(lo$fm_a,na.rm=TRUE), rm_,
    median(hi$f_assigned), median(lo$f_assigned),
    if (rt >= 1.5) "T-zone enriched: 'interfollicular' is supported"
    else if (rt >= 1.1) "mildly T-enriched: 'interfollicular' is weakly supported"
    else "NOT T-enriched: this is B-poor tissue, not a T-zone -- do NOT call it interfollicular"))
}

## Sweep the arch quantiles, so the choice is visible rather than asserted.
asw <- rbindlist(lapply(ARCH_SWEEP, function(q) {
  bb <- arch_call(B[bin_um == FOCUS_BIN], q[1], q[2])
  bb[, .(hi = q[1], lo = q[2], bins = .N, med_fb_a = round(median(fb_a, na.rm=TRUE),3),
         med_ft_a = round(median(ft_a, na.rm=TRUE),3)), by = .(samp, arch)]
}))
sc("arch_threshold_sweep.csv", asw[order(samp, arch, -hi)])

## Is B_dense spatially coherent, or just a noise threshold? Permute which bin POSITIONS
## carry the label and compare the mean nearest-neighbour distance -- the same test step 11
## used for its hotspots (11_script_hodgkin_niche.R:345). This is the plan's gate: if RB1's
## B_dense fails, RB1 has no resolvable follicular architecture and the whole cross-sample
## design has to be reconsidered rather than pushed through.
set.seed(SEED)
coh <- rbindlist(lapply(SAMPLES, function(S) {
  ib <- B[bin_um == FOCUS_BIN & samp == S]
  mnn <- function(sel) {
    if (sum(sel) < 3) return(NA_real_)
    xy <- as.matrix(ib[sel, .(x_mm, y_mm)])
    mean(RANN::nn2(xy, k = 2)$nn.dists[,2]) * 1000
  }
  sel <- ib$arch == "B_dense"; k <- sum(sel); obs <- mnn(sel)
  perm <- vapply(seq_len(N_PERM), function(i) {
    s <- rep(FALSE, nrow(ib)); s[sample.int(nrow(ib), k)] <- TRUE; mnn(s) }, numeric(1))
  p <- (1 + sum(perm <= obs)) / (N_PERM + 1)
  data.table(samp = S, B_dense_bins = k, obs_mnn_um = round(obs,1),
             perm_mnn_um = round(mean(perm),1), perm_sd = round(sd(perm),1), p = p,
             verdict = ifelse(is.finite(p) & p < 0.05, "clustered", "NOT clustered"))
}))
cat("\n=== is B_dense spatially coherent?\n"); print(coh)
sc("arch_coherence.csv", coh)
ARCH_OK <- coh[, all(verdict == "clustered")]
if (!ARCH_OK)
  cat("\n*** WARNING: B_dense is not spatially coherent in ",
      paste(coh[verdict != "clustered", samp], collapse=", "),
      ".\n*** The architecture axis is a threshold on noise there, not follicles.",
      "\n*** Downstream cross-sample contrasts on `arch` are NOT interpretable.\n", sep="")

## ---- 6. FL footprints (IB1 only) --------------------------------------------
## Members of KAPPA-RESTRICTED, non-artifact neighbourhoods define FL foci. Distance to
## ANY called neighbourhood is carried alongside, so "FL-specific" can be separated from
## "any light-chain cluster" -- that separation is the circularity control.
lc  <- fread(LC_F); nbh <- fread(NB_F)
fl_nb  <- nbh[class == "kappa-restricted (FL-like)" & artifact_suspect == FALSE, nb]
all_nb <- nbh$nb
mem_fl  <- lc[nb %in% fl_nb]
mem_any <- lc[nb %in% all_nb]
cat(sprintf("\nstep 07 r40um: %d called neighbourhoods, %d kappa-restricted non-artifact\n",
            nrow(nbh), length(fl_nb)))
cat(sprintf("  member cells: %d FL-like, %d in any called neighbourhood\n",
            nrow(mem_fl), nrow(mem_any)))
## Members come from the PRE-QC pool (step 07 ran on the step-01 object), so a member may
## not survive this tier. They still anchor the geometry -- that is the point of using the
## pre-QC pool -- but the overlap is reported rather than left implicit.
cat(sprintf("  of the FL-like members, %d (%.1f%%) are present at tier '%s'\n",
            sum(mem_fl$cell_id %in% md$cell_id),
            100*mean(mem_fl$cell_id %in% md$cell_id), TIER))
if (nrow(mem_fl) < 50) stop("too few FL-like member cells to define foci", call. = FALSE)

md[, `:=`(d_fl_um = NA_real_, d_nb_um = NA_real_, d_hot_um = NA_real_)]
ii <- md[, which(samp == "IB1")]
XY <- as.matrix(md[ii, .(x = x_slide_mm*1000, y = y_slide_mm*1000)])
md[ii, d_fl_um := RANN::nn2(as.matrix(mem_fl[,  .(x_slide_mm*1000, y_slide_mm*1000)]),
                            XY, k = 1)$nn.dists[,1]]
md[ii, d_nb_um := RANN::nn2(as.matrix(mem_any[, .(x_slide_mm*1000, y_slide_mm*1000)]),
                            XY, k = 1)$nn.dists[,1]]
HOT <- B[bin_um == FOCUS_BIN & niche == TRUE]
if (nrow(HOT) >= 1)
  md[ii, d_hot_um := RANN::nn2(as.matrix(HOT[, .(x_mm*1000, y_mm*1000)]), XY, k = 1)$nn.dists[,1]]
md[, is_fl_member := cell_id %in% mem_fl$cell_id]

## ---- 7. join bins to cells and label ----------------------------------------
md[, `:=`(bx = floor(x_slide_mm*1000/FOCUS_BIN), by_ = floor(y_slide_mm*1000/FOCUS_BIN))]
BJ <- B[bin_um == FOCUS_BIN, .(samp, bx, by_, bin_cells = cells, bin_depth = med_depth,
                               chemo_cp10k, f_B, f_T, f_myeloid, f_PC, f_none, f_amb,
                               f_assigned, fb_a, ft_a, fm_a, bin_low_info,
                               fov_bin = fov, seam_um, seam_flag, field, niche,
                               niche_seam_free, chemo_top, arch)]
n0 <- nrow(md)
md <- merge(md, BJ, by = c("samp","bx","by_"), all.x = TRUE, sort = FALSE)
if (nrow(md) != n0) stop("bin join changed the row count", call. = FALSE)
## Cells in bins that failed MIN_CELLS have no bin row. They are labelled `unbinned` rather
## than silently folded into background -- at the matched tier this is a large group.
md[, unbinned := is.na(bin_cells)]
cat(sprintf("\ncells with no qualifying %d um bin (< %d cells): %d (%.1f%%)\n",
            FOCUS_BIN, MIN_CELLS, md[, sum(unbinned)], 100*md[, mean(unbinned)]))
print(md[, .(cells=.N, unbinned=sum(unbinned), pct=round(100*mean(unbinned),1)), by=samp])

## DISEASE regions. Hotspot wins ties: a bin can be both chemokine-high and near an FL
## focus, and step 12 showed the two readouts are not spatially related (p = 0.23), so the
## overlap is small -- but the precedence is stated rather than left to row order.
md[, region := NA_character_]
md[samp == "IB1" & unbinned == TRUE,  region := "unbinned"]
md[samp == "IB1" & unbinned == FALSE, region := fifelse(niche == TRUE, "cHL_hotspot",
                                       fifelse(d_fl_um <= R_IN, "FL_focus",
                                       fifelse(d_fl_um <= R_BUF, "buffer", "background")))]
md[samp == "RB1" & unbinned == TRUE,  region := "unbinned"]
md[samp == "RB1" & unbinned == FALSE, region := fifelse(arch == "B_dense", "FL_follicular",
                                       fifelse(arch == "B_sparse", "interfollicular",
                                               "follicular_mid"))]
if (md[, any(is.na(region))]) stop("some cells were left unlabelled", call. = FALSE)

## BIN-LEVEL region label, which is a different thing from the per-cell one and is needed
## because the bin is the pseudobulk unit in 16d.
##
## `region` above is per CELL: it keys off d_fl_um, a per-cell distance. A 100 um bin is
## about the size of an FL focus (step 07 span_mm runs 0.05-0.13), so a bin routinely
## straddles a footprint boundary and its cells carry two different labels. 16d's own check
## caught 111 such bins. Aggregating a mixed bin into either group would quietly average
## focus tissue with its surroundings.
##
## So a bin is labelled by the SHARE of its defect-free cells that sit inside the footprint,
## with a stated threshold and a sweep, and partially-covered bins go to `buffer` where they
## are excluded from both sides rather than forced into one.
bshare <- md[is_defect == FALSE & unbinned == FALSE & samp == "IB1",
             .(f_in_fl = mean(d_fl_um <= R_IN), f_in_buf = mean(d_fl_um <= R_BUF)),
             by = .(samp, bx, by_)]
B <- merge(B, bshare, by = c("samp","bx","by_"), all.x = TRUE)
B[, bin_region := NA_character_]
B[samp == "IB1", bin_region := fifelse(niche == TRUE, "cHL_hotspot",
                              fifelse(!is.na(f_in_fl) & f_in_fl >= BIN_FL_FRAC, "FL_focus",
                              fifelse(is.na(f_in_buf) | f_in_buf > 0, "buffer", "background")))]
B[samp == "RB1", bin_region := fifelse(arch == "B_dense", "FL_follicular",
                              fifelse(arch == "B_sparse", "interfollicular",
                                      "follicular_mid"))]
cat(sprintf("\n=== BIN-level regions (the 16d pseudobulk unit), FL share threshold %.2f\n",
            BIN_FL_FRAC))
brs <- B[bin_um == FOCUS_BIN, .(bins = .N, cells = sum(cells), fovs = uniqueN(fov),
         med_cells = as.numeric(median(cells)), med_depth = as.numeric(median(med_depth))),
         by = .(samp, bin_region)][order(samp, -bins)]
print(brs); sc("bin_region_summary.csv", brs)
cat("\n  IB1 bins by FL coverage share:\n")
print(B[bin_um == FOCUS_BIN & samp == "IB1" & !is.na(f_in_fl),
        .(bins = .N), by = .(coverage = cut(f_in_fl, c(-.01,0,.25,.5,.75,1),
          labels = c("0","0-25%","25-50%","50-75%","75-100%")))][order(coverage)])
bsw <- rbindlist(lapply(c(0.25, 0.5, 0.75), function(th)
  B[bin_um == FOCUS_BIN & samp == "IB1" & niche == FALSE & !is.na(f_in_fl),
    .(threshold = th, fl_bins = sum(f_in_fl >= th),
      fl_cells = sum(cells[f_in_fl >= th]))]))
cat("\n  FL-bin count vs the coverage threshold:\n"); print(bsw)
sc("bin_fl_threshold_sweep.csv", bsw)

## Push the bin label back onto cells so 16d can group by it without recomputing.
md <- merge(md, B[bin_um == FOCUS_BIN, .(samp, bx, by_, bin_region, f_in_fl, f_in_buf)],
            by = c("samp","bx","by_"), all.x = TRUE)
md[unbinned == TRUE, bin_region := "unbinned"]
if (md[, any(is.na(bin_region))]) stop("some cells have no bin_region", call. = FALSE)
## One label per bin, by construction now -- assert it, because this is exactly what 16d
## checks and a failure here is cheaper to find.
mixed <- md[unbinned == FALSE, .(n = uniqueN(bin_region)), by = .(samp, bx, by_)][n > 1]
if (nrow(mixed)) stop(nrow(mixed), " bins still carry >1 bin_region", call. = FALSE)
cat("\n  every bin carries exactly one bin_region\n")

cat("\n=== cells by region\n")
rs <- md[, .(cells = .N, defect = sum(is_defect), fl_members = sum(is_fl_member),
             bins = uniqueN(paste(bx,by_)), fovs = uniqueN(fov_id),
             med_depth = as.numeric(median(nCount_RNA)),
             med_area = round(median(Area.um2),1),
             pct_B = round(100*mean(lineage=="B"),1),
             pct_T = round(100*mean(lineage=="T"),1),
             pct_myeloid = round(100*mean(lineage=="myeloid"),1),
             med_chemo_cp10k = round(median(chemo_cp10k, na.rm=TRUE),1)),
         by = .(samp, region)][order(samp, -cells)]
print(rs)
sc("region_summary.csv", rs)
cat("\n=== cells by architecture (the shared axis)\n")
as_ <- md[unbinned == FALSE, .(cells = .N, bins = uniqueN(paste(bx,by_)),
          med_fb_a = round(median(fb_a, na.rm=TRUE),3),
          pct_B = round(100*mean(lineage=="B"),1),
          pct_T = round(100*mean(lineage=="T"),1),
          pct_unassigned = round(100*mean(lineage %in%
            c("all_markers_0","below_background","ambiguous")),1),
          low_info = round(100*mean(bin_low_info),1),
          med_depth = as.numeric(median(nCount_RNA))), by = .(samp, arch)][order(samp, arch)]
print(as_)
sc("arch_cell_summary.csv", as_)
cat("\n=== region x architecture, IB1 (are the FL foci the B-dense tissue?)\n")
print(dcast(md[samp=="IB1" & unbinned==FALSE, .N, by=.(region, arch)],
            region ~ arch, value.var="N", fill=0L))

## MIN_CELLS sweep on what it costs in coverage, so the choice is visible.
msw <- rbindlist(lapply(MIN_SWEEP, function(mc) {
  b <- bin_at(FOCUS_BIN, mc)
  b[, .(min_cells = mc, bins = .N, cells_in_bins = sum(cells)), by = samp]
}))
cat("\n=== MIN_CELLS sweep at the focus scale\n"); print(dcast(msw, samp ~ min_cells,
     value.var = c("bins","cells_in_bins")))
sc("min_cells_sweep.csv", msw)

## ---- 8. write ----------------------------------------------------------------
setcolorder(md, c("cell_id","samp","fov_id","x_slide_mm","y_slide_mm","region","bin_region",
                  "arch","d_fl_um","d_nb_um","d_hot_um","is_fl_member","is_defect","unbinned"))
fwrite(md, file.path(OUT, "cell_regions.csv.gz"), compress = "gzip")
cat("\nwrote ", file.path(OUT, "cell_regions.csv.gz"), " (", nrow(md), " rows)\n", sep="")
fwrite(B, file.path(OUT, "bins_all_scales.csv.gz"), compress = "gzip")
cat("wrote ", file.path(OUT, "bins_all_scales.csv.gz"), " (", nrow(B), " rows)\n", sep="")
sc("thresholds.csv", data.table(
  tier = TIER, bin_um = FOCUS_BIN, min_cells = MIN_CELLS,
  field_cut_rb1_q99 = round(FIELD_CUT,2), hotspot_cut_ib1_q99 = round(CUT,2),
  rb1_own_q99 = round(RB_CUT,2), arch_hi_q = ARCH_HI, arch_lo_q = ARCH_LO,
  arch_on = "fb_a (B share among lineage-ASSIGNED cells)",
  arch_min_assigned_q = ARCH_MIN_ASSIGNED,
  r_in_um = R_IN, r_buf_um = R_BUF,
  ib1_hotspots = B[bin_um==FOCUS_BIN & niche, .N],
  ib1_hotspots_seam_free = B[bin_um==FOCUS_BIN & niche_seam_free, .N],
  rb1_bins_clearing_ib1_cut = Bf[samp=="RB1", sum(chemo_cp10k >= CUT)],
  arch_coherent = ARCH_OK))

## ---- 9. maps -----------------------------------------------------------------
## One map per sample per axis. Slide mm, y reversed so the orientation matches every
## other figure in the project (the trap is documented at 07_script_lightchain_views.R:576).
set.seed(SEED)
pl <- md[, .SD[sample.int(.N, min(.N, 200000L))], by = samp]
## SLIDE-MM PANELS REVERSE NEITHER AXIS. The fitted affine has BOTH slopes negative
## (ax = ay = -8313.9 px/mm), so slide mm runs opposite to canvas px on both axes. Canvas-px
## panels therefore reverse both (11_script_hodgkin_niche_views.R `.coord_px` + `.coord_win`),
## and mm panels reverse neither -- which is what step 12's locator does
## (12_script_fl_vs_hodgkin_overlay.R:147, plain coord_fixed()). A first version of this
## script applied scale_y_reverse() here, which mirrored every region map vertically against
## the image and against that locator. Asserted below rather than trusted.
mk <- function(d, col, cols, ttl) {
  ggplot(d[order(get(col) %in% c("background","B_mid","follicular_mid","unbinned"),
                 decreasing = TRUE)],
         aes(x_slide_mm, y_slide_mm, colour = get(col))) +
    geom_point(size = .12, stroke = 0) +
    scale_colour_manual(values = c(cols, unbinned = "#efeee8"), name = NULL,
                        na.value = "#efeee8") +
    coord_fixed() +
    ## Legend inside the panel, not above it: patchwork renders the two panels' legends
    ## on one shared row and they overlapped ("unbinned" landing on top of "B_dense").
    guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = 1)) +
    labs(x = "x (mm)", y = "y (mm)", title = ttl) +
    theme_viz() + theme(legend.position = "right", legend.justification = "top")
}
for (S in SAMPLES) {
  d <- pl[samp == S]
  p <- mk(d, "region", REG_COLS, sprintf("%s -- disease regions", S)) +
       mk(d, "arch", ARCH_COLS, sprintf("%s -- architecture (f_B, within-sample quantiles)", S)) +
       plot_annotation(title = sprintf("%s regions, tier '%s'", S, TIER),
         subtitle = sprintf(paste0("%d cells drawn of %d | %d um bins, min %d cells | ",
           "hotspot cut %.1f cp10k | B_dense coherence p = %s"),
           nrow(d), md[samp==S, .N], FOCUS_BIN, MIN_CELLS, CUT,
           signif(coh[samp==S, p], 3)))
  sf(sprintf("fig_16_regions_%s.png", S), p, 15, 8)
}

writeLines(format(Sys.time()), file.path(OUT, "regions.done"))
cat("============================================================\nOK\n")
