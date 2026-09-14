## ============================================================================
## 18_script_leakage_views.R -- the donor-aligned density maps.
##
## THE FIGURE. One panel per (argmax cell type) x (agrees / disagrees with the cluster).
## Within a panel, every focal cell has been rotated so its A-marker donor is straight UP,
## and rescaled so the donor centroid sits at (0, 1). Thousands of cells are then pooled
## into one map. If argmax is being fooled by transcripts bleeding in from the neighbour,
## the disagreeing panels should show A-marker density pulled toward +y and the agreeing
## panels should not.
##
## WHAT IS PLOTTED IS A RATIO, NOT A DENSITY. Plotting raw A-marker density would mostly
## redraw the cell's overall transcript cloud -- bigger cells, more transcripts, more
## markers, in every direction. Each bin carries an A-marker SHARE relative to an
## expectation, so cell size, cell shape and depth divide out.
##
## THE EXPECTATION IS PER RADIAL SHELL, NOT GLOBAL, AND THAT MATTERS. A first version
## normalised each bin against the panel-wide A-marker share and produced twelve panels
## that were blue in the middle and red at the rim -- an isotropic RADIAL pattern that
## swamped the directional one. The cause is real and is the same gene-dependent
## compartment bias that rules out a naive CellComp test: the panel's marker genes are
## cytoplasm- and membrane-biased, while the average transcript is 79% nuclear, so markers
## sit further from the centroid than average transcripts do IN EVERY DIRECTION.
## Normalising within radial shells divides that out exactly and leaves only the angular
## asymmetry -- which is the entire question. Red is enrichment, blue depletion.
##
## COORDINATES ARE IN CELL RADII, NOT DONOR DISTANCES. Scaling by the donor distance made
## the pooled map inhomogeneous: cells with a close donor have a large radius/distance
## ratio and fill the outer window, cells with a far donor occupy only the centre, so the
## rim of the map was a different subset of cells than its middle. Scaling by each cell's
## own r_eff makes every cell contribute one unit disk, and the question -- where WITHIN
## the cell do the markers sit -- is asked the same way of every cell.
##
## AGREE AND DISAGREE ARE BALANCED BY WEIGHTING, NOT BY SUBSAMPLING. Disagreeing cells are
## shallower and carry fewer markers than agreeing ones in every lineage (B: 445 vs 668
## counts, k_A 2 vs 3) -- that is the shot-noise explanation competing with the leakage
## one, and unbalanced panels would show it rather than leakage.
##
## An earlier version balanced by taking min(n_agree, n_disagree) from each (argmax, k_A,
## depth-quintile) stratum. That is exact, but it threw away 62% of the eligible cells and
## could never have used more: the disagree pool caps a 1:1 design at 36,140 cells and it
## was already returning 35,288, so no loosening of the strata could have helped. Agreeing
## cells now carry an inverse-probability weight that reproduces the disagreeing cells'
## joint (k_A, depth) distribution instead, which balances the same two covariates while
## keeping every eligible cell. Effective sample sizes are reported below, because a
## weighted n is not an n.
##
## THE X AXIS IS A BUILT-IN CONTROL. Rotation fixes only the y direction; x is arbitrary,
## so every panel MUST be left-right symmetric. Visible x asymmetry means the frame is
## broken, not that something was found.
##
## Usage: sbatch --mem=96G --time=2:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/18_script_leakage_views.R IB1
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); setDTthreads(2)
})
args   <- commandArgs(trailingOnly = TRUE)
SAMPLE <- if (length(args) >= 1) args[1] else "IB1"
OUT    <- "output/18_leakage"
FIG    <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
FRM_F  <- file.path(OUT, sprintf("%s_tx_frame_focal.csv.gz", SAMPLE))
POL_F  <- file.path(OUT, sprintf("%s_polarisation.csv.gz", SAMPLE))
for (f in c(FRM_F, POL_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
DIV_HI<-"#d03b3b"; DIV_MID<-"#f0efec"; DIV_LO<-"#184f95"
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
          strip.text=element_text(colour=INK2,size=base-1,face="bold"),
          legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left", legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF);cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

D_ADJ  <- 20      # adjacency cap actually used for the figure
NBIN   <- 13L     # bins per axis -- twice the bin WIDTH of the 25-bin version, so four times
                  # the area and roughly four times the transcripts behind each colour. Family-A
                  # markers are the scarce channel (162k across all twelve panels), and they set
                  # the noise, so bin area is the cheapest lever on it.
LIM    <- 2.0     # window half-width, in units of the focal cell radius r_eff
NSHELL <- 20L     # radial shells the expectation is computed within
MIN_BIN<- 300L    # transcripts in a bin before it is drawn (bins are 4x larger now)
CLAMP  <- 0.6     # colour clamp, log2 units
## k_A is capped before being used as a balancing stratum. The cap was 6, which left the
## myeloid panel unbalanced in the one direction that matters: agreeing cells averaged 8.41
## markers against 6.59 for disagreeing ones, because everything above 6 pooled into a
## single bin that the agreeing cells dominated. More of a cell's own markers means a
## smaller RELATIVE contribution from any leaked ones, so residual imbalance of that sign
## deflates the agreeing row and inflates the disagree-minus-agree difference -- i.e. it
## flatters the hypothesis. 15 covers the myeloid and fibro tails.
KCAP   <- 15L

cat("============================================================\n")
cat("step 18 -- donor-aligned density maps: ", SAMPLE, "\n", sep = "")
cat("============================================================\n")

pol <- fread(POL_F)
pol <- pol[scoreable == TRUE & adjacent == TRUE & donorA_dist <= D_ADJ & k_A >= 1 & n_tx >= 20]
cat(sprintf("cells eligible: %d\n", nrow(pol)))

## ---- balance agree to disagree by weighting, within argmax call ------------
## Disagreeing cells keep weight 1 and define the target distribution; agreeing cells are
## weighted by the ratio of the two stratum shares. Strata with no disagreeing cells get
## weight 0, which is the correct answer -- there is nothing there to compare against.
set.seed(1L)
pol[, nq := cut(nCount_RNA, quantile(nCount_RNA, 0:5/5), include.lowest = TRUE, labels = FALSE),
    by = A]
pol[, kbin := pmin(k_A, KCAP)]
## as.numeric(.N), not .N: normalising an integer column in place silently truncates every
## share to 0 and zeroes the whole agreeing row.
tgt <- pol[agree == FALSE, .(p_d = as.numeric(.N)), by = .(A, nq, kbin)][, p_d := p_d / sum(p_d), by = A]
src <- pol[agree == TRUE,  .(p_a = as.numeric(.N)), by = .(A, nq, kbin)][, p_a := p_a / sum(p_a), by = A]
pol <- merge(merge(pol, tgt, by = c("A","nq","kbin"), all.x = TRUE),
             src, by = c("A","nq","kbin"), all.x = TRUE)
pol[, w := 1]
pol[agree == TRUE, w := fifelse(is.na(p_d) | is.na(p_a) | p_a == 0, 0, p_d / p_a)]
## Rescale so the weighted agree total equals the disagree count within each argmax call,
## keeping the two rows of the figure on a comparable footing. Done via a join rather than
## a self-referencing `:=` -- querying `pol` inside a by-group assignment to `pol` reads a
## table that is being modified underneath it.
nd <- pol[agree == FALSE, .(n_d = as.numeric(.N)), by = A]
sw <- pol[agree == TRUE,  .(s_w = sum(w)),         by = A]
pol <- merge(merge(pol, nd, by = "A", all.x = TRUE), sw, by = "A", all.x = TRUE)
pol[agree == TRUE & s_w > 0, w := w * n_d / s_w]
matched <- pol[w > 0]
if (matched[, uniqueN(agree)] < 2L) stop("weighting removed one of the two groups", call. = FALSE)
ess <- function(w) sum(w)^2 / sum(w^2)
cat("\n=== cells per panel (raw n, and effective n after weighting) ===\n")
print(matched[, .(cells = .N, ess = round(ess(w))), by = .(A, agree)][order(A, agree)])
cat(sprintf("\ntotal cells used: %d (the 1:1 matched design used 35288)\n", nrow(matched)))
wm <- function(x, w) sum(x * w) / sum(w)
cat("\n=== covariate balance after weighting ===\n")
print(matched[, .(mean_nCount = round(wm(nCount_RNA, w)), mean_k_A = round(wm(k_A, w), 2)),
              by = .(A, agree)][order(A, agree)])

## ---- the frame, restricted to matched cells and normalised by donor distance ----
frm <- fread(FRM_F, select = c("cell_id","fam","is_A","CellComp","xr","yr"))
key <- matched[, .(cell_id, fov_id, A, agree, r_eff_um, w)]
frm <- frm[key, on = .(cell_id), nomatch = NULL]
frm[, `:=`(xn = xr / r_eff_um, yn = yr / r_eff_um)]
frm <- frm[abs(xn) <= LIM & abs(yn) <= LIM]
frm[, rad := sqrt(xn^2 + yn^2)]
cat(sprintf("\ntranscripts in window: %d (%d family-A markers, %.0f weighted)\n",
            nrow(frm), frm[, sum(is_A)], frm[, sum(w * is_A)]))

## ---- THE PRIMARY STATISTIC: up/down asymmetry, bootstrapped over FOVs ------
## Toward the donor vs away from it, at matched radius by construction -- the +y and -y
## half-planes have identical radial profiles, so the compartment/radial artefact cancels
## exactly rather than approximately. FOV is the resampling unit: region is the unit of
## every p-value in this project, and 22.9M transcripts are not 22.9M observations.
frm[, side := ifelse(yn > 0, "up", "dn")]
cellside <- frm[yn != 0, .(n_all = sum(w), n_A = sum(w * is_A)), by = .(A, agree, fov_id, side)]
## Returns one value per row of `base`, matched BY KEY. Trusting dcast to return rows in
## the same order every time is the kind of assumption that silently transposes a result
## when one bootstrap sample happens to lack a panel.
asym_of <- function(d, keys) {
  w <- dcast(d[, .(n_all = sum(n_all), n_A = sum(n_A)), by = .(A, agree, side)],
             A + agree ~ side, value.var = c("n_all","n_A"))
  w[, a := log2((n_A_up / n_all_up) / (n_A_dn / n_all_dn))]
  w[keys, on = .(A, agree)]$a
}
base <- dcast(cellside[, .(n_all = sum(n_all), n_A = sum(n_A)), by = .(A, agree, side)],
              A + agree ~ side, value.var = c("n_all","n_A"))
base[, asym := log2((n_A_up / n_all_up) / (n_A_dn / n_all_dn))]
set.seed(2L)
fovs <- unique(cellside$fov_id)
BOOT <- 400L
keys <- base[, .(A, agree)]
bs <- vapply(seq_len(BOOT), function(i) {
  f <- data.table(fov_id = sample(fovs, length(fovs), replace = TRUE))
  asym_of(f[cellside, on = .(fov_id), nomatch = NULL, allow.cartesian = TRUE], keys)
}, numeric(nrow(base)))
base[, `:=`(lo = apply(bs, 1, quantile, .025, na.rm = TRUE),
            hi = apply(bs, 1, quantile, .975, na.rm = TRUE))]
cat(sprintf("\n=== A-marker enrichment TOWARD the donor vs away (log2), %d FOVs, %d bootstraps ===\n",
            length(fovs), BOOT))
print(base[, .(A, agree, n_A_up = round(n_A_up), n_A_dn = round(n_A_dn), asym = round(asym, 4),
               ci = sprintf("[%.4f, %.4f]", lo, hi),
               excludes_0 = (lo > 0) | (hi < 0))][order(A, agree)])
sc(sprintf("%s_updown_asymmetry.csv.gz", SAMPLE), base)

## THE HYPOTHESIS IS A DIFFERENCE, so it gets tested as one. Comparing two confidence
## intervals for overlap is not the same test and is the less powerful of the two: the
## same bootstrap resample feeds both terms, so the FOV-level noise they share cancels.
di <- match(keys[agree == FALSE, paste(A)], keys$A[!keys$agree])
dif <- data.table(A = keys[agree == FALSE, A],
                  diff = bs[which(!keys$agree), ][, 1] * 0)   # placeholder, filled below
iF <- which(!keys$agree); iT <- which(keys$agree)
ordF <- order(keys$A[iF]); ordT <- order(keys$A[iT])
dbs  <- bs[iF[ordF], , drop = FALSE] - bs[iT[ordT], , drop = FALSE]
dif  <- data.table(A = keys$A[iF][ordF],
                   obs = base$asym[iF][ordF] - base$asym[iT][ordT],
                   lo  = apply(dbs, 1, quantile, .025, na.rm = TRUE),
                   hi  = apply(dbs, 1, quantile, .975, na.rm = TRUE))
dif[, p_boot := 2 * pmin(rowMeans(dbs <= 0, na.rm = TRUE), rowMeans(dbs >= 0, na.rm = TRUE))]
cat("\n=== disagree MINUS agree (log2), paired on the same FOV resamples ===\n")
print(dif[, .(A, diff = round(obs, 4), ci = sprintf("[%.4f, %.4f]", lo, hi),
              p_boot = round(p_boot, 4))][order(-diff)])
sc(sprintf("%s_updown_difference.csv.gz", SAMPLE), dif)

brk <- seq(-LIM, LIM, length.out = NBIN + 1L)
sbrk <- seq(0, LIM * sqrt(2), length.out = NSHELL + 1L)
frm[, `:=`(bx = as.integer(cut(xn, brk)), by_ = as.integer(cut(yn, brk)),
           shell = as.integer(cut(rad, sbrk)))]
grid <- frm[!is.na(bx) & !is.na(by_) & !is.na(shell),
            .(n_all = sum(w), n_A = sum(w * is_A), shell = shell[1]), by = .(A, agree, bx, by_)]
## the expectation: the A-marker share at this RADIUS, in this panel
shl <- frm[!is.na(shell), .(p_shell = sum(w * is_A) / sum(w)), by = .(A, agree, shell)]
grid <- shl[grid, on = .(A, agree, shell)]
grid[, `:=`(x = (brk[bx] + brk[bx+1L])/2, y = (brk[by_] + brk[by_+1L])/2)]
grid[, lr := log2(((n_A + 0.5) / (n_all + 1)) / p_shell)]
grid <- grid[n_all >= MIN_BIN & p_shell > 0]
sc(sprintf("%s_density_grid.csv.gz", SAMPLE), grid)

## x must be symmetric: rotation fixes y only. Counts, not a mean of per-bin logs -- the
## first version averaged logs over bins, which the low-count bins dominated and which is
## not constrained to zero even when the map is perfectly symmetric.
sym <- frm[xn != 0, .(n_all = sum(w), n_A = sum(w * is_A)), by = .(A, agree, side = ifelse(xn > 0, "R", "L"))]
sym <- dcast(sym, A + agree ~ side, value.var = c("n_all","n_A"))
sym[, lr_LR := log2((n_A_R / n_all_R) / (n_A_L / n_all_L))]
cat("\n=== x-symmetry control (log2 right/left; must be ~0 -- x is arbitrary) ===\n")
print(sym[, .(A, agree, lr_LR = round(lr_LR, 4))][order(A, agree)])

ORD <- c("B","T","myeloid","fibro","endo","PC")
grid[, A := factor(A, levels = ORD)]
grid[, panel := factor(ifelse(agree, "agrees", "DISAGREES"),
                       levels = c("DISAGREES","agrees"))]
## Where the donor actually sits, now that coordinates are in cell radii rather than donor
## distances. Drawn only if it falls inside the window.
dmed <- matched[, .(dy = median(donorA_dist / r_eff_um)), by = .(A, agree)]
dmed[, panel := factor(ifelse(agree, "agrees", "DISAGREES"), levels = c("DISAGREES","agrees"))]
dmed[, A := factor(A, levels = ORD)]
cat(sprintf("\nmedian donor distance in cell radii: %.2f (drawn where <= %.1f)\n",
            median(dmed$dy), LIM))
ann <- unique(matched[, .(A = factor(A, levels = ORD),
                          panel = factor(ifelse(agree, "agrees", "DISAGREES"),
                                         levels = levels(grid$panel)))])
ann <- merge(ann, matched[, .(N = .N), by = .(A, agree)][, .(A = factor(A, levels = ORD), agree, N)],
             by = "A", allow.cartesian = TRUE)
ann <- unique(ann[(panel == "argmax agrees with cluster") == agree, .(A, panel, N)])

p <- ggplot(grid[!is.na(A)], aes(x, y, fill = lr)) +
  geom_raster() +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = .2) +
  geom_vline(xintercept = 0, colour = MUTED, linewidth = .2) +
  annotate("path", x = cos(seq(0, 2*pi, length.out = 181)),
                  y = sin(seq(0, 2*pi, length.out = 181)),
           colour = INK2, linewidth = .3, linetype = "22") +
  geom_point(data = dmed[dy <= LIM], aes(x = 0, y = dy), inherit.aes = FALSE,
             colour = INK, size = 1.5) +
  annotate("point", x = 0, y = 0, colour = INK, shape = 3, size = 1.6) +
  geom_text(data = ann, aes(x = -LIM*0.94, y = LIM*0.88, label = paste0("n=", N)),
            inherit.aes = FALSE, hjust = 0, size = 2.6, colour = INK2) +
  facet_grid(panel ~ A, switch = "y") +
  scale_fill_gradient2(low = DIV_LO, mid = DIV_MID, high = DIV_HI, midpoint = 0,
                       limits = c(-CLAMP, CLAMP), oob = scales::squish,
                       name = "log2 enrichment of family-A marker transcripts") +
  coord_fixed(xlim = c(-LIM, LIM), ylim = c(-LIM, LIM), expand = FALSE) +
  labs(title = sprintf("%s -- are argmax's marker transcripts displaced toward the neighbour that could have donated them?", SAMPLE),
       subtitle = paste0("Each cell rotated so its nearest marker-rich neighbour is straight UP, and scaled so its own radius is 1. Rows: does argmax agree with the cluster label?\n",
                         "Cross = cell centroid, dashed circle = cell boundary, dot = median donor centroid. Agree and disagree matched on argmax call, marker count and depth."),
       x = "across the donor axis (cell radii)",
       y = "toward the donor  ->",
       caption = paste0("Each bin is the family-A marker share relative to the share at the SAME RADIUS, so the radial (compartment) gradient is divided out and only angular asymmetry can show. ",
                        "Bins under ", MIN_BIN, " transcripts dropped; colour clamped at +/- ", CLAMP,
                        " log2.\nx is arbitrary by construction, so each panel must be left-right symmetric -- the measured left/right imbalance is the noise floor for reading the up/down one.")) +
  theme_viz(11, xgrid = TRUE) +
  theme(panel.grid.major = element_blank(), panel.spacing = unit(6, "pt"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, size = 10))
sf(sprintf("%s_leakage_density_maps.png", SAMPLE), p, 15, 7.4)

## ---- the same thing as a 1D profile along the donor axis -------------------
## Collapsed across x from the shell-normalised grid, so the radial component is already
## divided out and what remains is purely the toward-vs-away asymmetry.
prof <- grid[, .(n_all = sum(n_all), n_A = sum(n_A),
                 p = weighted.mean(p_shell, n_all)), by = .(A, agree, y)]
prof[, lr := log2(((n_A + 0.5) / (n_all + 1)) / p)]
prof <- prof[n_all >= MIN_BIN]
p2 <- ggplot(prof[!is.na(A)], aes(y, lr, colour = agree)) +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = .3) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = .3, linetype = "22") +
  geom_line(linewidth = .6) + facet_wrap(~ A, nrow = 1) +
  scale_colour_manual(values = c(`TRUE` = "#2a78d6", `FALSE` = "#eb6834"),
                      labels = c(`TRUE` = "agrees", `FALSE` = "disagrees"), name = NULL) +
  labs(title = sprintf("%s -- marker enrichment along the donor axis", SAMPLE),
       subtitle = "Collapsed across x, with the radial gradient already divided out. Positive x is toward the donor.",
       x = "toward the donor (cell radii)", y = "log2 enrichment of family-A markers") +
  theme_viz(11, xgrid = TRUE)
sf(sprintf("%s_leakage_profile.png", SAMPLE), p2, 15, 3.6)
cat("\ndone\n")
