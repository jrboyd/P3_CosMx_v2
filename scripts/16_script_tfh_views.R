## ============================================================================
## 16_script_tfh_views.R -- representative Tfh cells, in the image.
##
## Three panels per window, and a tight and a context view of each cell:
##   1. the colorized multichannel image, with an inset locating the window on the slide
##   2. cell boundaries filled by CELL TYPE, so the Tfh call can be checked against morphology
##   3. transcript points -- PDCD1 in its own colour, the related checkpoint/Tfh probes
##      combined into a single second colour, per the brief
## A red box marks the centre Tfh cell in every panel so the same cell is findable across all
## three, and across both zooms.
##
## "REPRESENTATIVE" IS A RULE, NOT A PICK. Choosing the most extreme cell would make the
## figure a best case. The rule is: the cell is called Tfh, carries PDCD1 >= 2 AND CXCL13 >= 2
## (the calibrated threshold from qc_16_negprobe_floor.R -- below 2 counts these genes sit at
## the blank-probe floor, so a panel-3 with nothing in it would be uninformative rather than
## negative), its depth sits inside the interquartile range of Tfh cells in the same region
## (so it is neither a shallow nor an unusually deep cell), it is >= EDGE_UM from an FOV seam,
## and among the survivors it has the most Tfh neighbours within NBR_UM. That last term picks
## a cell sitting in Tfh-rich tissue rather than an isolated call, which is what makes it
## representative of the population rather than of the tail.
##
## COORDINATES: the mm->px affine is FITTED, not constant, and is asserted on the three things
## that make a bad fit look good (both slopes negative, |slope| within 1% of 1/COSMX_MM_PER_PX,
## residual under 1 px) -- see CLAUDE.md and 11_script_hodgkin_niche_views.R:128-160. Panels
## use scale_x_reverse(), never reversed coord limits, because ggplot sorts coord limits.
##
## Usage: Rscript scripts/16_script_tfh_views.R [n_per_group] [tier]
## ============================================================================

if (!file.exists("scripts/lib/cosmx_tiff.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(Matrix); library(data.table); setDTthreads(1)
})

args   <- commandArgs(trailingOnly = TRUE)
N_PER  <- if (length(args) >= 1) as.integer(args[1]) else 2L
TIER   <- if (length(args) >= 2) args[2] else "matched"

FLAT    <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
REG_F   <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT_F    <- "output/16_celltypes/cell_celltypes.csv.gz"
CNT_F   <- "output/16_regions/panel_counts_all.Rds"
OUT_DIR <- file.path("output/16_tfh_views", TIER)
for (f in c(REG_F, CT_F, CNT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TIGHT_UM <- 60; CONTEXT_UM <- 320
BOX_UM   <- 7          # half-width of the marker box around the centre cell
EDGE_UM  <- 40         # keep the window off an FOV seam
NBR_UM   <- 50         # radius for the "Tfh neighbours" tie-break
THR      <- 2L         # calibrated positivity threshold
LOC_N    <- 60000L
TARGETS  <- list(IB1 = c("FL_focus","cHL_hotspot"), RB1 = c("FL_follicular","interfollicular"))

## Panel 3: PD-1 alone, and everything else in one colour, as asked.
G_PD1   <- "PDCD1"
G_OTHER <- c("CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","ICOS","CXCL13","CXCR5")
GENES   <- c(G_PD1, G_OTHER)
PD1_COL <- "#ffd400"; OTH_COL <- "#00e5ff"
CT_COLS <- c(Tfh = "#eb6834", `other T` = "#1baf7a", B = "#2a78d6",
             myeloid = "#f2b705", other = "#d8d6ce")
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"
theme_v <- theme_minimal(base_size = 10) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA), panel.grid=element_blank(),
        axis.text=element_text(colour=MUTED,size=7), axis.title=element_text(colour=INK2,size=8),
        plot.title=element_text(colour=INK2,face="bold",size=9),
        legend.text=element_text(colour=INK2,size=7.5), legend.title=element_blank(),
        legend.position="bottom", legend.key.height=unit(8,"pt"))
UMPX <- COSMX_MM_PER_PX * 1000
.um  <- function(v) round(v * UMPX)
.coord_px <- function(rect) { q <- rect@coords[1L, ]
  list(scale_x_reverse(labels = .um), coord_fixed(xlim=c(q$xmin,q$xmax), ylim=c(q$ymin,q$ymax))) }

cat("============================================================\n")
cat("step 16 views -- representative Tfh cells, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

## ---- cells, types, counts ----------------------------------------------------
rg <- fread(REG_F, select=c("cell_id","samp","bin_region","x_slide_mm","y_slide_mm",
                            "is_defect","unbinned","fov_id","nCount_RNA"))
ct <- fread(CT_F, select=c("cell_id","celltype"))
d  <- merge(rg, ct, by="cell_id", sort=FALSE)[is_defect==FALSE & unbinned==FALSE]
M  <- readRDS(CNT_F)
jj <- match(d$cell_id, colnames(M)); stopifnot(!anyNA(jj))
for (g in intersect(GENES, rownames(M))) d[, (paste0("ct_", g)) := as.vector(M[g, jj])]
rm(M); invisible(gc())
d[, key := cosmx_seu_key(cell_id)]
cat(sprintf("cells %d | Tfh %d\n", nrow(d), d[celltype=="Tfh", .N]))

SAMPS <- names(TARGETS)
FLATS <- lapply(SAMPS, function(S) {
  f <- fread(sprintf(FLAT, S, S), select=c("cell","fov","CenterX_local_px","CenterY_local_px"))
  px <- cosmx_local_to_px(S, f$fov, f$CenterX_local_px, f$CenterY_local_px)
  f[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)][]
}); names(FLATS) <- SAMPS

## Fitted mm->px affine, asserted. Fit on CELLS, never on FOV origins (the origins are
## exactly linear in mm by construction, so a residual of 0 there confirms nothing).
FIT <- lapply(SAMPS, function(S) {
  m <- merge(d[samp==S, .(key, x_slide_mm, y_slide_mm)], FLATS[[S]][, .(key, X_px, Y_px)], by="key")
  stopifnot(nrow(m) > 1000)
  fx <- lm(X_px ~ x_slide_mm, data=m); fy <- lm(Y_px ~ y_slide_mm, data=m)
  a <- list(ax=unname(coef(fx)[2]), bx=unname(coef(fx)[1]),
            ay=unname(coef(fy)[2]), by=unname(coef(fy)[1]),
            res=c(sd(resid(fx)), sd(resid(fy))), n=nrow(m))
  ex <- 1/COSMX_MM_PER_PX
  if (a$ax > 0 || a$ay > 0) stop(S, ": mm->px slopes must both be negative (mirrored x); got ",
                                 round(a$ax), "/", round(a$ay), call.=FALSE)
  if (max(abs(abs(c(a$ax,a$ay)) - ex)/ex) > 0.01)
    stop(S, ": slope ", round(a$ax), "/", round(a$ay), " is >1% from ", round(ex), call.=FALSE)
  if (max(a$res) > 1) stop(S, ": affine residual SD ", signif(max(a$res),3),
                           " px exceeds 1 px -- do not trust it", call.=FALSE)
  cat(sprintf("%s affine: ax %.1f ay %.1f | resid %.2f/%.2f px | n=%d\n",
              S, a$ax, a$ay, a$res[1], a$res[2], a$n)); a
}); names(FIT) <- SAMPS
mm2px <- function(S,x,y){a<-FIT[[S]]; list(X_px=a$ax*x+a$bx, Y_px=a$ay*y+a$by)}
px2mm <- function(S,X,Y){a<-FIT[[S]]; list(x=(X-a$bx)/a$ax, y=(Y-a$by)/a$ay)}
CANVAS <- lapply(SAMPS, function(S){pl <- cosmx_fov_placement(S)
  c(x=max(pl$x0)+COSMX_FOV_PX-1L, y=max(pl$y0)+COSMX_FOV_PX-1L)}); names(CANVAS) <- SAMPS

## FOV seam distance, from empirical per-FOV mm boxes (as step 11 does)
FOVBOX <- d[, .(x1=quantile(x_slide_mm,.02), x2=quantile(x_slide_mm,.98),
                y1=quantile(y_slide_mm,.02), y2=quantile(y_slide_mm,.98), n=.N),
            by=.(samp,fov_id)][n>=50]

## ---- selection ---------------------------------------------------------------
sel <- rbindlist(lapply(SAMPS, function(S) {
 rbindlist(lapply(TARGETS[[S]], function(R) {
  pool <- d[samp==S & bin_region==R & celltype=="Tfh"]
  if (!nrow(pool)) { cat(sprintf("  %s / %s: no Tfh cells\n", S, R)); return(NULL) }
  n0 <- nrow(pool)
  if (all(c("ct_PDCD1","ct_CXCL13") %in% names(pool)))
    pool <- pool[ct_PDCD1 >= THR & ct_CXCL13 >= THR]
  n1 <- nrow(pool)
  qd <- quantile(d[samp==S & bin_region==R & celltype=="Tfh", nCount_RNA], c(.25,.75))
  pool <- pool[nCount_RNA >= qd[1] & nCount_RNA <= qd[2]]
  n2 <- nrow(pool)
  fb <- FOVBOX[samp==S]
  if (nrow(pool) && nrow(fb)) {
    nn <- RANN::nn2(as.matrix(fb[, .((x1+x2)/2,(y1+y2)/2)]),
                    as.matrix(pool[, .(x_slide_mm,y_slide_mm)]), k=1)$nn.idx[,1]
    ff <- fb[nn]
    pool[, seam_um := 1000*pmin(abs(x_slide_mm-ff$x1), abs(x_slide_mm-ff$x2),
                                abs(y_slide_mm-ff$y1), abs(y_slide_mm-ff$y2))]
    pool <- pool[seam_um >= EDGE_UM]
  }
  n3 <- nrow(pool)
  cat(sprintf("  %s / %-16s Tfh %5d -> PDCD1&CXCL13>=%d %4d -> depth IQR %4d -> off-seam %4d\n",
              S, R, n0, THR, n1, n2, n3))
  if (!n3) return(NULL)
  ## tie-break: most Tfh neighbours within NBR_UM -> a cell in Tfh-rich tissue
  allt <- d[samp==S & celltype=="Tfh", .(x_slide_mm, y_slide_mm)]
  cnt <- RANN::nn2(as.matrix(allt), as.matrix(pool[, .(x_slide_mm,y_slide_mm)]),
                   k=min(50L, nrow(allt)))$nn.dists*1000
  pool[, n_tfh_near := rowSums(cnt <= NBR_UM) - 1L]
  head(pool[order(-n_tfh_near, -nCount_RNA)], N_PER)[, `:=`(samp=S, bin_region=R)][]
 }), fill=TRUE)
}), fill=TRUE)
if (!nrow(sel)) stop("no representative Tfh cell passed the rule", call.=FALSE)
cat(sprintf("\nselected %d cells\n", nrow(sel)))
print(sel[, .(samp, bin_region, cell_id, nCount_RNA, ct_PDCD1, ct_CXCL13, n_tfh_near,
              seam_um=round(seam_um))])
fwrite(sel, file.path(OUT_DIR,"selected_tfh_cells.csv"))

## ---- per-cell type map for panel 2 ------------------------------------------
d[, ct5 := fifelse(celltype=="Tfh","Tfh",
           fifelse(celltype %in% c("Treg","CD8_T","CD4_conv","T_unresolved"),"other T",
           fifelse(celltype=="B","B",
           fifelse(celltype %in% c("Macrophage","Monocyte","cDC","myeloid_unresolved"),
                   "myeloid","other"))))]
CTMAP <- lapply(SAMPS, function(S) d[samp==S, .(key, ct5)]); names(CTMAP) <- SAMPS
LOCBG <- lapply(SAMPS, function(S) d[samp==S][sample.int(.N, min(.N, LOC_N)),
                                              .(x_slide_mm, y_slide_mm)]); names(LOCBG) <- SAMPS

GPCACHE <- new.env(parent=emptyenv())
gp_for_rect <- function(S, rect) {
  fv <- cosmx_fovs_in_rect(S, rect); if (is.null(fv) || !nrow(fv)) return(NULL)
  parts <- lapply(fv$fov, function(f) {
    k <- paste0(S,"_",f)
    if (is.null(GPCACHE[[k]])) {
      row <- fv[fv$fov==f, , drop=FALSE]
      fr <- cosmx_px_rect(row$x0, row$x0+COSMX_FOV_PX, row$y0, row$y0+COSMX_FOV_PX,
                          name=paste0("fov",f))
      g <- tryCatch(cosmx_gene_points(S, fr, GENES), error=function(e) NULL)
      GPCACHE[[k]] <- if (is.null(g)) data.frame() else as.data.frame(g)
    }
    GPCACHE[[k]]
  })
  parts <- parts[vapply(parts, nrow, 0L) > 0L]
  if (!length(parts)) NULL else do.call(rbind, parts)
}

## ---- panels ------------------------------------------------------------------
## Slide-mm panel: NEITHER axis reversed. The affine slopes are both negative, so mm runs
## opposite to canvas px on both axes -- px panels reverse both, mm panels reverse neither.
## Matches step 12's locator (12_script_fl_vs_hodgkin_overlay.R:147). The first version
## reversed y here, which mirrored the locator against the image panel it sits inside.
.locator <- function(S, cx_mm, cy_mm) {
  ggplot() +
    geom_point(data=LOCBG[[S]], aes(x_slide_mm,y_slide_mm), colour="#cfcec9", size=.04) +
    annotate("point", x=cx_mm, y=cy_mm, colour="#c1121f", size=2.2) +
    coord_fixed() + theme_void() +
    theme(plot.background=element_rect(fill="white", colour="#b5b3ac", linewidth=.3))
}
.boxann <- function(bx) annotate("rect", xmin=bx$x1, xmax=bx$x2, ymin=bx$y1, ymax=bx$y2,
                                 colour="#c1121f", fill=NA, linewidth=.8)
.p_image <- function(S, rect, bx, ttl) {
  img <- tryCatch(cosmx_image(S, rect), error=function(e) NULL)
  if (is.null(img)) return(ggplot()+theme_v+labs(title=paste0(ttl," (image unavailable)")))
  img@plots$colorized + .boxann(bx) + .coord_px(rect) + scale_y_reverse(labels=.um) +
    labs(x="x (um)", y="y (um)", title=ttl) + theme_v
}
.p_cells <- function(S, rect, bx, ttl) {
  poly <- tryCatch(cosmx_cells_in_rect(S, rect), error=function(e) NULL)
  if (is.null(poly)) return(list(p=ggplot()+theme_v+labs(title=paste0(ttl," (no cells)")), df=NULL))
  df <- poly@coords %>% dplyr::group_by(name) %>% dplyr::reframe(x=unlist(x), y=unlist(y))
  pk <- cosmx_poly_key(df$name); cm <- CTMAP[[S]]
  v  <- cm$ct5[match(pk, cm$key)]
  df$ct5 <- ifelse(is.na(v), "other", v)
  p <- ggplot(df, aes(x,y,group=name)) +
    geom_polygon(aes(fill=ct5), colour="#5a5a5a", linewidth=.18) +
    scale_fill_manual(values=CT_COLS, limits=names(CT_COLS), drop=FALSE, name=NULL) +
    .boxann(bx) + .coord_px(rect) + scale_y_reverse(labels=.um) +
    labs(x="x (um)", y="y (um)", title=ttl) + theme_v
  list(p=p, df=df)
}
.p_genes <- function(S, rect, bx, df, ttl, sz) {
  gp <- gp_for_rect(S, rect)
  p <- ggplot()
  if (!is.null(df) && nrow(df))
    p <- p + geom_polygon(data=df, aes(x,y,group=name), fill="#f4f3ee", colour="#33333333",
                          linewidth=.12)
  n1 <- 0L; n2 <- 0L
  if (!is.null(gp) && nrow(gp)) {
    g <- as.data.table(gp)
    q <- rect@coords[1L, ]
    g <- g[X_px>=q$xmin & X_px<=q$xmax & Y_px>=q$ymin & Y_px<=q$ymax]
    g[, grp := fifelse(target==G_PD1, "PD-1 (PDCD1)", "related probes")]
    n1 <- g[grp=="PD-1 (PDCD1)", .N]; n2 <- g[grp=="related probes", .N]
    if (n2) p <- p + geom_point(data=g[grp=="related probes"], aes(X_px,Y_px,colour=grp),
                                size=1.0*sz, stroke=0, alpha=.85)
    if (n1) p <- p + geom_point(data=g[grp=="PD-1 (PDCD1)"], aes(X_px,Y_px,colour=grp),
                                size=2.4*sz, stroke=0)
    p <- p + scale_colour_manual(values=c(`PD-1 (PDCD1)`=PD1_COL, `related probes`=OTH_COL),
                                 limits=c("PD-1 (PDCD1)","related probes"), drop=FALSE,
                                 name=NULL, guide=guide_legend(override.aes=list(size=3)))
  }
  p + .boxann(bx) + .coord_px(rect) + scale_y_reverse(labels=.um) +
    labs(x="x (um)", y="y (um)",
         title=sprintf("%s -- PDCD1 %d, related %d", ttl, n1, n2)) + theme_v
}

## ---- build --------------------------------------------------------------------
mkrect <- function(S, cx, cy, um, nm) {
  h <- (um/2)/UMPX; cv <- CANVAS[[S]]
  cosmx_px_rect(max(0, cx-h), min(cv[["x"]], cx+h),
                max(0, cy-h), min(cv[["y"]], cy+h), name=nm)
}
for (i in seq_len(nrow(sel))) {
  r <- sel[i]; S <- r$samp
  ctr <- mm2px(S, r$x_slide_mm, r$y_slide_mm)
  cx <- ctr$X_px; cy <- ctr$Y_px
  hb <- BOX_UM/UMPX
  bx <- list(x1=cx-hb, x2=cx+hb, y1=cy-hb, y2=cy+hb)
  tag <- sprintf("%s_%s_%s", S, r$bin_region, sub("^c_","",r$cell_id))
  cat("\n-- ", tag, "\n", sep="")
  panels <- list()
  for (zn in c("tight","context")) {
    um   <- if (zn=="tight") TIGHT_UM else CONTEXT_UM
    rect <- mkrect(S, cx, cy, um, zn)
    sz   <- if (zn=="tight") 1.6 else 0.8
    p1 <- .p_image(S, rect, bx, sprintf("1. image (%s, %g um)", zn, um))
    if (zn=="tight")
      p1 <- p1 + inset_element(.locator(S, r$x_slide_mm, r$y_slide_mm),
                               left=.62, bottom=.62, right=.995, top=.995, align_to="panel")
    cc <- .p_cells(S, rect, bx, sprintf("2. cell types (%s)", zn))
    p3 <- .p_genes(S, rect, bx, cc$df, sprintf("3. probes (%s)", zn), sz)
    panels <- c(panels, list(p1, cc$p, p3))
  }
  pg <- wrap_plots(panels, nrow=2, ncol=3) +
    plot_annotation(
      title = sprintf("%s -- representative Tfh cell in %s", S, r$bin_region),
      subtitle = sprintf(paste0("cell %s | nCount %d | PDCD1 %d, CXCL13 %d counts | %d Tfh ",
        "within %g um | %.0f um from an FOV seam\nred box = this cell (+/-%g um). ",
        "Selection rule is in the script header -- not a hand-picked best case."),
        r$cell_id, r$nCount_RNA, r$ct_PDCD1, r$ct_CXCL13, r$n_tfh_near, NBR_UM,
        r$seam_um, BOX_UM),
      theme = theme(plot.title=element_text(face="bold", colour=INK, size=13),
                    plot.subtitle=element_text(colour=INK2, size=8.5),
                    plot.background=element_rect(fill=SURF, colour=NA)))
  f <- file.path(OUT_DIR, sprintf("tfh_%s.png", tag))
  ggsave(f, pg, width=15, height=10.5, dpi=170, bg=SURF, limitsize=FALSE)
  cat("wrote ", f, "\n", sep="")
}
writeLines(format(Sys.time()), file.path(OUT_DIR,"tfh_views.done"))
cat("============================================================\nOK\n")
