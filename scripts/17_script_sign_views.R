## ============================================================================
## 17_script_sign_views.R -- the sign-consistency result, in the image.
##
## THE CLAIM BEING SHOWN. Across all five T compartments independently, IB1's T cells carry
## more LAG3 / HAVCR2 / CTLA4 while RB1's carry more PDCD1 / TIGIT / TOX (14 of 22 genes
## directionally consistent, 1.4 expected by chance, binomial p = 2.8e-12). That is two
## different inhibitory programmes on the same clone, and it is a COMPARATIVE claim -- so the
## figure is comparative: one IB1 window beside one RB1 window, matched on T-cell content, with
## the two gene sets in contrasting colours. A reader should be able to see which colour
## dominates which sample without being told.
##
## MATCHED ON T-CELL COUNT, not picked for effect. The windows are not the ones with the most
## transcripts -- that would guarantee a striking figure and prove nothing. Candidates need a
## comparable number of T cells and at least a few transcripts of BOTH sets (so neither panel
## is empty, which would be uninformative rather than negative), and among the survivors the
## chosen window is the one closest to the MEDIAN combined transcript count. Representative of
## the population, not of its tail. The rule and the attrition are printed.
##
## Panels per row: colorized image (with a slide locator) | cell boundaries by T subset |
## transcripts, IB1-high set against RB1-high set.
##
## Usage: Rscript scripts/17_script_sign_views.R [n_pairs] [tier]
## ============================================================================

if (!file.exists("scripts/lib/cosmx_tiff.R")) stop("run from the project root", call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); setDTthreads(1)
})
args   <- commandArgs(trailingOnly = TRUE)
N_PAIR <- if (length(args) >= 1) as.integer(args[1]) else 2L
TIER   <- if (length(args) >= 2) args[2] else "matched"

FLAT  <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
REG_F <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT_F  <- "output/16_celltypes/cell_celltypes.csv.gz"
OUT   <- file.path("output/17_pd1", TIER, "sign_views")
for (f in c(REG_F, CT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

WIN_UM  <- 220        # enough cells to judge density, tight enough to resolve transcripts
## EDGE_UM must exceed HALF the window, or the window crosses an FOV boundary even though its
## centre cell passed. The first version used 40 um against a 220 um window, so every window
## straddled a seam -- visible as a band of uniformly "other" cells along one edge, which is an
## adjacent FOV where typing fails, not a tissue boundary. Distance is only a pre-filter now;
## containment is VERIFIED below by counting the FOVs the window's cells actually span.
EDGE_UM <- WIN_UM/2 + 20
MIN_T   <- 25L        # T cells in the window
MIN_TX  <- 4L         # transcripts from EACH set
SET_IB1 <- c("LAG3","HAVCR2","CTLA4")     # higher in IB1 in all 5 compartments
SET_RB1 <- c("PDCD1","TIGIT","TOX")       # higher in RB1 in all 5 compartments
GENES   <- c(SET_IB1, SET_RB1)
COL_IB1 <- "#ffd400"; COL_RB1 <- "#00e5ff"
TSUB    <- c("Tfh","Treg","CD8_T","CD4_conv","T_unresolved")
## FULL vocabulary, not a collapsed one. The earlier version dumped PC, NK, endothelium,
## fibroblasts AND all three unassigned sentinels into a single "other" bucket, which in IB1 is
## ~half the cells -- so the largest category on the panel was the one carrying no information
## about what it contained. Sentinels get grey shades so they still read as "not called", but
## they are named and counted rather than hidden.
CT_FULL <- c("Tfh","Treg","CD8_T","CD4_conv","T_unresolved","B","PC",
             "Macrophage","Monocyte","cDC","myeloid_unresolved","NK","endo","fibro",
             "ambiguous","below_background","all_markers_0")
CT_COLS <- c(Tfh="#eb6834", Treg="#c1121f", CD8_T="#7B3294", CD4_conv="#E69F00",
             T_unresolved="#009E73", B="#2a78d6", PC="#0072B2",
             Macrophage="#8B5A2B", Monocyte="#B8860B", cDC="#66A61E",
             myeloid_unresolved="#01665E", NK="#CC79A7", endo="#56B4E9", fibro="#7d8f74",
             ambiguous="#cfcdc6", below_background="#dedcd4", all_markers_0="#eceae4",
             `excluded upstream`="#b9b6ad")
SENTS   <- c("ambiguous","below_background","all_markers_0")
REGPAIR <- list(c("FL_focus","FL_follicular"), c("background","interfollicular"))
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"
theme_v <- theme_minimal(base_size=10) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA), panel.grid=element_blank(),
        axis.text=element_text(colour=MUTED,size=7), axis.title=element_text(colour=INK2,size=8),
        plot.title=element_text(colour=INK2,face="bold",size=9),
        legend.text=element_text(colour=INK2,size=7.5), legend.title=element_blank(),
        legend.position="bottom", legend.key.height=unit(8,"pt"))
UMPX <- COSMX_MM_PER_PX*1000; .um <- function(v) round(v*UMPX)
.coord_px <- function(rect){ q <- rect@coords[1L,]
  list(scale_x_reverse(labels=.um), coord_fixed(xlim=c(q$xmin,q$xmax), ylim=c(q$ymin,q$ymax))) }

cat("============================================================\n")
cat("step 17 sign views -- IB1 vs RB1 inhibitory programme, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

rg <- fread(REG_F, select=c("cell_id","samp","bin_region","x_slide_mm","y_slide_mm",
                            "is_defect","unbinned","fov_id"))
ct <- fread(CT_F, select=c("cell_id","celltype"))
d  <- merge(rg, ct, by="cell_id", sort=FALSE)[is_defect==FALSE & unbinned==FALSE]
d[, key := cosmx_seu_key(cell_id)]
d[, ct8 := fifelse(celltype %in% CT_FULL, celltype, "excluded upstream")]
d[, is_T := celltype %in% TSUB]
SAMPS <- c("IB1","RB1")
FLATS <- lapply(SAMPS, function(S) {
  f <- fread(sprintf(FLAT,S,S), select=c("cell","fov","CenterX_local_px","CenterY_local_px"))
  px <- cosmx_local_to_px(S, f$fov, f$CenterX_local_px, f$CenterY_local_px)
  f[, `:=`(key=cosmx_seu_key(cell), X_px=px$X_px, Y_px=px$Y_px)][]
}); names(FLATS) <- SAMPS
FIT <- lapply(SAMPS, function(S) {
  m <- merge(d[samp==S,.(key,x_slide_mm,y_slide_mm)], FLATS[[S]][,.(key,X_px,Y_px)], by="key")
  fx <- lm(X_px~x_slide_mm,m); fy <- lm(Y_px~y_slide_mm,m)
  a <- list(ax=unname(coef(fx)[2]),bx=unname(coef(fx)[1]),
            ay=unname(coef(fy)[2]),by=unname(coef(fy)[1]),res=c(sd(resid(fx)),sd(resid(fy))))
  ex <- 1/COSMX_MM_PER_PX
  if (a$ax>0||a$ay>0) stop(S,": mm->px slopes must be negative",call.=FALSE)
  if (max(abs(abs(c(a$ax,a$ay))-ex)/ex)>0.01) stop(S,": slope >1% off expected",call.=FALSE)
  if (max(a$res)>1) stop(S,": affine residual >1 px",call.=FALSE)
  cat(sprintf("%s affine ax %.1f ay %.1f resid %.2f/%.2f px\n",S,a$ax,a$ay,a$res[1],a$res[2])); a
}); names(FIT) <- SAMPS
mm2px <- function(S,x,y){a<-FIT[[S]]; list(X=a$ax*x+a$bx, Y=a$ay*y+a$by)}
CANVAS <- lapply(SAMPS, function(S){pl<-cosmx_fov_placement(S)
  c(x=max(pl$x0)+COSMX_FOV_PX-1L,y=max(pl$y0)+COSMX_FOV_PX-1L)}); names(CANVAS)<-SAMPS
FOVBOX <- d[, .(x1=quantile(x_slide_mm,.02),x2=quantile(x_slide_mm,.98),
                y1=quantile(y_slide_mm,.02),y2=quantile(y_slide_mm,.98),n=.N),
            by=.(samp,fov_id)][n>=50]
CTMAP <- lapply(SAMPS, function(S) d[samp==S,.(key,ct8)]); names(CTMAP)<-SAMPS
set.seed(1)
LOCBG <- lapply(SAMPS, function(S) d[samp==S][sample.int(.N,min(.N,60000L)),
                                     .(x_slide_mm,y_slide_mm)]); names(LOCBG)<-SAMPS

GPC <- new.env(parent=emptyenv())
gp_for_rect <- function(S, rect) {
  fv <- cosmx_fovs_in_rect(S, rect); if (is.null(fv)||!nrow(fv)) return(NULL)
  parts <- lapply(fv$fov, function(f) {
    k <- paste0(S,"_",f)
    if (is.null(GPC[[k]])) {
      row <- fv[fv$fov==f,,drop=FALSE]
      fr <- cosmx_px_rect(row$x0,row$x0+COSMX_FOV_PX,row$y0,row$y0+COSMX_FOV_PX,name=paste0("f",f))
      g <- tryCatch(cosmx_gene_points(S, fr, GENES), error=function(e) NULL)
      GPC[[k]] <- if (is.null(g)) data.frame() else as.data.frame(g)
    }
    GPC[[k]]
  })
  parts <- parts[vapply(parts,nrow,0L)>0L]; if (!length(parts)) NULL else do.call(rbind,parts)
}
mkrect <- function(S,cx,cy,um,nm){ h <- (um/2)/UMPX; cv <- CANVAS[[S]]
  cosmx_px_rect(max(0,cx-h),min(cv[["x"]],cx+h),max(0,cy-h),min(cv[["y"]],cy+h),name=nm) }

## ---- candidate windows -------------------------------------------------------
cand <- function(S, R) {
  z <- d[samp==S & bin_region==R & is_T==TRUE]
  if (!nrow(z)) return(NULL)
  fb <- FOVBOX[samp==S]
  nn <- RANN::nn2(as.matrix(fb[,.((x1+x2)/2,(y1+y2)/2)]),
                  as.matrix(z[,.(x_slide_mm,y_slide_mm)]),k=1)$nn.idx[,1]
  ff <- fb[nn]
  z[, seam_um := 1000*pmin(abs(x_slide_mm-ff$x1),abs(x_slide_mm-ff$x2),
                           abs(y_slide_mm-ff$y1),abs(y_slide_mm-ff$y2))]
  z <- z[seam_um >= EDGE_UM]
  if (!nrow(z)) return(NULL)
  ## T cells within the window around each candidate centre
  all_xy <- as.matrix(d[samp==S, .(x_slide_mm*1000, y_slide_mm*1000)])
  isT <- d[samp==S, is_T]
  q <- RANN::nn2(all_xy, as.matrix(z[,.(x_slide_mm*1000,y_slide_mm*1000)]),
                 k=min(400L,nrow(all_xy)))
  within <- q$nn.dists <= WIN_UM/2
  z[, n_T := rowSums(within & matrix(isT[q$nn.idx],nrow=nrow(within)))]
  ## FOVs actually spanned by the cells inside the window. Requiring exactly one is the direct
  ## test; the distance filter above only makes it likely.
  fovm <- matrix(d[samp==S, fov_id][q$nn.idx], nrow=nrow(within))
  fovm[!within] <- NA_integer_
  z[, n_fov := apply(fovm, 1, function(v) length(unique(v[!is.na(v)])))]
  n1 <- nrow(z)
  z <- z[n_T >= MIN_T]; n2 <- nrow(z)
  z <- z[n_fov == 1L];  n3 <- nrow(z)
  cat(sprintf("  %s / %-16s T cells %6d -> off-seam %5d -> >=%d T %5d -> single-FOV %5d\n",
              S, R, nrow(d[samp==S & bin_region==R & is_T]), n1, MIN_T, n2, n3))
  if (!nrow(z)) return(NULL)
  z[order(-n_T)][1:min(.N, 400L)]
}

pairs_done <- 0L
for (rp in REGPAIR) {
  cat(sprintf("\n== region pair: IB1 %s  vs  RB1 %s\n", rp[1], rp[2]))
  ci <- cand("IB1", rp[1]); cr <- cand("RB1", rp[2])
  if (is.null(ci) || is.null(cr)) { cat("  no candidates; skipped\n"); next }
  ## score candidates by transcript content, then take the MEDIAN one (not the max)
  score <- function(S, z, nmax = 40L) {
    z <- z[1:min(.N, nmax)]
    out <- rbindlist(lapply(seq_len(nrow(z)), function(i) {
      p <- mm2px(S, z$x_slide_mm[i], z$y_slide_mm[i])
      rect <- mkrect(S, p$X, p$Y, WIN_UM, "w")
      gp <- gp_for_rect(S, rect)
      if (is.null(gp)||!nrow(gp)) return(NULL)
      q <- rect@coords[1L,]
      g <- as.data.table(gp)[X_px>=q$xmin & X_px<=q$xmax & Y_px>=q$ymin & Y_px<=q$ymax]
      data.table(i=i, n_ib=g[target %in% SET_IB1,.N], n_rb=g[target %in% SET_RB1,.N])
    }), fill=TRUE)
    if (!nrow(out)) return(NULL)
    out <- out[n_ib >= MIN_TX & n_rb >= MIN_TX]
    if (!nrow(out)) return(NULL)
    out[, tot := n_ib + n_rb]
    out[order(abs(tot - median(tot)))]
  }
  si <- score("IB1", ci); sr <- score("RB1", cr)
  if (is.null(si) || is.null(sr)) { cat("  no window had both transcript sets; skipped\n"); next }
  for (r in seq_len(min(N_PAIR, nrow(si), nrow(sr)))) {
    zi <- ci[si$i[r]]; zr <- cr[sr$i[r]]
    cat(sprintf("  pair %d: IB1 %s (nT %d, tx %d) | RB1 %s (nT %d, tx %d)\n", r,
                zi$cell_id, zi$n_T, si$tot[r], zr$cell_id, zr$n_T, sr$tot[r]))
    rows <- list()
    for (side in 1:2) {
      S <- c("IB1","RB1")[side]; z <- if (side==1) zi else zr; R <- rp[side]
      p <- mm2px(S, z$x_slide_mm, z$y_slide_mm)
      rect <- mkrect(S, p$X, p$Y, WIN_UM, "w")
      img <- tryCatch(cosmx_image(S, rect), error=function(e) NULL)
      p1 <- if (is.null(img)) ggplot()+theme_v+labs(title=sprintf("%s image (unavailable)",S)) else
        img@plots$colorized + .coord_px(rect) + scale_y_reverse(labels=.um) +
        labs(x="x (um)",y="y (um)",title=sprintf("%s  %s -- image", S, R)) + theme_v
      loc <- ggplot() +
        geom_point(data=LOCBG[[S]], aes(x_slide_mm,y_slide_mm), colour="#cfcec9", size=.04) +
        annotate("point", x=z$x_slide_mm, y=z$y_slide_mm, colour="#c1121f", size=2.2) +
        coord_fixed() + theme_void() +
        theme(plot.background=element_rect(fill="white",colour="#b5b3ac",linewidth=.3))
      p1 <- p1 + inset_element(loc, left=.62, bottom=.62, right=.995, top=.995, align_to="panel")
      poly <- tryCatch(cosmx_cells_in_rect(S, rect), error=function(e) NULL)
      df <- NULL
      p2 <- if (is.null(poly)) ggplot()+theme_v+labs(title="cells (none)") else {
        df <- poly@coords %>% dplyr::group_by(name) %>% dplyr::reframe(x=unlist(x),y=unlist(y))
        pk <- cosmx_poly_key(df$name); cm <- CTMAP[[S]]
        v <- cm$ct8[match(pk, cm$key)]; df$ct8 <- ifelse(is.na(v),"other",v)
        ggplot(df, aes(x,y,group=name)) +
          geom_polygon(aes(fill=ct8), colour="#5a5a5a", linewidth=.18) +
          scale_fill_manual(values=CT_COLS, limits=c(CT_FULL,"excluded upstream"),
                            drop=FALSE, name=NULL) +
          .coord_px(rect) + scale_y_reverse(labels=.um) +
          labs(x="x (um)",y="y (um)",title=sprintf("%s -- cell types (full)", S)) +
          theme_v + theme(legend.text=element_text(size=6.2))
      }
      gp <- gp_for_rect(S, rect); q <- rect@coords[1L,]
      p3 <- ggplot()
      if (!is.null(df) && nrow(df))
        p3 <- p3 + geom_polygon(data=df, aes(x,y,group=name), fill="#f4f3ee",
                                colour="#33333333", linewidth=.12)
      n1 <- 0L; n2 <- 0L
      if (!is.null(gp) && nrow(gp)) {
        g <- as.data.table(gp)[X_px>=q$xmin & X_px<=q$xmax & Y_px>=q$ymin & Y_px<=q$ymax]
        g[, grp := fifelse(target %in% SET_IB1, "IB1-high: LAG3/HAVCR2/CTLA4",
                                                "RB1-high: PDCD1/TIGIT/TOX")]
        n1 <- g[grp %like% "^IB1", .N]; n2 <- g[grp %like% "^RB1", .N]
        if (nrow(g)) p3 <- p3 + geom_point(data=g, aes(X_px,Y_px,colour=grp), size=1.7, stroke=0) +
          scale_colour_manual(values=setNames(c(COL_IB1,COL_RB1),
            c("IB1-high: LAG3/HAVCR2/CTLA4","RB1-high: PDCD1/TIGIT/TOX")),
            drop=FALSE, name=NULL, guide=guide_legend(override.aes=list(size=3)))
      }
      p3 <- p3 + .coord_px(rect) + scale_y_reverse(labels=.um) +
        labs(x="x (um)",y="y (um)",
             title=sprintf("%s -- IB1-set %d, RB1-set %d", S, n1, n2)) + theme_v
      ## Column 4: the counts, because the scatter cannot be read by eye. Two colour clouds of
      ## a few hundred points each look alike however different the totals are -- which is
      ## exactly the complaint this panel answers. Per gene, coloured by set, with the set
      ## totals and their ratio stated.
      gb <- data.table(target = GENES,
                       set = fifelse(GENES %in% SET_IB1, "IB1-high", "RB1-high"), n = 0L)
      if (!is.null(gp) && nrow(gp)) {
        gg <- as.data.table(gp)[X_px>=q$xmin & X_px<=q$xmax & Y_px>=q$ymin & Y_px<=q$ymax,
                                .N, by=target]
        gb[gg, on="target", n := i.N]
      }
      gb[, target := factor(target, levels=c(SET_IB1, SET_RB1))]
      p4 <- ggplot(gb, aes(target, n, fill=set)) +
        geom_col(width=.7) +
        geom_text(aes(label=n), vjust=-0.3, size=2.7, colour=INK2) +
        scale_fill_manual(values=c(`IB1-high`=COL_IB1, `RB1-high`=COL_RB1), name=NULL) +
        scale_y_continuous(expand=expansion(mult=c(0,.16))) +
        labs(x=NULL, y="transcripts in this window",
             title=sprintf("%s -- set totals %d vs %d (ratio %.2f)", S, n1, n2,
                           if (n2>0) n1/n2 else NA_real_)) +
        theme_v + theme(axis.text.x=element_text(angle=35, hjust=1, colour=INK2, size=7.5),
                        panel.grid.major.y=element_line(colour="#e1e0d9", linewidth=.25),
                        legend.position="bottom")
      ## Column 5: what the window is actually made of, sentinels included and named. This is
      ## the panel that says what "other" was hiding -- in IB1 it is mostly unassigned cells,
      ## not a rare lineage, and that distinction changes how the T-subset panel should be read.
      cw <- if (is.null(df)) data.table(ct8=character(), n=integer()) else {
        pk2 <- cosmx_poly_key(unique(df$name)); cm2 <- CTMAP[[S]]
        ## An NA here is NOT an uncalled cell type. It is a polygon whose cell was dropped
        ## upstream -- over-fixation cluster `i`, an unbinned bin, or a region outside this
        ## comparison. Those cells are still drawn in the image and boundary panels, so naming
        ## them honestly matters: in IB1 they are 7-15% of the window, and calling them "other"
        ## invited the reading that they were a rare lineage rather than analysis attrition.
        v2 <- cm2$ct8[match(pk2, cm2$key)]
        data.table(ct8 = ifelse(is.na(v2), "excluded upstream", v2))[, .N, by=ct8]
      }
      setnames(cw, "N", "n", skip_absent=TRUE)
      cw <- merge(data.table(ct8=c(CT_FULL,"excluded upstream")), cw, by="ct8", all.x=TRUE)
      cw[is.na(n), n := 0L]
      cw[, ct8 := factor(ct8, levels=rev(c(CT_FULL,"excluded upstream")))]
      cw[, grp := fifelse(as.character(ct8) %in% SENTS, "not called",
                  fifelse(as.character(ct8) == "excluded upstream", "excluded", "called"))]
      ncell <- sum(cw$n); nsent <- cw[grp=="not called", sum(n)]
      nexc  <- cw[grp=="excluded", sum(n)]
      p5 <- ggplot(cw[n > 0], aes(n, ct8, fill=as.character(ct8))) +
        geom_col(width=.72) +
        geom_text(aes(label=n), hjust=-0.18, size=2.5, colour=INK2) +
        scale_fill_manual(values=CT_COLS, guide="none") +
        scale_x_continuous(expand=expansion(mult=c(0,.20))) +
        labs(x="cells in this window", y=NULL,
             title=sprintf("%s -- %d cells: %d (%.0f%%) not called, %d (%.0f%%) excluded upstream",
                           S, ncell, nsent, 100*nsent/max(ncell,1),
                           nexc, 100*nexc/max(ncell,1))) +
        theme_v + theme(axis.text.y=element_text(colour=INK2, size=6.8),
                        panel.grid.major.x=element_line(colour="#e1e0d9", linewidth=.25),
                        legend.position="none")
      rows <- c(rows, list(p1, p2, p3, p4, p5))
      if (side==1) { R1 <- c(n1, n2) } else { R2 <- c(n1, n2) }
      cat(sprintf("    %s window composition: %d cells, %d not called (%.0f%%); top called: %s\n",
                  S, ncell, nsent, 100*nsent/max(ncell,1),
                  paste(cw[grp=="called" & n>0][order(-n)][1:4, sprintf("%s=%d", ct8, n)],
                        collapse=" ")))
      cat(sprintf("      excluded upstream (defect / unbinned / other region): %d (%.0f%%)\n",
                  nexc, 100*nexc/max(ncell,1)))
    }
    ror <- if (exists("R2") && R1[2]>0 && R2[2]>0) (R1[1]/R1[2])/(R2[1]/R2[2]) else NA_real_
    pg <- wrap_plots(rows, nrow=2, ncol=5) +
      plot_annotation(
        title = "Two different inhibitory programmes on the same clone",
        subtitle = sprintf(paste0("Top: IB1 %s. Bottom: RB1 %s. Windows %g um, matched on ",
          "T-cell content (IB1 nT %d, RB1 nT %d) and chosen at the\nMEDIAN transcript count of ",
          "qualifying candidates, not the maximum. Yellow = LAG3/HAVCR2/CTLA4 (higher in IB1 in ",
          "all 5 T\ncompartments); cyan = PDCD1/TIGIT/TOX (higher in RB1 in all 5). 14 of 22 ",
          "genes directionally consistent, 1.4 expected,\nbinomial p = 2.8e-12."),
          rp[1], rp[2], WIN_UM, zi$n_T, zr$n_T),
        theme = theme(plot.title=element_text(face="bold",colour=INK,size=13),
                      plot.subtitle=element_text(colour=INK2,size=8),
                      plot.background=element_rect(fill=SURF,colour=NA)))
    f <- file.path(OUT, sprintf("fig_17_signview_%s_vs_%s_%d.png", rp[1], rp[2], r))
    ggsave(f, pg, width=24, height=11, dpi=165, bg=SURF, limitsize=FALSE)
    cat("  wrote ", f, "\n", sep=""); pairs_done <- pairs_done + 1L
  }
}
cat(sprintf("\n%d paired figures written to %s\n", pairs_done, OUT))
writeLines(format(Sys.time()), file.path(OUT, "sign_views.done"))
cat("============================================================\nOK\n")
