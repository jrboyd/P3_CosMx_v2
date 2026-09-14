## ============================================================================
## 10_script_hrs_views.R -- image views of the HRS objects called by
## scripts/10_script_hrs_call.R.
##
## WHY THIS EXISTS AND WHY IT CARRIES CONTROLS.
## scripts/01_import_cosmx_seurat.R:25-26: "HRS calls are heuristic and MUST be
## confirmed visually / with pathology". A picture of a giant multinucleated CD30+ cell
## ringed by T cells is the only thing that settles whether the call is real -- but a
## plausible-looking picture proves nothing on its own, so the same view type is
## rendered for two control sets:
##   HRS          top IB1 objects, spanning the n_segments range rather than just the
##                largest, so single- and multi-segment calls can be compared
##   RB1control   top-scoring RB1 objects. RB1 is the FL recurrence with no Hodgkin
##                component, so these SHOULD NOT look like HRS
##   subthresh    IB1 objects scoring just below the cut
##
## Panels follow scripts/07_script_lightchain_views.R: locator / image / cells /
## transcripts. The window is centred on the OBJECT, and every member segment is
## outlined -- if those outlines carve up one visibly continuous giant cell, the merge
## in section A2 of the call script is doing what it is supposed to.
##
## Usage: Rscript scripts/10_script_hrs_views.R [sample] [n_hrs] [n_control]
##        scripts/10_submit_hrs_views.sh
## ============================================================================

if (!file.exists("scripts/lib/cosmx_tiff.R"))
  stop("run from the project root", call. = FALSE)
source("scripts/lib/cosmx_tiff.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table); setDTthreads(1)
})

args      <- commandArgs(trailingOnly = TRUE)
SAMPLE_ID <- if (length(args) >= 1) args[1] else "IB1"
N_HRS     <- if (length(args) >= 2) as.integer(args[2]) else 12L
N_CTRL    <- if (length(args) >= 3) as.integer(args[3]) else 6L

IN_DIR  <- "output/10_hrs"
OUT_DIR <- file.path("output/10_hrs_views", SAMPLE_ID)
FLAT    <- "/netfiles/volaric_research/DLBCL/CosMx/export1/flatFiles/%s/%s_metadata_file.csv.gz"
OBJ_F   <- file.path(IN_DIR, "hrs_objects.csv.gz")
SEG_F   <- file.path(IN_DIR, "hrs_segments.csv.gz")
for (f in c(OBJ_F, SEG_F)) if (!file.exists(f))
  stop("missing ", f, " -- run scripts/10_submit_hrs_call.sh first", call. = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAD       <- 1.30
MIN_SIDE  <- 700L      # px, ~84 um -- the legibility floor from 07_script_lightchain_views.R:162
ZOOMS     <- c(tight = 1, zoom3x = 3)
LOC_N     <- 60000L
PT_SIZE   <- 1.5

## Transcripts. COLOUR ENCODES FAMILY, NOT GENE -- three categories, not seven.
##
## Panel 3 is a scatter, so every pair of colours is simultaneously adjacent on screen
## and the palette has to clear the ALL-PAIRS bar, not the easier adjacent-pair one.
## The project palette validates all-pairs only for its first three slots (worst CVD
## dE 9.2, normal-vision dE 24.0, on this exact #fcfcfb surface); past three it fails.
## Seven distinguishable gene hues were never available. The first version of this
## panel tried anyway and collapsed: CCL17/IRF4/TRAF1 all landed on one pale cyan,
## TNFRSF8 sat a shade away from them, CD3E and CD4 were two near-identical ambers,
## and PDCD1 came out pink while the subtitle claimed the T-cell genes were amber.
##
## So: three validated hues by family, PDCD1 separated by SHAPE rather than colour
## (PD-1+ T cells rosetting an HRS cell is a specific cHL feature and worth keeping
## visible), and CD30 drawn last, larger, in slot 1 so the focal gene dominates.
GENE_FAMILY <- c(TNFRSF8 = "CD30",
                 CCL17 = "cHL", IRF4 = "cHL", TRAF1 = "cHL",
                 CD3E = "Tcell", CD4 = "Tcell", PDCD1 = "Tcell")
GENES     <- names(GENE_FAMILY)
FAM_LEVELS<- c("CD30", "cHL", "Tcell")
FAM_COLS  <- c(CD30 = "#2a78d6", cHL = "#eb6834", Tcell = "#1baf7a")  # slots 1/2/3
## Labels kept short: at 19in / 4 panels the first version clipped the third entry
## mid-word ("... = PDCD").
FAM_LABS  <- c(CD30  = "CD30 (TNFRSF8)",
               cHL   = "cHL (CCL17/IRF4/TRAF1)",
               Tcell = "T cell (CD3E/CD4, ▲PDCD1)")
FAM_SIZE  <- c(CD30 = 2.3, cHL = 1.2, Tcell = 1.2)
## cosmx_gene_points() wants a per-gene palette; it is only used to tag rows, since
## the layer below draws from `family` rather than from that column.
GENE_COLS <- FAM_COLS[GENE_FAMILY[GENES]]; names(GENE_COLS) <- GENES

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; AXIS<-"#c3c2b7"
## Panel 2 fills. Member and other-CD30+ are the SAME hue at two lightnesses -- they
## are two levels of one thing (CD30 positivity), not two unrelated categories, so a
## light/dark pair of slot 1 is the right encoding rather than a second hue.
CELL_FILL <- c(hrs = "#2a78d6", cd30 = "#a8c6ea", other = "#dedcd4")
theme_v <- theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid = element_blank(),
        axis.text = element_text(colour = MUTED, size = 7),
        axis.title = element_text(colour = INK2, size = 8),
        plot.title = element_text(colour = INK2, face = "bold", size = 9),
        legend.text = element_text(colour = INK2, size = 8),
        legend.title = element_blank(), legend.position = "bottom")

## Orientation. MUST be scale_x_reverse(), not reversed coord limits: ggplot sorts
## limits passed to coord_fixed, so coord_fixed(xlim = c(max, min)) silently does
## nothing. See scripts/07_script_lightchain_views.R:576-591.
.coord_px <- function(rect) {
  q <- rect@coords[1L, ]
  list(scale_x_reverse(),
       coord_fixed(xlim = c(q$xmin, q$xmax), ylim = c(q$ymin, q$ymax)))
}
## rect-clip gene points while PRESERVING attr("palette"), which cosmx_add_genes needs
## and plain `[` drops (07_script_lightchain_views.R:343)
.gp_subset <- function(gp, rect) {
  if (is.null(gp) || !nrow(gp)) return(NULL)
  q <- rect@coords[1L, ]
  pal <- attr(gp, "palette")
  out <- gp[gp$X_px >= q$xmin & gp$X_px <= q$xmax &
            gp$Y_px >= q$ymin & gp$Y_px <= q$ymax, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  attr(out, "palette") <- pal
  out
}

cat("============================================================\n")
cat("step 10 views -- HRS objects, ", SAMPLE_ID, "\n", sep = "")
cat("============================================================\n")

## ---- 1. selection -----------------------------------------------------------
OBJ <- fread(OBJ_F); SEG <- fread(SEG_F)
OBJ[, samp := as.character(samp)]; SEG[, samp := as.character(samp)]
cat(sprintf("objects: %d (%s)\n", nrow(OBJ),
            paste(OBJ[, .N, by = samp][, paste0(samp, "=", N)], collapse = " ")))

hrs <- OBJ[samp == SAMPLE_ID & is_hrs == TRUE][order(-obj_score)]
if (!nrow(hrs)) stop("no called HRS objects in ", SAMPLE_ID,
                     " -- nothing to render; read output/10_hrs/summary.csv", call. = FALSE)
## span the n_segments range rather than taking the top N by score alone: if every
## rendered object is a 5-segment giant, the figures cannot show what a single-segment
## call looks like, and that is exactly the comparison the merge needs
hrs[, seg_bin := pmin(n_segments, 4L)]
sel_hrs <- hrs[, head(.SD, ceiling(N_HRS / uniqueN(hrs$seg_bin))), by = seg_bin][
             order(-obj_score)][seq_len(min(N_HRS, .N))]
sel_hrs[, group := "HRS"]

sub <- OBJ[samp == SAMPLE_ID & is_hrs == FALSE][order(-obj_score)][seq_len(min(N_CTRL, .N))]
if (nrow(sub)) sub[, group := "subthresh"]
rb <- OBJ[samp != SAMPLE_ID][order(-obj_score)][seq_len(min(N_CTRL, .N))]
if (nrow(rb)) rb[, group := "RB1control"]

SEL <- rbindlist(list(sel_hrs[, !"seg_bin"], sub, rb), fill = TRUE)
SEL <- SEL[!is.na(obj_id)]
cat("\n=== selected for rendering\n")
print(SEL[, .(objects = .N, med_segments = as.numeric(median(n_segments)),
              med_area = round(median(area_sum), 1)), by = group])
fwrite(SEL, file.path(OUT_DIR, "selected_objects.csv"))

## ---- 2. exact canvas px for every cell --------------------------------------
## Route A: FOV-local px straight from the flat file. No Seurat object and no fitted
## affine, so no fit residual enters the geometry.
render_samples <- unique(SEL$samp)
FLATS <- lapply(render_samples, function(S) {
  f <- fread(sprintf(FLAT, S, S),
             select = c("cell", "fov", "CenterX_local_px", "CenterY_local_px"))
  px <- cosmx_local_to_px(S, f$fov, f$CenterX_local_px, f$CenterY_local_px)
  f[, `:=`(key = cosmx_seu_key(cell), X_px = px$X_px, Y_px = px$Y_px)]
  f
})
names(FLATS) <- render_samples

SEG[, key := cosmx_seu_key(cell_id)]
SEG <- rbindlist(lapply(render_samples, function(S)
  merge(SEG[samp == S], FLATS[[S]][, .(key, fov, X_px, Y_px)], by = "key", all.x = TRUE)))
n_lost <- SEG[obj_id %in% SEL$obj_id & is.na(X_px), .N]
if (n_lost) warning(n_lost, " selected-object segments did not join the flat file")

## ---- 3. per-sample locator background ---------------------------------------
PX2MM <- lapply(render_samples, function(S) {
  d <- SEG[samp == S & !is.na(X_px)]
  fx <- lm(x_slide_mm ~ X_px, data = d); fy <- lm(y_slide_mm ~ Y_px, data = d)
  r <- c(sd(resid(fx)), sd(resid(fy)))
  cat(sprintf("%s px->mm for the locator: resid SD %.4f / %.4f mm over %d cells\n",
              S, r[1], r[2], nrow(d)))
  if (max(r) > 0.01) warning(S, ": px->mm residual over 10 um; locator may be off")
  list(ax = unname(coef(fx)[2]), bx = unname(coef(fx)[1]),
       ay = unname(coef(fy)[2]), by = unname(coef(fy)[1]))
})
names(PX2MM) <- render_samples
set.seed(1)
LOCBG <- lapply(render_samples, function(S) {
  f <- FLATS[[S]]; b <- f[sample.int(.N, min(.N, LOC_N)), .(X_px, Y_px)]
  a <- PX2MM[[S]]
  b[, `:=`(x_mm = a$ax * X_px + a$bx, y_mm = a$ay * Y_px + a$by)]; b
})
names(LOCBG) <- render_samples

.locator_plot <- function(S, rect) {
  q <- rect@coords[1L, ]; a <- PX2MM[[S]]
  wx <- range(a$ax * c(q$xmin, q$xmax) + a$bx)
  wy <- range(a$ay * c(q$ymin, q$ymax) + a$by)
  ggplot() +
    geom_point(data = LOCBG[[S]], aes(x_mm, y_mm), colour = "#cfcec9", size = 0.04) +
    annotate("rect", xmin = wx[1], xmax = wx[2], ymin = wy[1], ymax = wy[2],
             colour = "#c1121f", fill = NA, linewidth = 0.4) +
    coord_fixed() +
    labs(title = sprintf("0. slide position, mm\n(%.2f, %.2f)", mean(wx), mean(wy))) +
    theme_void() +
    theme(plot.title = element_text(size = 8),
          plot.background = element_rect(fill = "white", colour = "#b5b3ac", linewidth = 0.3))
}

CANVAS <- lapply(render_samples, function(S) {
  pl <- cosmx_fov_placement(S)
  c(x = max(pl$x0) + COSMX_FOV_PX - 1L, y = max(pl$y0) + COSMX_FOV_PX - 1L)
})
names(CANVAS) <- render_samples

## ---- 4. panel composers ------------------------------------------------------
## Panel 2 as a standalone polygon plot rather than an overlay: a solid fill over the
## multichannel image washes both out (07_script_lightchain_views.R:501-503).
.cell_plot <- function(S, rect, member_keys, cd30_keys) {
  poly <- tryCatch(cosmx_cells_in_rect(S, rect), error = function(e) NULL)
  if (is.null(poly)) return(list(p = ggplot() + theme_v + labs(title = "2. cells (none)"),
                                 members = character(0)))
  df <- poly@coords %>% dplyr::group_by(name) %>%
        dplyr::reframe(x = unlist(x), y = unlist(y))
  k  <- cosmx_poly_key(df$name)
  df$fill <- ifelse(k %in% member_keys, CELL_FILL[["hrs"]],
             ifelse(k %in% cd30_keys,   CELL_FILL[["cd30"]], CELL_FILL[["other"]]))
  p <- ggplot(df, aes(x = x, y = y, group = name)) +
    geom_polygon(aes(fill = fill), colour = "black", linewidth = 0.3) +
    geom_polygon(fill = NA, colour = "white", linewidth = 0.12) +
    scale_fill_identity()
  ## every member segment outlined, so fragmentation of one giant cell is visible
  fo <- df[k %in% member_keys, , drop = FALSE]
  if (nrow(fo))
    p <- p + geom_polygon(data = fo, fill = NA, colour = "#ffffff", linewidth = 0.9)
  list(p = p + .coord_px(rect) + scale_y_reverse() + labs(x = "pixel", y = "pixel",
         title = "2. cells -- member (dark blue, white outline), other CD30+ (pale)") + theme_v,
       poly = poly, df = df)
}

.gene_plot <- function(gp, df, rect) {
  p <- ggplot()
  if (!is.null(df) && nrow(df))
    p <- p + geom_polygon(data = df, aes(x = x, y = y, group = name),
                          fill = "#f2f1ec", colour = NA)
  n_fam <- setNames(integer(length(FAM_LEVELS)), FAM_LEVELS)
  if (!is.null(gp) && nrow(gp)) {
    d <- as.data.table(gp)
    d[, family := factor(GENE_FAMILY[target], levels = FAM_LEVELS)]
    d <- d[!is.na(family)]
    d[, pd1 := target == "PDCD1"]
    n_fam[FAM_LEVELS] <- as.integer(table(factor(d$family, levels = FAM_LEVELS)))
    ## legend order comes from the factor levels; DRAW order from row order, so the
    ## focal gene is plotted last and is never buried under the supporting families
    ## Two layers rather than one with size mapped to `family`: mapping a discrete
    ## variable to size warns on every panel, and splitting the layers also guarantees
    ## the focal gene is drawn over the supporting families rather than relying on
    ## row order.
    p <- p +
      geom_point(data = d[family != "CD30"],
                 aes(X_px, Y_px, colour = family, shape = pd1),
                 size = FAM_SIZE[["cHL"]], stroke = 0) +
      geom_point(data = d[family == "CD30"],
                 aes(X_px, Y_px, colour = family, shape = pd1),
                 size = FAM_SIZE[["CD30"]], stroke = 0) +
      scale_colour_manual(values = FAM_COLS, labels = FAM_LABS, name = NULL,
                          drop = FALSE,
                          guide = guide_legend(override.aes = list(size = 3, shape = 16))) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none")
  }
  if (!is.null(df) && nrow(df))
    p <- p + geom_polygon(data = df, aes(x = x, y = y, group = name),
                          fill = NA, colour = "#33333355", linewidth = 0.15)
  ## Counts in the title are the relief for slot 3 sitting under 3:1 contrast on a
  ## light surface -- the reader gets the numbers without having to judge dot colour.
  p + .coord_px(rect) + scale_y_reverse() +
    labs(x = "pixel", y = "pixel",
         title = sprintf("3. transcripts -- CD30 %d, cHL %d, T cell %d",
                         n_fam[["CD30"]], n_fam[["cHL"]], n_fam[["Tcell"]])) +
    theme_v
}

## ---- 5. render, batched by FOV ----------------------------------------------
## cosmx_gene_points reads ~120 MB of transcript csv per intersecting FOV regardless
## of window size, so it is paid once per FOV and clipped per object. Cell boundaries
## are read per object via cosmx_cells_in_rect, which is cheaper and handles the FOV
## seam case correctly.
SEL <- merge(SEL, SEG[obj_id %in% SEL$obj_id,
                      .(fov = fov[1], n_seg_px = sum(!is.na(X_px))), by = obj_id],
             by = "obj_id", all.x = TRUE)
setorder(SEL, samp, fov, -obj_score)
cat(sprintf("\nrendering %d objects across %d FOVs\n", nrow(SEL), uniqueN(SEL[, .(samp, fov)])))

n_done <- 0L
for (S in render_samples) {
  for (FV in sort(unique(SEL[samp == S, fov]))) {
    if (is.na(FV)) next
    objs <- SEL[samp == S & fov == FV]
    pl   <- cosmx_fov_placement(S)
    row  <- pl[pl$fov == FV, ]
    if (!nrow(row)) { warning("FOV ", FV, " absent from the image set for ", S); next }
    fov_rect <- cosmx_px_rect(row$x0, row$x0 + COSMX_FOV_PX,
                              row$y0, row$y0 + COSMX_FOV_PX, name = sprintf("fov%d", FV))
    gp_fov <- tryCatch(cosmx_gene_points(S, fov_rect, GENES, palette = GENE_COLS),
                       error = function(e) { message("gene points failed for FOV ", FV,
                                                     ": ", conditionMessage(e)); NULL })
    cd30_keys <- SEG[samp == S & fov == FV & cd30_pos == TRUE, key]

    for (r in seq_len(nrow(objs))) {
      o    <- objs[r]
      mem  <- SEG[obj_id == o$obj_id & !is.na(X_px)]
      if (!nrow(mem)) next
      cx <- mean(range(mem$X_px)); cy <- mean(range(mem$Y_px))
      side <- max(diff(range(mem$X_px)), diff(range(mem$Y_px))) * PAD
      side <- max(side, MIN_SIDE)
      for (z in names(ZOOMS)) {
        s2   <- side * ZOOMS[[z]] / 2
        rect <- cosmx_px_rect(max(0, cx - s2), min(CANVAS[[S]][["x"]], cx + s2),
                              max(0, cy - s2), min(CANVAS[[S]][["y"]], cy + s2),
                              name = sprintf("%s_%s", o$obj_id, z))
        img <- tryCatch(cosmx_image(S, rect) %>% add_um_scale_bar(),
                        error = function(e) { message("image failed: ",
                                                      conditionMessage(e)); NULL })
        p_img <- if (is.null(img)) ggplot() + theme_v + labs(title = "1. image (unavailable)")
                 else img@plots$colorized + .coord_px(rect) + theme_v + labs(title = "1. image")
        cw  <- .cell_plot(S, rect, mem$key, cd30_keys)
        p_g <- .gene_plot(.gp_subset(gp_fov, rect), cw$df, rect)
        lab <- sprintf(paste0("%s %s | %s | %d segment(s), pooled area %.0f um2 | ",
                              "CD30 %d counts | score %.2f"),
                       S, o$group, o$obj_id, o$n_segments, o$area_sum, o$cd30_sum, o$obj_score)
        pg <- .locator_plot(S, rect) + p_img + cw$p + p_g +
          plot_layout(nrow = 1, widths = c(0.5, 1, 1, 1)) +
          plot_annotation(title = lab,
            subtitle = sprintf(paste0("%s -- %.0f px (%.0f um) window, FOV %d. White outline ",
                                      "= member segments of the object; if they carve up one ",
                                      "continuous giant cell, the merge is correct."),
                               z, 2*s2, 2*s2*COSMX_MM_PER_PX*1000, FV))
        out_f <- file.path(OUT_DIR, sprintf("%s_%s_%s.png", o$group, o$obj_id, z))
        ggsave(out_f, plot = pg, width = 19, height = 6, dpi = 140, limitsize = FALSE)
        cat("wrote ", out_f, "\n", sep = "")
        n_done <- n_done + 1L
      }
    }
    rm(gp_fov); invisible(gc())
  }
}

writeLines(format(Sys.time()), file.path(OUT_DIR, "views.done"))
cat("\n============================================================\n")
cat(sprintf("%d panels written to %s\n", n_done, OUT_DIR))
cat("Read these before believing the call. An HRS cell should be visibly LARGE and\n")
cat("multinucleated, carry CD30 transcripts, and sit in a ring of T cells. The\n")
cat("RB1control panels should look clearly different -- RB1 has no Hodgkin component,\n")
cat("so anything that looks like HRS there is a false positive in the method.\n")
cat("============================================================\nOK\n")
