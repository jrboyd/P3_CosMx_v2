## ============================================================================
## v1 REFERENCE ONLY -- DO NOT RUN. Byte-identical copy (below this banner) of
##   /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx/11b_HRS_and_FL_image_validation.R
## It cannot run here: its three source() paths and its 11_bcell_compartments_v5/
## inputs exist only in the v1 project, and it needs metadata columns v2 does not
## have (compartment, celltype, kappa_dom, cd30_counts, ...).
## The helpers it used now live in scripts/lib/cosmx_tiff.R (renamed with a
## cosmx_ prefix); scripts/06_tiff_demo.R is the worked example.
## Kept as the provenance record for the v1 validation_images/{HRS,FL}/*.png.
## ============================================================================

# fetch views from Tiff around FL enriched domains and putative HRS cells


# =============================================================================
# CosMx 6K — composite classic Hodgkin + follicular lymphoma
# B/HRS compartment: malignant-cell calling (FL + HRS)
#
# Case: same patient, PRIMARY (IB1, composite cHL + FL) vs RECURRENCE (RB1, FL).
# This object is the weakly-filtered B/HRS compartment: a cell was kept if ANY
# B marker (MS4A1/CD19/CD79A/CD79B) OR TNFRSF8 (CD30) was detected. It is
# therefore high-recall / low-precision — precision is recovered here.
#
# KEY FINDINGS driving this script (established from full-scale diagnostics):
#   * RB1 is FL throughout: every cluster is kappa-restricted (log2 k:l ~3-4.7,
#     i.e. 8:1-25:1) EXCEPT cluster 11 (~1.76 log2 = polyclonal reactive B).
#   * IB1 is polyclonal at the cluster level: NO coarse cluster exceeds ~1.8
#     log2. The FL clone is a sparse minority that only surfaces at the
#     sub-cluster / per-cell level. Cluster-based FL calling therefore CANNOT
#     work in IB1 and must be sample-aware + higher-resolution.
#   * BCL2 counts cap at ~5 in IB1 but reach ~18 in RB1; BCL2 magnitude is not a
#     reliable FL discriminator in IB1, so restriction (a ratio) is the anchor,
#     with BCL2 used only as supporting evidence.
#   * Light-chain restriction is POPULATION-level. Per-cell double-positivity is
#     spillover/threshold noise, so we use count-based DOMINANCE, never ">0".
#
# n = 1: descriptive / hypothesis-generating only. HRS calls MUST be confirmed
# on tissue (rare, large, multinucleated, segmentation-prone).
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(SpatialExperiment)
  library(hoodscanR); library(dplyr); library(tidyr); library(ggplot2)
})

library(patchwork)

source("assemble_cosmx/16_cosmx_tiff_mapping.R")
source("20_cell_boundaries.R")
source("21_gene_points.R")



# ------------------------------------------------------------------ CONFIG ----
OUT_DIR      <- "11_bcell_compartments_v5"
PRIMARY_ID   <- "IB1"
RECURR_ID    <- "RB1"

seu = readRDS(file.path(OUT_DIR, "seurat_labeled.rds"))

meta_df.full = as.data.frame(seu@meta.data)
meta_df.full$cell_id = rownames(meta_df.full)

meta_df.full$NucArea
meta_df = meta_df.full %>% select(cell_id, celltype, sample_id, bcell_clusters, bcell_clusters_high, joint_cluster_high, lineage, celltype_coarse, kappa_dom, lambda_dom, compartment, gc_pos, bcl2_count, cd30_counts)
area_df = meta_df.full %>% select(cell_id, Area, Area.um2, NucArea)

xy_df.l = lapply(seu@images, function(img){
  xy_df = img$centroids@coords %>% as.data.frame %>% tibble
  xy_df$cell_id = img$centroid@cells
  xy_df
})

xy_df = xy_df.l %>% bind_rows(.id = "sample_id")
anno_df = merge(meta_df, xy_df, by = c("cell_id", "sample_id"))

anno_df$celltype %>% table
anno_df.hrs = anno_df %>% filter(compartment == "HRS_malignant")

area_df = area_df %>% pivot_longer(-cell_id)
area_df = area_df %>% mutate(is_HRS = cell_id %in% anno_df.hrs$cell_id)
ggplot(area_df, aes(x = value, fill = is_HRS)) +
  geom_histogram(color = "black") +
  facet_grid(is_HRS~name, scales = "free")

ggplot(area_df, aes(x = is_HRS, y = value, fill = is_HRS)) +
  geom_boxplot(color = "black") +
  facet_wrap(~name, scales = "free")

anno_df.hrs = anno_df.hrs %>% arrange(-cd30_counts)

ggplot(anno_df.hrs, aes(x = cd30_counts, fill = sample_id)) +
  facet_wrap(~sample_id) +
  geom_histogram(color = "black")

library(TiffPlotR)
mm_view_size = .1

mm_rect = with(anno_df.hrs, {
  TiffRect(x - mm_view_size/2, x + mm_view_size/2, y - mm_view_size/2, y+mm_view_size/2, name = cell_id)
})

mm_rect



cosmx_tiffs = getCosMxTiffPaths()
channel_names = getCosMxChannelNames()
channel_colors = getCosMxChannelColors()
global_precalc_max = getCosMxPrecalcMax()

# snapshots sent to UW spatial core to figure out channel G
if(FALSE){
  test_view = TiffRect(40e3, 45e3, 45e3, 50e3)
  test_view = test_view %>% shape_shift(-1e3, 2e3)
  fetchTiffData(cosmx_tiffs$IB1, precalc_max = global_precalc_max$IB1, channel_names = channel_names, rect = test_view)

  test_view = test_view %>% shape_resize_mult(.2)
  test_img.ib1 = fetchTiffData(cosmx_tiffs$IB1, precalc_max = global_precalc_max$IB1, channel_names = channel_names, rect = test_view)


  img_dat = fetchTiffData(cosmx_tiffs$RB1, precalc_max = global_precalc_max$RB1, channel_names = channel_names)

  test_view = TiffRect(40e3, 45e3, 45e3, 50e3)
  test_view = test_view %>% shape_shift(5e3, -5e3)

  img_dat %>% shape_annotate(test_view)

  fetchTiffData(cosmx_tiffs$RB1, precalc_max = global_precalc_max$RB1, channel_names = channel_names, rect = test_view)

  test_view = test_view %>% shape_resize_mult(.2)
  test_img.rb1 = fetchTiffData(cosmx_tiffs$RB1, precalc_max = global_precalc_max$RB1, channel_names = channel_names, rect = test_view)
  fetchTiffData(cosmx_tiffs$RB1, channel_names = channel_names, rect = test_view)

  library(patchwork)
  test_img.ib1@plots$normalized + labs(title = "IB1") +
    test_img.rb1@plots$normalized + labs(title = "RB1")

  ggsave("test_channel_names.png", width = 16, height = 7, dpi = 150)
}
aff_ib1 = fit_mm_to_px(seu, "IB1")
aff_rb1 = fit_mm_to_px(seu, "RB1")

aff_list = list(
  IB1 = aff_ib1,
  RB1 = aff_rb1
)


anno_df.hrs.by_id = split(anno_df.hrs, anno_df.hrs$sample_id)

px_xy.ib1 = mm_to_tiff_px(anno_df.hrs.by_id$IB1 %>% select(x, y), aff_ib1)
px_xy.rb1 = mm_to_tiff_px(anno_df.hrs.by_id$RB1 %>% select(x, y), aff_rb1)

#expect perfect line
plot(px_xy.ib1$X_px, anno_df.hrs.by_id$IB1$x)
plot(px_xy.rb1$X_px, anno_df.hrs.by_id$RB1$x)

plot(px_xy.ib1$Y_px, anno_df.hrs.by_id$IB1$y)
plot(px_xy.rb1$Y_px, anno_df.hrs.by_id$RB1$y)


anno_df.hrs.by_id$IB1 = cbind(anno_df.hrs.by_id$IB1, px_xy.ib1)
anno_df.hrs.by_id$RB1 = cbind(anno_df.hrs.by_id$RB1, px_xy.rb1)

anno_df.hrs_px = rbind(anno_df.hrs.by_id$IB1, anno_df.hrs.by_id$RB1) %>% tibble
anno_df.hrs_px = anno_df.hrs_px %>% arrange(-cd30_counts)

anno_df.hrs_px %>% select(cell_id, X_px, Y_px, cd30_counts)

px_view_size = 500

px_rect.all = with(anno_df.hrs_px, {
  TiffRect(X_px - px_view_size/2, X_px + px_view_size/2, Y_px - px_view_size/2, Y_px +px_view_size/2, name = cell_id)
})
px_rect.all@coords %>% head

i = 1

grps = sub("_.+", "", px_rect.all$name)
px_rect.ib1 = px_rect.all[grps == "IB1"]
px_rect.rb1 = px_rect.all[grps == "RB1"]

celltype_cols = seqsetvis::safeBrew(seu$celltype %>% unique %>% sort)

plot(seq_along(celltype_cols), rep(1, length(celltype_cols)), col = celltype_cols, pch = 16, cex = 3)
text(seq_along(celltype_cols), rep(1.1, length(celltype_cols)), labels = names(celltype_cols))

celltype_cols["B cell"] = "dodgerblue2"
celltype_cols["Mast cell"] =celltype_cols["NK cell"]
celltype_cols["NK cell"] = "firebrick2"

plot(seq_along(celltype_cols), rep(1, length(celltype_cols)), col = celltype_cols, pch = 16, cex = 3)
text(seq_along(celltype_cols), rep(1.1, length(celltype_cols)), labels = names(celltype_cols))


.plot_hrs_i = function(px_rect, i){
  message(i)
  i_rect = px_rect[i]
  i_cell_id = px_rect$name[i]
  i_anno_df = anno_df.hrs_px %>% filter(cell_id == i_cell_id)
  i_x = i_anno_df$X_px[i]
  i_y = i_anno_df$Y_px[i]
  i_sample_id = i_rect$name %>% sub("_.+", "", .)
  tf = cosmx_tiffs[[i_sample_id]]



  i_precalc_max = global_precalc_max[[i_sample_id]]
  i_precalc_max$channel = as.numeric(i_precalc_max$channel)


  img_data = fetchTiffData.colorized(
    tf,
    i_rect,
    channel_names = channel_names,
    precalc_max = i_precalc_max,
    channel_colors = channel_colors
  )

  cells_in_rect(i_rect, i_sample_id, return_fovs = TRUE)   # uncoloured polygons
  cell_poly = cells_in_rect(i_rect, i_sample_id, color = FALSE)   # uncoloured polygons

  palette <- celltype_cols



  cell_poly = TiffPlotR::assign_adjacent_colors(cell_poly, fill = TRUE, border = TRUE, palette = palette)

  theme_set(theme(panel.background = element_blank(), panel.grid = element_blank(), legend.position = "bottom"))

  img_data = img_data %>% add_um_scale_bar(label_offset_fraction = .02)
  img_data@plots$colorized = img_data@plots$colorized + annotate("point", x = i_x, y = i_y, color = "white", shape = 23, size = 3)
  p_img = img_data@plots$colorized + annotate("point", x = i_x, y = i_y, color = "white", shape = 23, size = 3)
  p_cells = p_img %>% shape_annotate(cell_poly, fill = "#00000000")


  hrs_genes = c("TNFRSF8")

  gene_points_in_rect(i_rect, i_sample_id, hrs_genes, endogenous_only = FALSE)

  gp <- gene_points_in_rect(i_rect, i_sample_id, hrs_genes,
                            palette = c(TNFRSF8 = "cyan", FUT4 = "yellow"))

  p_genes <- add_gene_points(p_img, gp, size = 1.6)
  # p_genes


  p_cells = p_cells + guides(color = "none")
  # cowplot::plot_grid(p_img, p_cells, p_genes, nrow = 1)

  pg = p_img + p_cells + p_genes
  pg
}

px_rect = px_rect.rb1
i = 1

fig_dir = file.path(OUT_DIR, "validation_images")
dir.create(fig_dir)

if(FALSE){
  for(px_rect in list(px_rect.ib1, px_rect.rb1)){
    for(i in seq(1:20)){
      i_rect = px_rect[i]
      pg = .plot_hrs_i(px_rect, i)
      dir.create(paste0(fig_dir, "/HRS"))
      # ggsave(paste0("quick_figs/validate_HRS_", i, "_", i_rect$name, ".pdf"), plot = pg, width = 12, height = 5.5)
      ggsave(paste0(fig_dir, "/HRS/validate_HRS_", i, "_", i_rect$name, ".png"), plot = pg, width = 12, height = 5.5)
    }
  }
}else{
  warning("skip HRS plots")
}

seu@meta.data$compartment %>% table
bc = LayerData(seu, "counts")

.goi_counts = function(M, goi){
  if(!goi %in% rownames(M)) stop(goi, " missing")
  bc[goi, ] %>% as.matrix %>% as.vector
}


k <- .goi_counts(bc, "IGKC")
l <- .goi_counts(bc, "IGLC1/2")
pseudo = 2
kl_lgr = log2((k + pseudo) / (l + pseudo))
hist(kl_lgr, breaks = seq(-10, 10, by= .5)+.25)

seu$kl_lgr = kl_lgr

idx = sample(ncol(seu), 5e3)

plot(k[idx], seu@meta.data$x_FOV_px[idx])
plot(l[idx], seu@meta.data$x_FOV_px[idx])
plot(kl_lgr[idx], seu@meta.data$x_FOV_px[idx])

plot(k[idx], seu@meta.data$y_FOV_px[idx])
plot(l[idx], seu@meta.data$y_FOV_px[idx])
plot(kl_lgr[idx], seu@meta.data$y_FOV_px[idx])

df = data.frame(l = l, k = k, kl_lgr = kl_lgr, x_fov = seu@meta.data$x_FOV_px, y_fov = seu@meta.data$y_FOV_px, sample_id = seu@meta.data$sample_id)
ggplot(df[idx,], aes(x = x_fov, y = l)) +
  geom_jitter(width = 15, height = 5)

ggplot(df[idx,], aes(x = x_fov, y = k)) +
  geom_jitter(width = 15, height = 5)

ggplot(df[idx,], aes(x = x_fov, y = kl_lgr)) +
  geom_jitter(width = 15, height = .8)


ggplot(df[idx,], aes(x = y_fov, y = l)) +
  geom_jitter(width = 15, height = 5)

ggplot(df[idx,], aes(x = y_fov, y = k)) +
  geom_jitter(width = 15, height = 5)

ggplot(df[idx,], aes(x = y_fov, y = kl_lgr)) +
  geom_jitter(width = 15, height = .8)

df$y_fov %>% range

df = df %>% mutate(y_bin = ceiling(y_fov/50), x_bin = ceiling(x_fov/50))
df.agg = df %>% group_by(sample_id, y_bin, x_bin) %>% summarise(l = sum(l), k = sum(k), kl_lgr = mean(kl_lgr))

df.agg

all_plots = list()

ggplot(df.agg %>% filter(sample_id == "IB1"), aes(x = x_bin, y = y_bin, fill = l)) +
  geom_raster() +
  facet_wrap(~sample_id)

ggplot(df.agg %>% filter(sample_id == "IB1"), aes(x = x_bin, y = y_bin, fill = k)) +
  geom_raster() +
  facet_wrap(~sample_id)

ggplot(df.agg %>% filter(sample_id == "IB1"), aes(x = x_bin, y = y_bin, fill = kl_lgr)) +
  geom_raster() +
  facet_wrap(~sample_id)

ggplot(df.agg %>% filter(sample_id == "RB1"), aes(x = x_bin, y = y_bin, fill = l)) +
  geom_raster() +
  facet_wrap(~sample_id)

ggplot(df.agg %>% filter(sample_id == "RB1"), aes(x = x_bin, y = y_bin, fill = k)) +
  scale_x_continuous(labels = function(x)x*50) +
  scale_y_continuous(labels = function(x)x*50) +
  geom_raster() +
  facet_wrap(~sample_id)

ggplot(df.agg %>% filter(sample_id == "RB1"), aes(x = x_bin, y = y_bin, fill = kl_lgr)) +
  geom_raster() +
  facet_wrap(~sample_id)

#### FL enriched regions ####
region_files = dir("11_bcell_compartments_v5/", pattern = "region.+csv$", full.names = TRUE)
names(region_files) = region_files %>% basename %>% str_split_i("[_\\.]", 4)
region_files = as.list(region_files)

rf = region_files$IB1

zoom_factor = 1



.plot_sample_fl_i = function(rf, i, zoom_factor = 1, plot_indexes = c(1,2,3,4), point_cex = 1, show_rect = FALSE){
  i_sample_id = strsplit(basename(rf), "[_\\.]")[[1]][4]
  message(i_sample_id)
  i_tiff = cosmx_tiffs[[i_sample_id]]
  regions_df = readr::read_csv(rf)

  message(i)
  region_info = regions_df[i,] %>% as.list
  mm_xy_list = list(
    xy_min = c(region_info$x_min_mm, region_info$y_min_mm),
    xy_max = c(region_info$x_max_mm, region_info$y_max_mm)
  )

  i_aff = aff_list[[i_sample_id]]

  px_xy_list = lapply(mm_xy_list, function(x){
    mm_to_tiff_px(data.frame(x = x[1], y = x[2]), i_aff)
  })
  px_rect_region = TiffRect(px_xy_list$xy_max$X_px, px_xy_list$xy_min$X_px, px_xy_list$xy_max$Y_px, px_xy_list$xy_min$Y_px)

  # resize to > of 2x shortest side or 1.2x the longest side
  dx = px_rect_region$xmax - px_rect_region$xmin
  dy = px_rect_region$ymax - px_rect_region$ymin
  side_long = max(dx, dy)
  side_short = min(dx, dy)
  w = max(2*side_short, 1.2*side_long)
  w = zoom_factor*w

  px_rect = shape_resize_abs(px_rect_region, width = w, height = w)

  i_precalc_max = global_precalc_max[[i_sample_id]]
  i_precalc_max$channel = as.numeric(i_precalc_max$channel)

  if(1 %in% plot_indexes){
    img_data = fetchTiffData.colorized(i_tiff, px_rect, channel_names = channel_names, precalc_max = i_precalc_max,
                                       channel_colors = channel_colors)
    img_data = img_data %>% add_um_scale_bar(pad_fraction = .04)

    p_img = img_data@plots$colorized
    if(show_rect)
      p_img = p_img %>% shape_annotate(px_rect_region, color = "black", linewidth = 2)
  }else{
    p_img = NULL
  }

  if(any(plot_indexes %in% c(2,3,4))){
    poly <- cells_in_rect(px_rect, i_sample_id, color = FALSE)   # uncoloured polygons
    poly <- color_cells_by_metadata(poly, seu, "celltype", border = TRUE, palette = celltype_cols)     # categorical clusters

    poly_kl <- color_cells_by_metadata(poly, seu, "kl_lgr", border = TRUE)     # categorical clusters
    poly_kl_df = poly_kl@coords %>% group_by(name) %>% reframe(x = unlist(x), y = unlist(y))

    poly_kl_df = merge(poly_kl_df, poly_kl@meta %>% select(name, value), by = "name")


    igl_genes = c("IGKC", "IGLC1/2")
    if(2 %in% plot_indexes){
      gp = gene_points_in_rect(px_rect, i_sample_id, igl_genes, endogenous_only = FALSE, palette = c("red", "blue"))

      small_side = min(abs(px_rect$xmin - px_rect$xmax), abs(px_rect$ymin - px_rect$ymax))

      psize = 500/small_side*point_cex

      p_genes = ggplot() +
        coord_fixed(xlim = c(px_rect$xmin, px_rect$xmax), ylim = c(px_rect$ymin, px_rect$ymax)) +
        scale_y_reverse() +
        labs(x = "px", y = "px")


      p_genes <- add_gene_points(p_genes, gp, size = psize)
      p_genes = p_genes + geom_polygon(data= poly_kl_df, aes(x = x, y = y, group = name), color = "#33333333", fill = NA)

      if(show_rect)
        p_genes = p_genes %>% shape_annotate(px_rect_region, color = "black", linewidth = 2)
    }else{
      p_genes = NULL
    }
    celltypes_df = poly@meta %>% select(value, fill) %>% unique

    p_cell_types = ggplot(data = celltypes_df,
                          aes(fill = value,
                              x = value %>% factor %>% as.numeric,
                              y = value %>% factor %>% as.numeric
                          )) +
      # geom_point(shape = 21) +
      geom_tile() +
      scale_fill_manual(values =celltype_cols) +
      labs(fill = "cell type")


    p_cell_types = p_cell_types %>% shape_annotate(poly, color = "gray", alpha = 1)

    p_cell_types = p_cell_types +
      coord_fixed(xlim = c(px_rect$xmin, px_rect$xmax), ylim = c(px_rect$ymin, px_rect$ymax)) +
      labs(x = "px", y = "px") +
      scale_y_reverse()



    lims = poly_kl_df$value %>% range(na.rm = TRUE)
    lims = max(abs(lims))
    lims = c(-lims, lims)

    p_kl = ggplot(poly_kl_df, aes(x = x, y = y, group = name, fill = value)) +
      geom_polygon(color = "gray") +
      scale_fill_gradientn(limits = lims, colors = c("blue", "white", "red"), na.value = "gray95") +
      coord_fixed(xlim = c(px_rect$xmin, px_rect$xmax), ylim = c(px_rect$ymin, px_rect$ymax)) +
      scale_y_reverse() +
      labs(x = "px", y = "px", fill = "log2 ratio\nof kappa vs lambda")


    if(show_rect)
      p_cell_types = p_cell_types %>% shape_annotate(px_rect_region, color = "black", linewidth = 2)
    if(show_rect)
      p_kl = p_kl %>% shape_annotate(px_rect_region, color = "black", linewidth = 2)
  }

  if(all(plot_indexes == c(1,2,3,4))){
    pg = (p_img + p_genes) /
      (p_cell_types + p_kl)
  }else if(length(plot_indexes) == 1){
    pg = switch (LETTERS[plot_indexes],
                 A = p_img,
                 B = p_genes,
                 C = p_cell_types,
                 D = p_kl
    )
  }else{
    pg = do.call(plot_layout, c(list(p_img, p_genes, p_cell_types, p_kl)[plot_indexes], list(nrow = 1)))
  }

  pg
}

rf = region_files$RB1
rf = region_files$IB1
i = 17
for(rf in region_files){
  i_sample_id = strsplit(basename(rf), "[_\\.]")[[1]][4]
  regions_df = readr::read_csv(rf)
  for(i in seq(20)){
    region_info = regions_df[i,] %>% as.list
    region_info$x_center_mm
    region_info$y_center_mm

    pg = .plot_sample_fl_i(rf, i)

    dir.create(file.path(fig_dir, "FL"))
    out_f = file.path(fig_dir, "FL", paste0("validate_FL_", i_sample_id, "_", i, "_", region_info$region_id, ".png"))
    ggsave(out_f, plot = pg, width = 14, height = 12)
  }
}

# FL regions with high signal
# 2, 5, 10, 11, 12, 18, 19

.plot_sample_fl_i(region_files$IB1, 2, plot_indexes = 1)
.plot_sample_fl_i(region_files$IB1, 2, plot_indexes = 2)
.plot_sample_fl_i(region_files$IB1, 2, zoom_factor = 3, plot_indexes = 2)

.plot_sample_fl_i(region_files$IB1, 5, zoom_factor = 3, plot_indexes = 2)

.plot_sample_fl_i(region_files$IB1, 10, zoom_factor = 3, plot_indexes = 2)


.plot_sample_fl_i(region_files$IB1, 5)
.plot_sample_fl_i(region_files$IB1, 10)
.plot_sample_fl_i(region_files$IB1, 18)
