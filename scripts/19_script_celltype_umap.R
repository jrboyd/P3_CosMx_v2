## ============================================================================
## 19_script_celltype_umap.R -- the step-19 cell types on the high-confidence embedding.
##
## Three panels on one embedding, because the project keeps TWO concepts of cell type and a
## figure that showed only one would hide the disagreement that matters:
##   1. celltype_cluster  guilt by association -- the cluster's name, inherited by every
##                        cell in it including the 43% that resolved no marker at all
##   2. lineage_argmax    direct per-cell marker evidence; cells that resolved nothing are
##                        drawn grey, which is most of the tissue and should look like it
##   3. by sample         IB1 and RB1 on the same axes, since the object is 62.6% RB1 and
##                        any structure that is really "which sample" must be visible as such
##
## `mixed` and `junk` clusters are drawn in grey and named as such. They carry no cell-type
## label on purpose (19d), and a figure that quietly coloured them would undo that.
##
## Usage: sbatch --mem=32G --time=1:00:00 scripts/sbatch/run_rscript.sh \
##          scripts/19_script_celltype_umap.R [high|low]
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); setDTthreads(2)
})
args <- commandArgs(trailingOnly = TRUE)
RES  <- if (length(args) >= 1) args[1] else "high"
OUT  <- "output/19_celltypes"; FIG <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
## hc_cells_final.csv.gz carries the refined labels (the Plasma cell cluster split into
## Plasma cell / Plasmablast / dropped); fall back to the pre-refinement table if absent.
F_IN <- file.path(OUT, "hc_cells_final.csv.gz")
if (!file.exists(F_IN)) F_IN <- file.path(OUT, sprintf("hc_cells_labelled_%s.csv.gz", RES))
cat("reading ", F_IN, "\n", sep="")
if (!file.exists(F_IN)) stop("missing ", F_IN, call. = FALSE)
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=170,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
PAL <- c("B cell" = "#184f95",
         ## Plasmablast: rose, deliberately adjacent to Plasma-cell magenta because it is the
         ## B-to-plasma intermediate, and distinct enough not to be read as the same label.
         "Plasmablast"                            = "#f2789f",
         "B cell, naive/mantle"                   = "#5b9bd5",
         "T cell"                                 = "#d03b3b",
         "NK cell"                                = "#7b3f9e",
         "Myeloid"                                = "#e08214",
         "Endothelial"                            = "#1b7837",
         "Stromal/FDC"                            = "#8c6d31",
         "Plasma cell"                            = "#c51b7d",
         "pDC"                                    = "#00857a",
         ## PALETTE CLASH, FLAGGED RATHER THAN SILENTLY RESOLVED. The project's fixed hue
         ## for Hodgkin is #0b3d7a, but on a cell-type figure that is indistinguishable from
         ## B cell (#184f95) -- and B is 54% of the object, so the niche would vanish into
         ## it. Near-black is used here instead, which also reads correctly: this is the one
         ## entry that is a REGION and not a cell type. The fixed hue still governs the
         ## disease figures in steps 12 and 16.
         "Hodgkin chemokine niche"                = "#151a2b",
         "unnamed (mixed)"                        = "#b9b8b2",
         "unnamed (junk)"                         = "#dcdcd6")
LPAL <- c(B="#184f95", T="#d03b3b", NK="#7b3f9e", myeloid="#e08214", endo="#1b7837",
          fibro="#8c6d31", PC="#c51b7d", `no marker call`="#dcdcd6")
th <- function(base=11) theme_minimal(base_size=base) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.major=element_line(colour=GRID,linewidth=.2), panel.grid.minor=element_blank(),
        axis.text=element_text(colour=MUTED,size=base-3), axis.title=element_text(colour=INK2,size=base-2),
        plot.title=element_text(colour=INK,face="bold",size=base+1),
        plot.subtitle=element_text(colour=INK2,size=base-2),
        plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
        strip.text=element_text(colour=INK2,face="bold"),
        legend.text=element_text(colour=INK2,size=base-2), legend.title=element_blank())

## fread returns an unwritten character field as "", not NA, so an is.na() test silently
## misses every unnamed cluster and they fall through to a bare NA level on the legend.
d <- fread(F_IN, na.strings = c("NA",""))
d[, ct := fifelse(!is.na(celltype_cluster), celltype_cluster,
                  paste0("unnamed (", cluster_status, ")"))]
if (d[, any(is.na(ct))]) stop("ct is NA for ", d[, sum(is.na(ct))], " cells", call. = FALSE)
d[, lin := fifelse(lineage_argmax %in% names(LPAL), lineage_argmax, "no marker call")]
n_by <- d[, .N, by = ct][order(-N)]
cat("=== cells per label ===\n"); print(n_by)
LEV <- c(intersect(names(PAL), n_by$ct))
d[, ct := factor(ct, levels = LEV)]
d[, lin := factor(lin, levels = names(LPAL))]
lg <- setNames(sprintf("%s  (%s, %.1f%%)", LEV, format(n_by[match(LEV, ct), N], big.mark=","),
                       100*n_by[match(LEV, ct), N]/nrow(d)), LEV)

## AXIS RANGE EXCLUDES THE JUNK CLUSTERS, though they are still drawn. 159 cells of the
## 688-cell junk cluster sit at U1 ~ -18.3 while 99.99% of the object lies between -4.42
## and +4.19 -- left alone they stretch x to 4.5x the width the data needs and squash two
## million cells into a fifth of the panel. Clipping is done with coord_cartesian, which
## zooms; scale limits would DROP the points instead of hiding them.
rl <- d[cluster_status != "junk"]
XL <- range(rl$U1); YL <- range(rl$U2)
pad <- 0.04 * c(diff(XL), diff(YL))
XL <- XL + c(-pad[1], pad[1]); YL <- YL + c(-pad[2], pad[2])
n_out <- d[U1 < XL[1] | U1 > XL[2] | U2 < YL[1] | U2 > YL[2], .N]
cat(sprintf("\naxis %.2f..%.2f x %.2f..%.2f -- %d cells outside it, all junk\n",
            XL[1], XL[2], YL[1], YL[2], n_out))
## coord_fixed, not coord_cartesian: a UMAP is a distance embedding, so unequal x and y
## scaling misrepresents it and makes panels of different widths incomparable. Limits are
## passed in ascending order -- ggplot SORTS the limits given to coord_fixed, so a reversed
## pair silently does nothing (see CLAUDE.md).
CO <- coord_fixed(ratio = 1, xlim = XL, ylim = YL, expand = FALSE)

set.seed(1L); s <- d[sample(.N, min(.N, 400000L))]
## grey first so the named types are never hidden under the unnamed ones. Panels 1 and 2
## colour by different columns, so each needs its own draw order -- reusing one ordering
## buries the signal in whichever panel it was not built for.
s1 <- s[order(!(ct  %like% "unnamed"))]
s2 <- s[order(lin != "no marker call")]
lab <- s1[!(ct %like% "unnamed"), .(U1 = median(U1), U2 = median(U2)), by = ct]

p1 <- ggplot(s1, aes(U1, U2, colour = ct)) + geom_point(size=.06, alpha=.45) +
  scale_colour_manual(values = PAL, labels = lg, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1), ncol=2)) +
  ggrepel::geom_text_repel(data = lab, aes(label = ct), colour = INK, size = 3,
                           fontface = "bold", show.legend = FALSE, seed = 1,
                           bg.color = SURF, bg.r = 0.12, min.segment.length = 0) +
  labs(title = sprintf("step 19 cell types -- %d high-confidence cells", nrow(d)),
       subtitle = "cluster label (guilt by association): every cell in a cluster inherits its name,\nincluding the 43% that resolve no marker of their own",
       x = "UMAP 1", y = "UMAP 2") + CO + th() + theme(legend.position = "bottom")

p2 <- ggplot(s2, aes(U1, U2, colour = lin)) + geom_point(size=.06, alpha=.45) +
  scale_colour_manual(values = LPAL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1), ncol=1)) +
  labs(title = "per-cell marker evidence (argmax)",
       subtitle = "the same cells, coloured by their OWN markers. Grey = no call:\nthis is what the clustering is filling in.",
       x = "UMAP 1", y = NULL) + CO + th() + theme(legend.position = "bottom")

p3 <- ggplot(s1, aes(U1, U2, colour = ct)) + geom_point(size=.05, alpha=.4, show.legend=FALSE) +
  scale_colour_manual(values = PAL, drop = FALSE) + facet_wrap(~ merge_sample) +
  labs(title = "by sample", subtitle = "the object is 62.6% RB1; IB1 is the composite node",
       x = "UMAP 1", y = NULL,
       caption = sprintf(paste0("axis clipped to the non-junk range; %d cells (all junk cluster) fall outside.\n",
                                "Hodgkin chemokine niche is a REGION, not a cell type, and is drawn near-black rather than in the project's Hodgkin blue,\n",
                                "which is indistinguishable from B cell at this scale. Plasma cell is halo_suspect (12.7x) -- use celltype_argmax for it."), n_out)) +
  CO + th()

sf(sprintf("celltype_umap_%s.png", RES), p1 + p2 + p3 + plot_layout(widths=c(1,1,1.4)), 21, 8)

## agreement between the two concepts, where argmax has an opinion
a <- d[lin != "no marker call" & !is.na(celltype_cluster) & celltype_cluster != ""]
a[, base := sub(",.*$", "", celltype_cluster)]
map <- c(B="B cell", T="T cell", NK="NK cell", myeloid="Myeloid", endo="Endothelial",
         fibro="Stromal/FDC", PC="Plasma cell")
cat(sprintf("\nargmax agrees with the cluster label for %.1f%% of the %s cells where argmax resolves\n",
            100*a[, mean(map[as.character(lin)] == base)], format(nrow(a), big.mark=",")))
cat("\n=== cluster label x argmax lineage ===\n")
print(dcast(a[, .N, by=.(base, lin)], base ~ lin, value.var="N", fill=0L))
cat("\ndone\n")
