## ============================================================================
## 19_script_composition_figure.R -- the step-19 composition figure, three panels.
##
##   A  the two sections on the same embedding, coloured by cell type
##   B  relative abundance of each type WITHIN each sample (each sample sums to 100%)
##   C  for each (type, sample), what the per-cell argmax evidence says about those cells
##
## Panel C is the one that earns the figure. A cluster label is guilt by association: every
## cell in a cluster takes its name, including the ~43% with no marker call of their own.
## Panel C opens each bar up and shows how much independent marker support is actually
## underneath -- the matching-colour segment is agreement, everything else is a cell the
## label and the markers disagree about, and grey is a cell the markers had no opinion on.
## Overall that agreement runs at 66.2% across the high-confidence set.
##
## Percentages in B are WITHIN sample, not of the object. The object is 63.6% RB1, so
## object-level shares are dominated by RB1 throughout and would hide IB1 entirely.
##
## Usage: sbatch --mem=48G --cpus-per-task=2 --time=1:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/19_script_composition_figure.R
## ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); setDTthreads(2)
})
OUT <- "output/19_celltypes"; FIG <- file.path(OUT,"figures")
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=170,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
SAMP_COLS <- c(IB1="#2a78d6", RB1="#eb6834")   # project-wide, CLAUDE.md
PAL <- c("B cell"="#184f95", "Plasmablast"="#f2789f", "Plasma cell"="#c51b7d", "T cell"="#d03b3b",
         "NK cell"="#7b3f9e", "Myeloid"="#e08214", "pDC"="#00857a",
         "Endothelial"="#1b7837", "Stromal/FDC"="#8c6d31",
         "Hodgkin chemokine niche"="#151a2b")
LPAL <- c(B="#184f95", PC="#c51b7d", T="#d03b3b", NK="#7b3f9e", myeloid="#e08214",
          endo="#1b7837", fibro="#8c6d31", `no marker call`="#dcdcd6")
ORD <- c("B cell","Plasmablast","Plasma cell","T cell","NK cell","Myeloid","pDC",
         "Endothelial","Stromal/FDC","Hodgkin chemokine niche")
th <- function(base=11) theme_minimal(base_size=base) +
  theme(plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour=GRID,linewidth=.25),
        axis.text=element_text(colour=MUTED), axis.title=element_text(colour=INK2,size=base-1),
        plot.title=element_text(colour=INK,face="bold",size=base+1),
        plot.subtitle=element_text(colour=INK2,size=base-2),
        plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
        strip.text=element_text(colour=INK2,face="bold"),
        legend.text=element_text(colour=INK2,size=base-2), legend.title=element_blank())

F_IN <- file.path(OUT,"hc_cells_final.csv.gz")
if (!file.exists(F_IN)) F_IN <- file.path(OUT,"hc_cells_labelled_high.csv.gz")
cat("reading ", F_IN, "\n", sep=""); d <- fread(F_IN, na.strings=c("NA",""))
d <- d[!is.na(celltype_cluster)]
d[, ct := factor(celltype_cluster, levels=ORD)]
d[, lin := factor(fifelse(lineage_argmax %in% names(LPAL), lineage_argmax, "no marker call"),
                  levels=names(LPAL))]
cat(sprintf("cells %d | IB1 %d | RB1 %d\n", nrow(d),
            d[merge_sample=="IB1",.N], d[merge_sample=="RB1",.N]))

## ---- A: the two sections ----------------------------------------------------
## Axis range from the non-junk cells and coord_fixed: a UMAP is a distance embedding, and
## 159 near-empty cells previously stretched x to 4.5x the width the data needed.
XL <- range(d$U1); YL <- range(d$U2)
pad <- .04*c(diff(XL), diff(YL)); XL <- XL+c(-pad[1],pad[1]); YL <- YL+c(-pad[2],pad[2])
set.seed(1L); s <- d[sample(.N, min(.N, 400000L))]
pA <- ggplot(s, aes(U1,U2,colour=ct)) + geom_point(size=.05, alpha=.4) +
  scale_colour_manual(values=PAL, drop=FALSE) +
  guides(colour=guide_legend(override.aes=list(size=2.5,alpha=1), ncol=3)) +
  facet_wrap(~merge_sample) +
  coord_fixed(ratio=1, xlim=XL, ylim=YL, expand=FALSE) +
  labs(title="A. the two sections on one embedding", x="UMAP 1", y="UMAP 2",
       subtitle=sprintf("%s cells shown of %s high-confidence", format(nrow(s),big.mark=","), format(nrow(d),big.mark=","))) +
  th() + theme(legend.position="bottom", panel.grid.major.x=element_line(colour=GRID,linewidth=.2))

## ---- B: abundance within sample ---------------------------------------------
ab <- d[, .(cells=.N), by=.(ct, merge_sample)][, pct := 100*cells/sum(cells), by=merge_sample]
sc("composition_by_sample.csv", ab[order(ct, merge_sample)])
pB <- ggplot(ab, aes(ct, pct, fill=merge_sample)) +
  geom_col(position=position_dodge(width=.72), width=.68) +
  geom_text(aes(label=sprintf("%.1f", pct)), position=position_dodge(width=.72),
            vjust=-0.35, size=2.5, colour=INK2) +
  scale_fill_manual(values=SAMP_COLS) +
  scale_y_continuous(expand=expansion(mult=c(0,.12))) +
  labs(title="B. relative abundance within each sample",
       subtitle="Each sample sums to 100%. The object is 63.6% RB1, so object-level shares would hide IB1.",
       x=NULL, y="% of that sample's high-confidence cells") +
  th() + theme(legend.position=c(.88,.85),
               axis.text.x=element_text(angle=35,hjust=1,colour=INK2))

## ---- C: what the markers say about each label -------------------------------
cc <- d[, .(cells=.N), by=.(ct, merge_sample, lin)][, frac := cells/sum(cells), by=.(ct,merge_sample)]
sc("composition_argmax_within_label.csv", cc[order(ct, merge_sample, lin)])
ag <- d[, .(agree_pct = 100*mean(as.character(lin) ==
        c("B cell"="B","Plasma cell"="PC","T cell"="T","NK cell"="NK","Myeloid"="myeloid",
          "Endothelial"="endo","Stromal/FDC"="fibro","Plasmablast"="PC","pDC"="__none__",
          "Hodgkin chemokine niche"="__none__")[as.character(ct)])), by=.(ct, merge_sample)]
pC <- ggplot(cc, aes(merge_sample, frac, fill=lin)) +
  geom_col(width=.75) +
  geom_text(data=ag[agree_pct>0], aes(x=merge_sample, y=1.04, label=sprintf("%.0f%%", agree_pct)),
            inherit.aes=FALSE, size=2.4, colour=INK2) +
  facet_grid(~ ct, switch="x") +
  scale_fill_manual(values=LPAL) +
  scale_y_continuous(labels=scales::percent, expand=expansion(mult=c(0,.1))) +
  guides(fill=guide_legend(nrow=1)) +
  labs(title="C. per-cell marker evidence underneath each cluster label",
       subtitle="Each bar is one (cell type, sample). Colour = what that cell's OWN markers say; grey = no marker call.\nThe matching colour is agreement; the number above each bar is that agreement %. pDC and the Hodgkin niche have no argmax family, so none is marked.",
       x=NULL, y="fraction of the label's cells",
       caption="pDC and Hodgkin chemokine niche were named by cluster-level module rules because LINEAGE_SETS_COARSE has no family for them --\nwhich is why their bars carry no matching colour. Plasmablast is scored against the PC family, argmax having no plasmablast set.\nThe earlier halo_suspect flag on Plasma cell is withdrawn: it measured spatial aggregation, not leakage (pc_halo_test.csv).") +
  th() + theme(legend.position="top", axis.text.x=element_text(size=7,colour=INK2),
               strip.text=element_text(size=7,face="bold"), panel.spacing.x=unit(2,"pt"))

sf("composition_figure.png", (pA | pB) / pC + plot_layout(heights=c(1,1.05), widths=c(1.25,1)), 17, 12)
cat("\n=== abundance within sample (%) ===\n")
print(dcast(ab, ct ~ merge_sample, value.var="pct")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,2) else x)])
cat("\n=== agreement between label and argmax, per sample ===\n")
print(dcast(ag[agree_pct>0], ct ~ merge_sample, value.var="agree_pct")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,1) else x)])
cat("\n=== share of each label with NO marker call ===\n")
print(dcast(cc[lin=="no marker call"], ct ~ merge_sample, value.var="frac")[, lapply(.SD, function(x) if(is.numeric(x)) round(100*x,1) else x)])
cat("\ndone\n")
