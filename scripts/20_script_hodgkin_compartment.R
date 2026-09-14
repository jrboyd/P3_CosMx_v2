## ============================================================================
## 20_script_hodgkin_compartment.R -- the Hodgkin compartment, re-read with step-19 labels.
##
## WHY REPEAT IT. Step 11 tested whether IB1's CCL17/CCL22 elevation is Hodgkin biology or a
## cellular-composition artefact, and concluded the field is real. But it could only test
## WITHIN COARSE LINEAGE, and its own table shows where the signal actually sits:
##       all_markers_0   10.81x IB1/RB1        myeloid   2.49x
##       below_background 7.53x               fibro     2.50x
## The elevation is STRONGEST in the cells argmax could not call at all -- which is precisely
## the population step 19 can now name, and which step 11 had to leave as two buckets called
## "no marker call".
##
## THE QUESTION THIS STEP MUST ANSWER FIRST. CCL17, CCL22 and BATF3 are dendritic-cell genes,
## and step 19 named an 82,228-cell cluster "Hodgkin chemokine niche" on exactly those three.
## Step 11 already warned that "IB1 has more CCL17 than RB1" may just mean "IB1 has more
## T/DC-zone tissue left, because RB1's is replaced by FL". So before the label is used for
## anything: IS IT A DENDRITIC-CELL POPULATION? Discriminators, none of which named it:
##      cDC        CLEC9A, CD1C, FCER1A, LAMP3      -- if high, it is a DC population
##      pDC        CLEC4C, LILRA4, IL3RA            -- should be LOW (pDC is its own cluster)
##      HRS/cHL    TNFRSF8, IL13, CD274, PDCD1LG2, IL6, IL22RA2
##      T zone     CD3E, CCR7                       -- if high, it is interfollicular tissue
##
## AND THE SHARPEST CONTROL IS INTERNAL. The niche label is 80.5% IB1 but RB1 still holds
## 16,029 niche-labelled cells. If those RB1 cells carry the chemokines too, the label is a
## cell type present in both samples and the "Hodgkin" name is wrong. If they do not, the
## label is capturing something specific to IB1. RB1 is the no-Hodgkin control by
## construction, which is what makes this worth more than the between-sample ratio.
##
## Per-cell HRS calling stays closed (TNFRSF8 maxes at 6 counts, 78% singletons). Everything
## here is a pooled, label-level readout.
##
## Usage: sbatch --mem=200G --cpus-per-task=2 --time=3:00:00 \
##          scripts/sbatch/run_rscript.sh scripts/20_script_hodgkin_compartment.R
## ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); library(RANN)
  setDTthreads(2)
})
OUT <- "output/20_hodgkin_v2"; FIG <- file.path(OUT,"figures")
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=160,bg="#fcfcfb");cat("wrote",file.path(FIG,n),"\n")}

## TWO CIRCULARITY LEAKS CLOSED. The step-19 module that named this label is
## CCL17/CCL22/BATF3/IL22RA2. A first pass put BATF3 in the cDC set and IL22RA2 in the HRS
## set, so both of the tests that matter were partly scoring their own defining genes -- the
## cDC score came out highest for the niche (3.58) and the HRS score too (2.76), and neither
## number could be trusted. Both are removed here and kept only in `defining`, which is
## labelled as circular on the figure.
DEFINING <- c("CCL17","CCL22","BATF3","IL22RA2")
SETS <- list(
  `defining (circular)` = c("CCL17","CCL22","BATF3","IL22RA2"),
  cDC       = c("CLEC9A","CD1C","FCER1A","LAMP3"),          # BATF3 removed
  pDC       = c("CLEC4C","LILRA4","IL3RA","IRF7"),
  `HRS/cHL` = c("TNFRSF8","IL13","CD274","PDCD1LG2","IL6"), # IL22RA2 removed
  `T zone`  = c("CD3E","CCR7","CCL19","CCL21"),
  B         = c("MS4A1","CD79A"))
ORD <- c("Hodgkin chemokine niche","pDC","Myeloid","T cell","B cell","Plasmablast",
         "Plasma cell","NK cell","Endothelial","Stromal/FDC")

lab <- fread("output/19_celltypes/hc_cells_final.csv.gz", na.strings=c("NA",""))
md  <- fread("output/16_regions/cells_meta_all.csv.gz",
             select=c("cell_id","x_slide_mm","y_slide_mm","nCount_RNA"))
lab <- md[lab, on=.(cell_id)]
seu <- readRDS("output/19_celltypes/hc_02_clustered.Rds")
ca  <- if (!is.null(seu[["RNA"]]$counts)) seu[["RNA"]]$counts else seu[["RNA"]]@counts
setkey(lab, cell_id); lab <- lab[J(colnames(ca))]
G <- intersect(unlist(SETS, use.names=FALSE), rownames(ca))
miss <- setdiff(unlist(SETS, use.names=FALSE), G)
if (length(miss)) cat("not on panel: ", paste(miss, collapse=", "), "\n", sep="")
sub <- as(as.matrix(ca[G,,drop=FALSE]),"dgCMatrix"); rownames(sub) <- G
tot <- seu$nCount_RNA; rm(seu, ca); invisible(gc())
cat(sprintf("cells %d | genes %d\n", nrow(lab), length(G)))

## ---- 1. what is the niche cluster made of? ---------------------------------
prof <- rbindlist(lapply(intersect(ORD, unique(lab$celltype_cluster)), function(L)
  rbindlist(lapply(c("IB1","RB1"), function(S) {
    j <- which(lab$celltype_cluster == L & lab$merge_sample == S)
    if (length(j) < 300) return(NULL)
    data.table(label=L, samp=S, cells=length(j), gene=G,
               cp10k=as.numeric(Matrix::rowSums(sub[,j,drop=FALSE]))/sum(tot[j])*1e4,
               det=as.numeric(Matrix::rowSums(sub[,j,drop=FALSE]>0))/length(j))
  }))))
prof[, set := {m <- rep(NA_character_, .N)
  for (nm in names(SETS)) m[gene %in% SETS[[nm]]] <- nm; m}]
sc("hodgkin_label_profile.csv", prof)

mod <- prof[!is.na(set), .(cp10k=mean(cp10k), det=mean(det)), by=.(label, samp, set)]
cat("\n=== 1. IS THE NICHE A DENDRITIC-CELL POPULATION? module cp10k, IB1 ===\n")
print(dcast(mod[samp=="IB1"], label ~ set, value.var="cp10k")[
        , lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)][
        order(match(label, ORD))])

## ---- 2. the internal control: the niche label in RB1 -----------------------
cat("\n=== 2. THE INTERNAL CONTROL: the same labels in RB1 (no Hodgkin component) ===\n")
print(dcast(mod[samp=="RB1"], label ~ set, value.var="cp10k")[
        , lapply(.SD, function(x) if (is.numeric(x)) round(x,2) else x)][
        order(match(label, ORD))])
cat("\n=== IB1 / RB1 ratio per label, per module ===\n")
rr <- dcast(mod, label + set ~ samp, value.var="cp10k")
rr[, ratio := round(IB1/RB1, 2)]
print(dcast(rr, label ~ set, value.var="ratio")[order(match(label, ORD))])
sc("hodgkin_label_ratios.csv", rr)
cat("\nStep 11's coarse-lineage IB1/RB1 chemokine ratios ran 2.49x (myeloid) to 10.81x\n")
cat("(all_markers_0). A label whose ratio is near 1 is the same population in both samples.\n")

## ---- 3. spatial: does the niche label sit on step 11's called bins? --------
nb <- fread("output/11_hodgkin_niche/niche_bins_called.csv")[samp=="IB1" & niche==TRUE]
cat(sprintf("\n=== 3. step-11 called niche bins: %d ===\n", nrow(nb)))
ib <- lab[merge_sample=="IB1" & !is.na(x_slide_mm)]
r <- RANN::nn2(nb[, .(x_mm, y_mm)], ib[, .(x_slide_mm, y_slide_mm)], k=1L)
ib[, d_hot_mm := r$nn.dists[,1]]
ib[, dbin := cut(d_hot_mm, c(-Inf,.1,.25,.5,1,2,Inf),
                 labels=c("<0.1","0.1-0.25","0.25-0.5","0.5-1","1-2",">2"))]
en <- ib[, .(cells=.N, pct_niche=round(100*mean(celltype_cluster=="Hodgkin chemokine niche",
                                                na.rm=TRUE),2)), by=dbin][order(dbin)]
en[, fold_vs_far := round(pct_niche/pct_niche[.N], 2)]
cat("share of IB1 cells carrying the niche label, by distance to a step-11 hotspot:\n")
print(en)
sc("hodgkin_label_vs_step11_bins.csv", en)

SURF<-"#fcfcfb"
hm <- prof[!is.na(set) & samp=="IB1"]
hm[, z := (cp10k-mean(cp10k))/pmax(sd(cp10k),1e-9), by=gene]
hm[, label := factor(label, levels=ORD)]
hm[, set := factor(set, levels=names(SETS))]
p <- ggplot(hm, aes(label, gene, fill=z)) + geom_tile(colour=SURF, linewidth=.3) +
  scale_fill_gradient2(low="#184f95", mid="#f0efec", high="#d03b3b", midpoint=0,
                       limits=c(-1.5,1.5), oob=scales::squish, name="z across labels (IB1)") +
  facet_grid(set ~ ., scales="free_y", space="free_y", switch="y") +
  labs(title="Is the Hodgkin chemokine niche a dendritic-cell population?",
       subtitle="IB1 only. Top block is the defining module and is circular. Everything below it is independent of it:\nBATF3 is excluded from cDC and IL22RA2 from HRS, because both are defining genes and would otherwise score themselves.",
       x=NULL, y=NULL) +
  theme_minimal(10) + theme(plot.background=element_rect(fill=SURF,colour=NA),
    panel.grid=element_blank(), axis.text.x=element_text(angle=35,hjust=1),
    strip.text.y.left=element_text(angle=0,face="bold",size=7), strip.placement="outside",
    plot.title=element_text(face="bold"), plot.subtitle=element_text(size=8),
    legend.position="bottom")
sf("hodgkin_identity.png", p, 9, 8)
cat("\ndone\n")
