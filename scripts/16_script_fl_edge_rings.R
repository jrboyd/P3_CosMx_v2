## ============================================================================
## 16_script_fl_edge_rings.R -- expression as a function of distance from the FL edge.
##
## WHY RINGS RATHER THAN IN/OUT. At n = 1 patient a two-group p-value is weak evidence: it
## says two arbitrary groups differ. A MONOTONE gradient across five distance bands is much
## harder to produce by accident, and it is the shape a real boundary effect should have. This
## replaces step 14's binary in/out with a dose-response, in both samples.
##
## THE EDGE IS DEFINED SYMMETRICALLY, which is the only way the two nodes are comparable.
## 16c established that IB1's FL foci and RB1's follicles do not overlap in composition, so an
## IB1 light-chain edge and an RB1 architecture edge would not be the same object. Instead the
## boundary in BOTH samples is the B_dense / not-B_dense bin boundary from the shared `arch`
## axis, and the signed distance is:
##     negative  = inside B_dense tissue, distance to the nearest non-B_dense cell
##     positive  = outside it, distance to the nearest B_dense cell
## IB1's light-chain FL edge is run alongside as a sensitivity, because it is the definition
## the rest of the project uses -- but it exists for IB1 only.
##
## THE UNIT IS THE FOV, NOT THE CELL. Cells in a ring are not independent, and rings are
## nested inside FOVs. So expression is pseudobulked per (FOV, ring), a per-FOV Spearman is
## taken across rings, and the test is a signed-rank of those per-FOV rho values against zero.
## An FOV therefore contributes one number per gene, which is what makes the p-value mean
## something. Genes are BH-corrected within sample.
##
## Usage: Rscript scripts/16_script_fl_edge_rings.R [tier]
## ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(RANN); library(data.table); library(ggplot2); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "matched"

REG_F <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
PB    <- file.path("output/16_pseudobulk", TIER)     # for the gene list only
CNT   <- "output/16_regions/panel_counts_all.Rds"
OUT   <- file.path("output/16_rings", TIER)
FIG   <- file.path(OUT, "figures")
for (f in c(REG_F, CNT)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

BREAKS <- c(-Inf, -40, 0, 40, 100, 200, 500, Inf)
LABS   <- c("in >40","in 0-40","out 0-40","out 40-100","out 100-200","out 200-500","out >500")
MIN_CELLS_FR <- 20L   # cells for a (FOV, ring) pseudobulk to be usable
MIN_RINGS    <- 4L    # rings an FOV must cover to contribute a rho
MIN_FOV      <- 12L   # FOVs contributing before a gene is testable
SEED <- 1L
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
SAMP_COLS <- c(IB1="#2a78d6", RB1="#eb6834")

cat("============================================================\n")
cat("step 16r -- FL-edge distance rings, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

d <- fread(REG_F, select=c("cell_id","samp","bin_region","arch","x_slide_mm","y_slide_mm",
                           "is_defect","unbinned","fov_id","nCount_RNA","d_fl_um"))
d <- d[is_defect==FALSE & unbinned==FALSE & !is.na(arch)]
M <- readRDS(CNT); jj <- match(d$cell_id, colnames(M)); stopifnot(!anyNA(jj))
M <- M[, jj, drop=FALSE]; GENES <- rownames(M)
cat(sprintf("cells %d | genes %d (step-16 panel)\n", nrow(d), length(GENES)))

## ---- signed distance to the B_dense boundary --------------------------------
d[, inside := arch == "B_dense"]
for (S in unique(d$samp)) {
  ii <- d[, which(samp == S)]
  z  <- d[ii]
  xin  <- as.matrix(z[inside == TRUE,  .(x_slide_mm*1000, y_slide_mm*1000)])
  xout <- as.matrix(z[inside == FALSE, .(x_slide_mm*1000, y_slide_mm*1000)])
  if (!nrow(xin) || !nrow(xout)) next
  xy <- as.matrix(z[, .(x_slide_mm*1000, y_slide_mm*1000)])
  ## inside cells: how far to the nearest OUTSIDE cell (negative); outside: to nearest INSIDE
  dn <- numeric(nrow(z))
  dn[z$inside]  <- -RANN::nn2(xout, xy[z$inside, , drop=FALSE],  k=1)$nn.dists[,1]
  dn[!z$inside] <-  RANN::nn2(xin,  xy[!z$inside, , drop=FALSE], k=1)$nn.dists[,1]
  d[ii, d_edge_um := dn]
}
d[, ring := cut(d_edge_um, BREAKS, labels=LABS, right=FALSE)]
cat("\n=== cells per ring (signed distance to the B_dense boundary)\n")
print(dcast(d[!is.na(ring), .N, by=.(samp, ring)], ring ~ samp, value.var="N", fill=0L))

## ---- per (FOV, ring) pseudobulk ---------------------------------------------
ring_profile <- function(dist_col, tag, samps) {
  z <- d[samp %in% samps & !is.na(get(dist_col))]
  z[, rg := cut(get(dist_col), BREAKS, labels=LABS, right=FALSE)]
  z <- z[!is.na(rg)]
  z[, grp := paste(samp, fov_id, rg, sep="|")]
  g <- z[, .(cells=.N, samp=samp[1], fov=fov_id[1], rg=rg[1], nCount=sum(nCount_RNA)), by=grp]
  g <- g[cells >= MIN_CELLS_FR]
  if (!nrow(g)) return(NULL)
  z <- z[grp %in% g$grp]
  gi <- match(z$grp, g$grp)
  IND <- sparseMatrix(i=seq_len(nrow(z)), j=gi, x=1, dims=c(nrow(z), nrow(g)))
  SUM <- as.matrix(M[, match(z$cell_id, colnames(M)), drop=FALSE] %*% IND)
  CP  <- sweep(SUM, 2, pmax(g$nCount,1), "/") * 1e4
  dimnames(CP) <- list(GENES, g$grp)
  cat(sprintf("\n[%s] (FOV x ring) groups with >= %d cells: %d across %d FOVs\n",
              tag, MIN_CELLS_FR, nrow(g), uniqueN(g$fov)))
  print(g[, .(groups=.N, fovs=uniqueN(fov), med_cells=as.numeric(median(cells))), by=.(samp,rg)][order(samp,rg)])
  ## per-FOV Spearman of cp10k against ring index, then signed-rank across FOVs
  g[, ridx := as.integer(rg)]
  res <- rbindlist(lapply(unique(g$samp), function(S) {
    gg <- g[samp == S]
    fv <- gg[, .N, by=fov][N >= MIN_RINGS, fov]
    if (length(fv) < MIN_FOV) {
      cat(sprintf("  %s: only %d FOVs cover >= %d rings; not testable\n", S, length(fv), MIN_RINGS))
      return(NULL) }
    rho <- vapply(fv, function(f) {
      k <- gg[fov == f, which = TRUE]; cols <- gg$grp[k]; ri <- gg$ridx[k]
      suppressWarnings(apply(CP[, cols, drop=FALSE], 1, function(v) cor(v, ri, method="spearman")))
    }, numeric(length(GENES)))
    med <- apply(rho, 1, median, na.rm=TRUE)
    p   <- apply(rho, 1, function(v) { v <- v[is.finite(v)]
      if (length(v) < MIN_FOV) NA_real_ else
        suppressWarnings(wilcox.test(v, mu=0)$p.value) })
    data.table(set=tag, samp=S, gene=GENES, n_fov=length(fv),
               med_rho=med, p=p)[, q := p.adjust(p, "BH")]
  }), fill=TRUE)
  list(res=res, cp=CP, meta=g)
}

## primary: the symmetric B_dense boundary, both samples
A <- ring_profile("d_edge_um", "B_dense boundary", c("IB1","RB1"))
## sensitivity: IB1's light-chain FL edge, the definition the rest of the project uses
B <- ring_profile("d_fl_um", "IB1 light-chain FL edge", "IB1")
RES <- rbindlist(list(A$res, B$res), fill=TRUE)
if (!nrow(RES)) stop("no set was testable", call.=FALSE)
sc("ring_trends.csv", RES)

for (k in unique(paste(RES$set, RES$samp))) {
  z <- RES[paste(set,samp)==k & q < 0.05]
  cat(sprintf("\n=== %s: %d of %d genes with a monotone trend across rings (q<0.05)\n",
              k, nrow(z), RES[paste(set,samp)==k, .N]))
  if (nrow(z)) {
    cat("  rising with distance FROM the FL edge (higher outside):\n")
    print(head(z[med_rho>0][order(-med_rho), .(gene, med_rho=round(med_rho,3), n_fov, q=signif(q,2))], 8))
    cat("  falling with distance (higher INSIDE the FL):\n")
    print(head(z[med_rho<0][order(med_rho), .(gene, med_rho=round(med_rho,3), n_fov, q=signif(q,2))], 8))
  }
}

## Do the two samples agree on which genes are edge-graded?
w <- dcast(RES[set=="B_dense boundary"], gene ~ samp, value.var=c("med_rho","q"))
if (all(c("med_rho_IB1","med_rho_RB1") %in% names(w))) {
  w[, both := q_IB1 < 0.05 & q_RB1 < 0.05 & sign(med_rho_IB1)==sign(med_rho_RB1)]
  cat(sprintf("\n=== genes edge-graded in BOTH nodes, same direction: %d\n", w[both==TRUE, .N]))
  print(head(w[both==TRUE][order(-abs(med_rho_IB1)),
        .(gene, IB1=round(med_rho_IB1,3), RB1=round(med_rho_RB1,3))], 15))
  cat(sprintf("\ncorrelation of med_rho between samples across all genes: %.3f\n",
              w[, cor(med_rho_IB1, med_rho_RB1, use="complete.obs")]))
  sc("ring_trend_concordance.csv", w)
}

## ---- figure: the checkpoint/Tfh panel across rings ---------------------------
SHOW <- intersect(c("PDCD1","CTLA4","LAG3","CXCL13","CXCR5","ICOS","FOXP3","CD274",
                    "CCL17","CCL22","BCL2","MS4A1"), GENES)
if (!is.null(A) && length(SHOW)) {
  cp <- A$cp; mt <- A$meta
  pl <- rbindlist(lapply(SHOW, function(g)
    data.table(gene=g, samp=mt$samp, rg=mt$rg, cells=mt$cells, v=cp[g, mt$grp])))
  pm <- pl[, .(mean=weighted.mean(v, cells), se=sd(v)/sqrt(.N), n=.N), by=.(gene,samp,rg)]
  pm[, gene := factor(gene, levels=SHOW)]
  p <- ggplot(pm, aes(rg, mean, colour=samp, group=samp)) +
    geom_vline(xintercept=2.5, linetype="22", colour=MUTED, linewidth=.4) +
    geom_line(linewidth=.65) + geom_point(size=1.7) +
    geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se), width=.15, linewidth=.35) +
    facet_wrap(~ gene, scales="free_y", ncol=4) +
    scale_colour_manual(values=SAMP_COLS, name=NULL) +
    labs(x="signed distance from the B_dense (FL) boundary", y="cp10k per (FOV, ring)",
         title="Expression across the FL boundary, both nodes",
         subtitle=paste0("Dashed line is the boundary; left of it is inside FL tissue. Rings ",
           "defined identically in both samples from the shared\narchitecture axis, because ",
           "the two nodes' FL territory does not overlap in composition (16c). FOV is the unit."),
         caption="Data: output/16_rings/<tier>/ring_trends.csv | 95% CI") +
    theme_minimal(base_size=10) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_line(colour=GRID, linewidth=.25),
          axis.text.x=element_text(colour=MUTED, angle=45, hjust=1, size=6.5),
          axis.text.y=element_text(colour=MUTED, size=7),
          axis.title=element_text(colour=INK2, size=9),
          strip.text=element_text(colour=INK2, face="bold", size=8),
          plot.title=element_text(colour=INK, face="bold", size=12),
          plot.subtitle=element_text(colour=INK2, size=8),
          plot.caption=element_text(colour=MUTED, size=7, hjust=0),
          legend.position="top", legend.justification="left")
  ggsave(file.path(FIG,"fig_16_fl_edge_rings.png"), p, width=12, height=8, dpi=170, bg=SURF)
  cat("wrote ", file.path(FIG,"fig_16_fl_edge_rings.png"), "\n", sep="")
}
writeLines(format(Sys.time()), file.path(OUT,"rings.done"))
cat("============================================================\nOK\n")
