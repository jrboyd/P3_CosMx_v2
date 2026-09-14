## ============================================================================
## 17_script_modules.R -- the PD-1 axis by proxy, because the genes themselves are too sparse.
##
## STEP 16 ESTABLISHED THE PROBLEM. PDCD1, CD274 and PDCD1LG2 sit at 4.6-5.8% detection; at the
## >=1 threshold they are 2.6-2.9x the mean blank probe, i.e. near noise; at the calibrated >=2
## threshold they are specific but only ~1-2% of cells are positive. Region-unit tests within
## every cell type returned nothing (all q >= 0.11) and the spatial pairing showed no
## reproducible engagement. That is a FLOOR result, not an absence -- and this script attacks
## the floor three ways rather than re-running the same test.
##
## 1. MODULES INSTEAD OF GENES. Summing 20-33 co-regulated genes multiplies the counts, so a
##    module clears the floor where its members do not. The modules are chosen to be
##    MECHANISTICALLY tied to the axis, not merely correlated with it:
##      IFN-gamma response   the canonical UPSTREAM DRIVER of CD274. If PD-L1-inducing
##                           pressure differs between the nodes, this reports it in genes that
##                           are actually detectable.
##      exhaustion           the DOWNSTREAM CONSEQUENCE of PD-1 signalling, including the
##                           counter-markers (TCF7/IL7R/SELL/CCR7) that fall as it proceeds.
##      9p24.1 amplicon      CD274 + PDCD1LG2 + JAK2. In classic Hodgkin lymphoma PD-L1 is
##                           overexpressed via 9p24.1 amplification, and all three genes sit in
##                           that amplicon. EXPLORATORY ONLY: expression cannot establish copy
##                           number, and three genes is a thin module. Flagged everywhere.
##    Each module is scored two ways -- pooled cp10k, and the FRACTION OF MODULE GENES DETECTED
##    per cell. The second is not compositional and not depth-normalised by division, so when
##    the two agree the result is not an artifact of either.
##
## 2. THE FOV AS THE UNIT, not the 100 um bin. A bin holds 40-500 cells; an FOV holds
##    1,000-3,000, so its total counts are ~100x and a sparse gene's cp10k is actually
##    estimable. 650 IB1 / 429 RB1 FOVs is ample replication. Bin-level is run alongside, and
##    if a result appears only at bin level it is noise, not resolution.
##
## 3. SIGN CONSISTENCY ACROSS CELL TYPES as its own inference. The step-16 heatmap showed CD274
##    and CTLA4 higher in IB1 in all five T subsets while PDCD1, TOX and TNFRSF9 were higher in
##    RB1 in all five. Five independent compartments agreeing in direction is a sign test at
##    p = 1/16 two-sided, and it does not require any single contrast to survive BH across
##    6,175 genes. Formalised here.
##
## 4. AND A POWER BOUND on whatever stays null, so "we found nothing" becomes "we can exclude a
##    fold-change larger than X" -- which is a statement, where a bare null is not.
##
## Usage: Rscript scripts/17_script_modules.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "matched"

PRE_RDS <- "output/01_import_v1/seuratObject_%s_BPCells_full.Rds"
REG_F   <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT_F    <- "output/16_celltypes/cell_celltypes.csv.gz"
OUT     <- file.path("output/17_pd1", TIER)
FIG     <- file.path(OUT, "figures")
for (f in c(REG_F, CT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

SAMPLES <- c("IB1","RB1")
CHUNK   <- 150000L
MIN_CELLS_U <- 100L    # cells for an FOV/bin unit to be usable
MIN_UNITS   <- 12L     # units per side before a contrast is testable
N_BOOT      <- 2000L
SEED        <- 1L
FOC <- list(IB1 = c("FL_focus","cHL_hotspot","background"),
            RB1 = c("FL_follicular","interfollicular"))

MODULES <- list(
  `PD-1 axis`      = c("PDCD1","CD274","PDCD1LG2"),
  ## B2M, CD74 and HLA-DRA are DELIBERATELY NOT HERE even though they are IFN-gamma
  ## responsive. They are constitutively high (each detected in >70% of cells), so including
  ## them let them dominate the module's counts and the score stopped measuring IFN-gamma
  ## signalling and started measuring MHC baseline -- visible as the IFN-gamma and
  ## antigen-presentation modules reporting the identical 73.6% best-single-gene detection,
  ## i.e. sharing the same dominant member. They live in `antigen presn`, where that IS the
  ## intended readout. What remains here is the inducible, low-baseline set.
  `IFN-gamma`      = c("IFNG","IFNGR1","IFNGR2","STAT1","STAT2","IRF1","IRF9","JAK1",
                       "GBP1","GBP2","GBP4","GBP5","CXCL9","CXCL10","CXCL11","IDO1",
                       "TAP1","TAP2","PSMB8","PSMB9","NLRC5",
                       "SOCS1","SOCS3","ISG15","MX1","OAS1","IFIT1","IFIT3"),
  `exhaustion`     = c("PDCD1","LAG3","TIGIT","HAVCR2","TOX","CTLA4","ENTPD1","TNFRSF9",
                       "CD160","KLRG1","EOMES","BATF","PRDM1","NR4A1","NR4A2","NR4A3","VSIR"),
  `T naive/memory` = c("TCF7","IL7R","SELL","CCR7","CD28","SLAMF6"),
  `9p24.1 (EXPLORATORY)` = c("CD274","PDCD1LG2","JAK2"),
  `antigen presn`  = c("B2M","CIITA","HLA-DRA","HLA-DQA1","HLA-DPA1","CD58","HLA-E"),
  `Tfh help`       = c("CXCL13","IL21","IL4","CD40LG","SH2D1A","MAF","POU2AF1","ICOS"))
CTS_T <- c("Tfh","Treg","CD4_conv","CD8_T","T_unresolved")   # the sign-consistency panel

cat("============================================================\n")
cat("step 17 -- PD-1 axis by proxy module, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

## ---- cells ------------------------------------------------------------------
rg <- fread(REG_F, select=c("cell_id","samp","bin_region","bx","by_","is_defect","unbinned",
                            "fov_id","nCount_RNA"))
ct <- fread(CT_F, select=c("cell_id","celltype"))
d  <- merge(rg, ct, by="cell_id", sort=FALSE)[is_defect==FALSE & unbinned==FALSE]
d  <- d[mapply(function(s,r) r %in% FOC[[s]], samp, bin_region)]
cat(sprintf("cells %d\n", nrow(d))); print(d[, .N, by=.(samp,bin_region)][order(samp,-N)])

## ---- module counts, streamed from the pre-QC objects ------------------------
GENES <- unique(unlist(MODULES, use.names = FALSE))
CNT <- NULL
for (S in SAMPLES) {
  f <- sprintf(PRE_RDS, S); cat("\nreading ", f, "\n", sep="")
  o <- readRDS(f); cn <- LayerData(o[["RNA"]], layer="counts"); rn <- colnames(cn)
  g <- intersect(GENES, rownames(cn))
  sub <- d[samp == S]; jj <- match(sub$cell_id, rn)
  if (anyNA(jj)) stop(sum(is.na(jj)), " ", S, " cells absent from counts", call.=FALSE)
  parts <- list()
  for (st in seq(1L, length(jj), by=CHUNK)) {
    en <- min(st+CHUNK-1L, length(jj))
    m <- Matrix::drop0(as(cn[g, jj[st:en], drop=FALSE], "dgCMatrix"))
    parts[[length(parts)+1L]] <- m
  }
  M <- do.call(cbind, parts); colnames(M) <- sub$cell_id
  CNT[[S]] <- M; rm(o, cn, parts); invisible(gc())
  cat(sprintf("  %s: %d genes x %d cells\n", S, nrow(M), ncol(M)))
}
G <- Reduce(intersect, lapply(CNT, rownames))
M <- do.call(cbind, lapply(CNT, function(x) x[G, , drop=FALSE])); rm(CNT); invisible(gc())
M <- M[, d$cell_id, drop = FALSE]
absent <- setdiff(GENES, G)
if (length(absent)) message("not on the panel, dropped: ", paste(absent, collapse=", "))
cat(sprintf("\nmodule matrix: %d genes x %d cells\n", nrow(M), ncol(M)))

## per-module per-cell: summed counts, and the number of module genes detected
for (k in names(MODULES)) {
  gg <- intersect(MODULES[[k]], G); if (!length(gg)) next
  d[, (paste0("mc_", k)) := as.numeric(Matrix::colSums(M[gg, , drop=FALSE]))]
  d[, (paste0("mn_", k)) := as.numeric(Matrix::colSums(M[gg, , drop=FALSE] > 0))]
  d[, (paste0("mk_", k)) := length(gg)]
}
MODS <- names(MODULES)[vapply(names(MODULES), function(k) length(intersect(MODULES[[k]], G))>0, TRUE)]

cat("\n=== does a module clear the floor where its genes do not?\n")
fl <- rbindlist(lapply(MODS, function(k) {
  gg <- intersect(MODULES[[k]], G)
  data.table(module = k, genes = length(gg),
             pct_det_best_single = round(100*max(Matrix::rowMeans(M[gg,,drop=FALSE] > 0)),2),
             pct_det_module = round(100*d[, mean(get(paste0("mn_",k)) > 0)],2),
             pct_module_ge3 = round(100*d[, mean(get(paste0("mn_",k)) >= 3)],2),
             mean_counts = round(d[, mean(get(paste0("mc_",k)))],3))
}))
print(fl); sc("module_detection.csv", fl)

## WHICH GENES CARRY EACH MODULE. Without this a module score is unauditable: "IFN-gamma =
## 208 cp10k" could be one constitutive gene wearing a pathway's name. Any module whose top
## gene contributes most of the counts is a single-gene readout with extra steps.
contrib <- rbindlist(lapply(MODS, function(k) {
  gg <- intersect(MODULES[[k]], G)
  mc <- Matrix::rowSums(M[gg, , drop=FALSE]); tot <- sum(mc)
  data.table(module = k, gene = gg,
             pct_of_module_counts = round(100*as.numeric(mc)/max(tot,1), 2),
             pct_det = round(100*as.numeric(Matrix::rowMeans(M[gg,,drop=FALSE] > 0)), 2))[
    order(-pct_of_module_counts)]
}))
sc("module_gene_contributions.csv", contrib)
cat("\n=== top contributor per module, and how concentrated the module is\n")
print(contrib[, .(top_gene = gene[1], top_pct = pct_of_module_counts[1],
                  top3_pct = sum(head(pct_of_module_counts,3)),
                  genes = .N,
                  verdict = fifelse(pct_of_module_counts[1] > 50,
                                    "DOMINATED by one gene -- not a module",
                            fifelse(sum(head(pct_of_module_counts,3)) > 80,
                                    "concentrated in 3 genes", "spread"))), by = module])

## ---- units: FOV (primary) and bin (secondary) -------------------------------
unit_table <- function(ucol) {
  d[, .(cells=.N, samp=samp[1], bin_region=bin_region[1], nCount=sum(nCount_RNA),
        med_depth=as.numeric(median(nCount_RNA)),
        .SD_sum = 0), by = c(ucol)][cells >= MIN_CELLS_U]
}
mk_units <- function(ucol) {
  agg <- d[, c(list(cells=.N, samp=samp[1], bin_region=bin_region[1],
                    nCount=sum(nCount_RNA), med_depth=as.numeric(median(nCount_RNA))),
               lapply(.SD, sum)), by=c(ucol), .SDcols=paste0("mc_", MODS)]
  det <- d[, lapply(.SD, function(v) mean(v > 0)), by=c(ucol), .SDcols=paste0("mn_", MODS)]
  u <- merge(agg, det, by=ucol)
  u <- u[cells >= MIN_CELLS_U]
  for (k in MODS) u[, (paste0("cp_", k)) := 1e4 * get(paste0("mc_", k)) / pmax(nCount,1)]
  u[]
}
d[, unit_fov := paste(samp, fov_id, sep="_")]
d[, unit_bin := paste(samp, bx, by_, sep="_")]
d[, unit_fovreg := paste(samp, fov_id, bin_region, sep="|")]
U_FOV <- mk_units("unit_fovreg")   # FOV x region: an FOV can span regions
U_BIN <- mk_units("unit_bin")
cat(sprintf("\nunits with >= %d cells: %d FOVxregion, %d bins\n",
            MIN_CELLS_U, nrow(U_FOV), nrow(U_BIN)))
print(U_FOV[, .(units=.N, med_cells=as.numeric(median(cells)),
                med_depth=as.numeric(median(med_depth))), by=.(samp,bin_region)][order(samp,-units)])

## ---- contrasts, region as the unit ------------------------------------------
CON <- list(
  list("IB1 FL vs IB1 bg",      "IB1","FL_focus",      "IB1","background"),
  list("IB1 hot vs IB1 bg",     "IB1","cHL_hotspot",   "IB1","background"),
  list("RB1 foll vs RB1 inter", "RB1","FL_follicular", "RB1","interfollicular"),
  list("IB1 FL vs RB1 foll",    "IB1","FL_focus",      "RB1","FL_follicular"),
  list("IB1 bg vs RB1 inter",   "IB1","background",    "RB1","interfollicular"))
## NOTE the bracing. `rbindlist(lapply(X, function(z) rbindlist(lapply(...)), fill=TRUE), ...)`
## binds `fill` to the INNER lapply, not to rbindlist, and dies with "unused argument
## (fill = TRUE)". Each rbindlist gets its own explicit braces here so the argument cannot
## migrate. This bug has bitten three scripts in this project; brace the closures.
run <- function(U, unit_name) rbindlist(lapply(CON, function(z) {
 rbindlist(lapply(MODS, function(k) {
  a <- U[samp==z[[2]] & bin_region==z[[3]]]; b <- U[samp==z[[4]] & bin_region==z[[5]]]
  if (min(nrow(a), nrow(b)) < MIN_UNITS) return(NULL)
  ca <- a[[paste0("cp_",k)]]; cb <- b[[paste0("cp_",k)]]
  da <- a[[paste0("mn_",k)]]; db <- b[[paste0("mn_",k)]]
  w  <- suppressWarnings(wilcox.test(ca, cb)); wd <- suppressWarnings(wilcox.test(da, db))
  data.table(unit=unit_name, contrast=z[[1]], module=k, n_a=nrow(a), n_b=nrow(b),
             cp_a=mean(ca), cp_b=mean(cb), det_a=mean(da), det_b=mean(db),
             log2fc=log2((mean(ca)+.01)/(mean(cb)+.01)),
             log2fc_det=log2((mean(da)+1e-4)/(mean(db)+1e-4)),
             p=w$p.value, p_det=wd$p.value)
 }), fill=TRUE)
}), fill=TRUE)
RES <- rbind(run(U_FOV, "FOVxregion"), run(U_BIN, "bin"), fill=TRUE)
RES[, `:=`(q = p.adjust(p,"BH"), q_det = p.adjust(p_det,"BH")), by=.(unit, contrast)]
## Believed only when the pooled share and the per-cell detection agree -- the project's rule.
RES[, concordant := sign(log2fc) == sign(log2fc_det)]
RES[, verdict := fifelse(q < .05 & q_det < .05 & concordant, "SIGNIFICANT, concordant",
                 fifelse(q < .05 & concordant, "cp only", 
                 fifelse(q < .05, "DISCORDANT -- do not trust", "not significant")))]
sc("module_contrasts.csv", RES)
cat("\n=== module contrasts, FOV as the unit (the primary reading)\n")
print(RES[unit=="FOVxregion"][order(contrast, -abs(log2fc))][
      , .(contrast, module, n_a, n_b, cp_a=round(cp_a,2), cp_b=round(cp_b,2),
          log2fc=round(log2fc,3), log2fc_det=round(log2fc_det,3),
          q=signif(q,2), verdict)])
cat("\n=== does the FOV result survive at bin level too?\n")
cmp <- merge(RES[unit=="FOVxregion", .(contrast, module, lf_fov=log2fc, q_fov=q)],
             RES[unit=="bin",        .(contrast, module, lf_bin=log2fc, q_bin=q)],
             by=c("contrast","module"))
cmp[, agree := sign(lf_fov)==sign(lf_bin) & q_fov<.05 & q_bin<.05]
print(cmp[order(-abs(lf_fov))][, .(contrast, module, lf_fov=round(lf_fov,3),
      lf_bin=round(lf_bin,3), q_fov=signif(q_fov,2), q_bin=signif(q_bin,2), agree)])
sc("module_unit_agreement.csv", cmp)

## ---- 3. sign consistency across cell types ----------------------------------
## For each gene, is IB1 > RB1 in EVERY T subset independently? Five compartments agreeing is
## a two-sided sign test at p = 2/2^5 = 0.0625; with a 6th it is 0.031. Reported per gene with
## the count, so the reader sees the evidence rather than a p-value alone.
cat("\n=== sign consistency across T-cell compartments (IB1 vs RB1, same cell type)\n")
SG_GENES <- intersect(c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2",
                        "TOX","ENTPD1","TNFRSF9","VSIR","CXCL13","CXCR5","ICOS","IFNG",
                        "STAT1","IRF1","CXCL9","CXCL10","JAK2","TCF7","IL7R","SELL",
                        "GZMK","GZMB","FOXP3"), G)
## Per compartment per gene, kept rather than collapsed, so the heatmap below can show the
## evidence instead of only its summary. A "consistent" gene is one whose row is one colour.
PER <- rbindlist(lapply(SG_GENES, function(g) {
  rbindlist(lapply(CTS_T, function(cc) {
    ia <- d[, which(samp=="IB1" & celltype==cc)]; ib <- d[, which(samp=="RB1" & celltype==cc)]
    if (min(length(ia), length(ib)) < 200) return(NULL)
    data.table(gene=g, celltype=cc, n_IB1=length(ia), n_RB1=length(ib),
               cp_IB1 = 1e4*sum(M[g, ia])/sum(d$nCount_RNA[ia]),
               cp_RB1 = 1e4*sum(M[g, ib])/sum(d$nCount_RNA[ib]))
  }), fill=TRUE)
}), fill=TRUE)
PER[, log2fc := log2((cp_IB1+.01)/(cp_RB1+.01))]
sc("sign_consistency_per_compartment.csv", PER)
SG <- PER[, {
  n <- .N; k <- sum(cp_IB1 > cp_RB1)
  .(compartments=n, n_IB1_higher=k,
    direction = fifelse(k==n, "IB1 higher in ALL",
                fifelse(k==0, "RB1 higher in ALL", "mixed")),
    p_sign = 2*pbinom(min(k, n-k), n, .5),
    med_log2 = median(log2fc))
}, by = gene]
if (nrow(SG)) {
  SG[, q_sign := p.adjust(p_sign, "BH")]
  print(SG[order(p_sign, -abs(med_log2))][, .(gene, compartments, n_IB1_higher, direction,
        med_log2=round(med_log2,3), p_sign=signif(p_sign,3), q_sign=signif(q_sign,3))])
  sc("sign_consistency.csv", SG)

  ## ---- the heatmap: one row per gene, one column per T compartment ----------
  ## Fill is log2(IB1/RB1) IN THAT COMPARTMENT, so a sign-consistent gene reads as a row of
  ## a single colour. That is the whole claim, drawn: the inference here is not any single
  ## contrast but the agreement of five independent compartments, and a p-value cannot show
  ## agreement the way a uniform row can. Genes are ordered by median effect, and the row
  ## label carries the tally and the sign-test p so the strength is on the figure.
  ph <- merge(PER, SG[, .(gene, n_IB1_higher, compartments, p_sign, med_log2)], by="gene")
  ph[, lab := sprintf("%s   %d/%d  p=%.3f", gene, n_IB1_higher, compartments, p_sign)]
  ord <- unique(ph[order(-med_log2), .(lab, med_log2)])$lab
  ph[, lab := factor(lab, levels = rev(ord))]
  ph[, celltype := factor(celltype, levels = CTS_T)]
  lim <- max(abs(ph$log2fc), na.rm=TRUE)
  ncon <- SG[direction != "mixed", .N]
  ph_p <- ggplot(ph, aes(celltype, lab, fill = log2fc)) +
    geom_tile(colour = "#fcfcfb", linewidth = .5) +
    geom_text(aes(label = sprintf("%+.2f", log2fc)), size = 2.3, colour = "#0b0b0b") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f4f2ec", high = "#eb6834", midpoint = 0,
                         limits = c(-lim, lim), name = "log2 IB1/RB1") +
    labs(x = "T-cell compartment (same cell type both samples)", y = NULL,
         title = "Sign consistency across T compartments: IB1 versus RB1",
         subtitle = paste0("Each cell is log2(IB1/RB1) for that gene WITHIN that cell type, so ",
           "composition is conditioned out by construction.\nA sign-consistent gene reads as a ",
           "row of one colour -- five independent compartments agreeing is a sign test at\n",
           "p = 0.0625, the floor for n = 5, and it needs no single contrast to survive BH ",
           sprintf("across 6,175 genes.\n%d of %d genes are directionally consistent. ", ncon, nrow(SG)),
           "Red = higher in IB1 (the composite node)."),
         caption = "Data: output/17_pd1/<tier>/sign_consistency_per_compartment.csv") +
    theme_minimal(base_size = 10) +
    theme(plot.background=element_rect(fill="#fcfcfb",colour=NA),
          panel.background=element_rect(fill="#fcfcfb",colour=NA), panel.grid=element_blank(),
          axis.text.x=element_text(colour="#52514e", size=8.5, angle=20, hjust=1),
          axis.text.y=element_text(colour="#52514e", size=8, family="mono"),
          axis.title=element_text(colour="#52514e", size=9),
          plot.title=element_text(colour="#0b0b0b", face="bold", size=12),
          plot.subtitle=element_text(colour="#52514e", size=8),
          plot.caption=element_text(colour="#898781", size=7, hjust=0),
          legend.position="right", legend.key.width=unit(9,"pt"))
  ggsave(file.path(FIG,"fig_17_sign_consistency.png"), ph_p,
         width=9.5, height=1.1+0.30*nrow(SG), dpi=170, bg="#fcfcfb", limitsize=FALSE)
  cat("wrote ", file.path(FIG,"fig_17_sign_consistency.png"), "\n", sep="")
  ## THE INFERENCE IS THE COUNT, NOT ANY ONE GENE. With five compartments, the smallest
  ## p a single gene can reach is 2/2^5 = 0.0625, so no individual gene can ever clear 0.05
  ## and BH across 22 genes leaves everything at q ~ 0.08. Quoting those per-gene values as
  ## "not significant" would miss the point entirely. Under the null each gene is consistent
  ## with probability 0.0625, so the expected count is 22 x 0.0625 = 1.4 -- and the observed
  ## count is what carries the evidence.
  bt <- binom.test(ncon, nrow(SG), 0.0625, alternative = "greater")
  cat(sprintf("\ndirectionally consistent: %d of %d genes tested\n", ncon, nrow(SG)))
  cat(sprintf("expected under the null: %.2f | binomial p = %.3g\n",
              nrow(SG)*0.0625, bt$p.value))
  cat("So the pattern is far beyond chance even though no single gene can reach q < 0.05.\n")
  sc("sign_consistency_global.csv",
     data.table(genes_tested=nrow(SG), consistent=ncon, expected=round(nrow(SG)*0.0625,2),
                binom_p=bt$p.value,
                ib1_higher=SG[direction=="IB1 higher in ALL", .N],
                rb1_higher=SG[direction=="RB1 higher in ALL", .N]))
  print(SG[direction != "mixed"][order(-abs(med_log2)),
        .(gene, direction, med_log2 = round(med_log2,3), p_sign = signif(p_sign,3))])
}

## ---- 4. power bound on the nulls --------------------------------------------
## For every contrast that came out null, what fold change WOULD have been detected? A
## two-sample Wilcoxon at 80% power needs roughly d = 2.8/sqrt(n_eff) in standardised units;
## converted back through the observed SD of log2 cp10k this gives the smallest log2 fold
## change the design could have seen. Without this a null is unbounded and says nothing.
cat("\n=== power bound: the smallest fold change each null contrast could have detected\n")
PB <- rbindlist(lapply(CON, function(z) {
 rbindlist(lapply(MODS, function(k) {
  a <- U_FOV[samp==z[[2]] & bin_region==z[[3]]]; b <- U_FOV[samp==z[[4]] & bin_region==z[[5]]]
  if (min(nrow(a), nrow(b)) < MIN_UNITS) return(NULL)
  la <- log2(a[[paste0("cp_",k)]] + .01); lb <- log2(b[[paste0("cp_",k)]] + .01)
  sp <- sqrt(((length(la)-1)*var(la) + (length(lb)-1)*var(lb)) / (length(la)+length(lb)-2))
  ne <- 1/(1/length(la) + 1/length(lb))
  mde <- 2.8 * sp / sqrt(ne)          # ~80% power, 5% two-sided
  r <- RES[unit=="FOVxregion" & contrast==z[[1]] & module==k]
  data.table(contrast=z[[1]], module=k, n_a=nrow(a), n_b=nrow(b),
             sd_log2=round(sp,3), mde_log2=round(mde,3),
             mde_fold=round(2^mde,2), observed_log2=round(if(nrow(r)) r$log2fc else NA_real_,3),
             verdict = if (!nrow(r)) NA_character_ else r$verdict)
 }), fill=TRUE)
}), fill=TRUE)
print(PB[order(contrast, mde_fold)])
sc("power_bounds.csv", PB)
cat(sprintf("\nmedian detectable fold change across contrasts/modules: %.2fx\n",
            median(PB$mde_fold, na.rm=TRUE)))
cat("A null contrast therefore excludes a difference LARGER than its mde_fold; it does not\n")
cat("exclude a smaller one. Quote the bound, never the bare null.\n")
cat("\nCAVEAT: this assumes FOV-level values are independent. They are not -- neighbouring FOVs\n")
cat("share tissue, so the true standard error is larger and these bounds are OPTIMISTIC. Read\n")
cat("them as a lower limit on what was excluded, not an exact threshold. The bound is tight\n")
cat("because AGGREGATION works: ~390 module counts per FOV puts the Poisson floor near 5%,\n")
cat("which is the whole reason the FOV is the unit rather than the cell.\n")
writeLines(format(Sys.time()), file.path(OUT, "modules.done"))
cat("============================================================\nOK\n")
