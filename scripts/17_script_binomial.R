## ============================================================================
## 17_script_binomial.R -- does region predict PD-1/PD-L1 positivity, once depth is modelled?
##
## WHAT THIS REPLACES AND WHY. Step 16 tested this by DEPTH-MATCHING regions and then running a
## Wilcoxon. That works but throws away most of the data -- matching IB1 FL to RB1 follicular
## bins kept 27 of 80 groups -- and it removes the confound rather than measuring it. A
## binomial model keeps every cell and puts depth in the model instead.
##
## THE OUTCOME is per-cell positivity at >= 2 counts. That threshold is not a choice: at >= 1
## count PDCD1/CD274/PDCD1LG2 sit 2.6-2.9x the mean blank probe, i.e. near noise, and at >= 2
## the blank rate collapses to 0.05% and they are 25-31x it (qc_16_negprobe_floor.R). Below 2 we
## would be modelling background.
##
## DEPTH ENTERS AS A DECILE STRATUM, not as a linear term in log counts. Assuming log-linearity
## would be an untested modelling choice, and depth predicts positivity at AUC ~0.70 for EVERY
## gene here including IGKC and CD3E, so getting it wrong matters. Deciles are agnostic about
## the shape. It also makes the fit cheap: cells collapse to (FOV x region x decile) counts and
## a binomial GLM on those aggregated counts is mathematically identical to the per-cell fit,
## which turns a 2-million-row regression into a few thousand rows.
##
## FOV IS THE CLUSTER. Cells in one field of view are not independent observations, so the
## standard errors come from resampling FOVS with replacement, not cells. `sandwich` is not
## installed and no install was made for this; a cluster bootstrap needs no package and makes
## the same correction. lme4 IS available, so a mixed-model cross-check is run on the primary
## contrasts to confirm the bootstrap intervals are not wildly off.
##
## THE OVERLAP DIAGNOSTIC IS REPORTED FIRST, because it decides whether the model can answer
## the question at all: if two regions barely share depth deciles then region and depth are
## collinear and no amount of modelling separates them. That is the same identifiability issue
## 16c found for composition, in a different guise.
##
## Usage: Rscript scripts/17_script_binomial.R [tier]
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

THR    <- 2L
N_DEC  <- 10L
N_BOOT <- 400L
MIN_FOV<- 12L
MIN_POS<- 30L
SEED   <- 1L
GOI    <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","TOX","CXCL13","CXCL9")
STRATA <- list(all = NULL, B = "B", T = c("Tfh","Treg","CD4_conv","CD8_T","T_unresolved"),
               myeloid = c("Macrophage","Monocyte","cDC","myeloid_unresolved"))
FOC <- list(IB1 = c("FL_focus","cHL_hotspot","background"),
            RB1 = c("FL_follicular","interfollicular"))
CON <- list(
  list("IB1 FL vs IB1 bg",      "IB1","FL_focus",      "IB1","background"),
  list("IB1 hot vs IB1 bg",     "IB1","cHL_hotspot",   "IB1","background"),
  list("RB1 foll vs RB1 inter", "RB1","FL_follicular", "RB1","interfollicular"),
  list("IB1 FL vs RB1 foll",    "IB1","FL_focus",      "RB1","FL_follicular"),
  list("IB1 bg vs RB1 inter",   "IB1","background",    "RB1","interfollicular"))
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"

cat("============================================================\n")
cat("step 17 -- binomial model of per-cell positivity, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

## d_hot_um is needed to clean the reference: field_extent_verdict.csv shows the Hodgkin
## chemokine field decays 2.71x [1.95, 3.91] out to ~400-800um, so IB1 `background` within
## 800um of a hotspot is NOT disease-free and must not act as reference tissue.
rg <- fread(REG_F, select=c("cell_id","samp","bin_region","is_defect","unbinned",
                            "fov_id","nCount_RNA","d_hot_um"))
ct <- fread(CT_F, select=c("cell_id","celltype"))
d  <- merge(rg, ct, by="cell_id", sort=FALSE)[is_defect==FALSE & unbinned==FALSE]
d  <- d[mapply(function(s,r) r %in% FOC[[s]], samp, bin_region)]
cat(sprintf("cells %d\n", nrow(d)))

CNT <- NULL
for (S in c("IB1","RB1")) {
  f <- sprintf(PRE_RDS, S); cat("reading ", f, "\n", sep="")
  o <- readRDS(f); cn <- LayerData(o[["RNA"]], layer="counts"); rn <- colnames(cn)
  g <- intersect(GOI, rownames(cn)); sub <- d[samp==S]
  jj <- match(sub$cell_id, rn); stopifnot(!anyNA(jj))
  M <- as(cn[g, jj, drop=FALSE], "dgCMatrix"); colnames(M) <- sub$cell_id
  CNT[[S]] <- M; rm(o, cn); invisible(gc())
}
G <- Reduce(intersect, lapply(CNT, rownames))
M <- do.call(cbind, lapply(CNT, function(x) x[G,,drop=FALSE])); rm(CNT); invisible(gc())
M <- M[, d$cell_id, drop=FALSE]
cat(sprintf("genes %d x cells %d | threshold >= %d counts\n", nrow(M), ncol(M), THR))
for (g in G) d[, (paste0("pos_", g)) := as.vector(M[g, ]) >= THR]
cat("\npositivity rate per gene:\n")
print(d[, lapply(.SD, function(v) round(100*mean(v),3)), .SDcols=paste0("pos_",G)])

## ---- the model ---------------------------------------------------------------
## Aggregated binomial: rows are (fov, region, depth decile), response is (n_pos, n_neg).
## Cluster bootstrap resamples FOVS, so a whole field's cells move together.
fit_one <- function(g, z, lab_a, lab_b, stratum) {
  z <- z[!is.na(get(paste0("pos_", g)))]
  if (!nrow(z)) return(NULL)
  z[, dec := cut(nCount_RNA, quantile(nCount_RNA, seq(0,1,length.out=N_DEC+1)),
                 include.lowest=TRUE, labels=FALSE)]
  agg <- z[, .(n = .N, k = sum(get(paste0("pos_", g)))),
           by = .(fov_key = paste(samp, fov_id, sep="_"), grp = bin_region, dec)]
  ## OVERLAP: deciles in which BOTH regions appear. Without overlap, region and depth are
  ## collinear and the coefficient is not identified -- report, do not fit.
  ov <- agg[, .(regions = uniqueN(grp)), by = dec][regions == 2L, .N]
  nfa <- agg[grp == lab_a, uniqueN(fov_key)]; nfb <- agg[grp == lab_b, uniqueN(fov_key)]
  npos <- agg[, sum(k)]
  if (ov < 3L || min(nfa, nfb) < MIN_FOV || npos < MIN_POS)
    return(data.table(gene=g, stratum=stratum, contrast=paste(lab_a,"vs",lab_b),
                      deciles_overlapping=ov, n_fov_a=nfa, n_fov_b=nfb, n_pos=npos,
                      or=NA_real_, lo=NA_real_, hi=NA_real_, p=NA_real_,
                      verdict="not identified (no depth overlap or too few FOVs/positives)"))
  agg[, grp := factor(grp, levels=c(lab_b, lab_a))]   # reference = b
  agg[, dec := factor(dec)]
  cf <- function(a) {
    m <- suppressWarnings(glm(cbind(k, n-k) ~ grp + dec, family=binomial, data=a))
    unname(coef(m)[2])
  }
  est <- tryCatch(cf(agg), error=function(e) NA_real_)
  if (!is.finite(est)) return(NULL)
  fovs <- unique(agg$fov_key)
  set.seed(SEED)
  bs <- vapply(seq_len(N_BOOT), function(i) {
    fs <- sample(fovs, length(fovs), replace=TRUE)
    a  <- agg[data.table(fov_key=fs), on="fov_key", allow.cartesian=TRUE, nomatch=0L]
    if (uniqueN(a$grp) < 2L) return(NA_real_)
    tryCatch(cf(a), error=function(e) NA_real_)
  }, 0)
  bs <- bs[is.finite(bs)]
  if (length(bs) < 50) return(NULL)
  se <- sd(bs)
  data.table(gene=g, stratum=stratum, contrast=paste(lab_a,"vs",lab_b),
             deciles_overlapping=ov, n_fov_a=nfa, n_fov_b=nfb, n_pos=npos,
             or=exp(est), lo=exp(quantile(bs,.025)), hi=exp(quantile(bs,.975)),
             p=2*pnorm(-abs(est/se)), boot_ok=length(bs))
}

## The fit is the expensive part (~27 min). It is cached to binomial_model_raw.Rds the
## moment it completes, so a crash in the *reporting* code below never costs the fit again.
## Set REFIT=1 in the environment to force a refit.
RAW <- file.path(OUT, "binomial_model_raw.Rds")
if (file.exists(RAW) && !nzchar(Sys.getenv("REFIT"))) {
  cat("\n=== loading cached fit from", RAW, "-- set REFIT=1 to refit\n")
  RES <- readRDS(RAW)
} else {
## Two reference definitions, fitted side by side rather than one replacing the other, so the
## effect of cleaning the reference is visible instead of asserted. `far` drops IB1 background
## cells within FIELD_UM of a hotspot; every other group is untouched, and RB1 has no hotspots
## so it is unaffected either way.
FIELD_UM <- 800
REFDEFS <- list(all_background = NULL, far_background = FIELD_UM)

RES <- rbindlist(lapply(names(REFDEFS), function(rdf) {
 cut_um <- REFDEFS[[rdf]]
 dd <- if (is.null(cut_um)) d else
       d[!(samp == "IB1" & bin_region == "background" & d_hot_um <= cut_um)]
 cat(sprintf("\n########## reference definition '%s': %d cells (IB1 background %d)\n",
             rdf, nrow(dd), dd[samp=="IB1" & bin_region=="background", .N]))
 out <- rbindlist(lapply(names(STRATA), function(st) {
  keep <- if (is.null(STRATA[[st]])) rep(TRUE, nrow(dd)) else dd$celltype %in% STRATA[[st]]
  dz <- dd[keep]
  cat(sprintf("\n--- stratum '%s': %d cells\n", st, nrow(dz)))
  rbindlist(lapply(CON, function(z) {
    zz <- dz[(samp==z[[2]] & bin_region==z[[3]]) | (samp==z[[4]] & bin_region==z[[5]])]
    if (!nrow(zz)) return(NULL)
    rbindlist(lapply(G, function(g) fit_one(g, zz, z[[3]], z[[5]], st)), fill=TRUE)
  }), fill=TRUE)
 }), fill=TRUE)
 if (nrow(out)) out[, refdef := rdf]
 out
}), fill=TRUE)
RES[, q := p.adjust(p, "BH"), by=.(refdef, stratum, contrast)]
## `verdict` only exists if at least one contrast returned the not-identified stub, and
## rbindlist(fill=TRUE) cannot invent a column no row supplied. Every contrast being
## identifiable therefore DELETED the column and the next line failed -- after 26 minutes of
## completed fitting. Create it unconditionally instead of assuming a branch was taken.
if (!"verdict" %in% names(RES)) RES[, verdict := NA_character_]
RES[, verdict := fifelse(!is.na(verdict), verdict,
                 fifelse(is.na(q), "not fitted",
                 fifelse(q < .05, "significant", "not significant")))]
saveRDS(RES, RAW)   # so reformatting never refits
}
sc("binomial_model.csv", RES)

cat("\n=== identifiability first: contrasts where depth and region are collinear\n")
ni <- RES[grepl("not identified", verdict)]
if (nrow(ni)) print(unique(ni[, .(stratum, contrast, deciles_overlapping, n_fov_a, n_fov_b)])) else
  cat("  none -- every contrast had depth overlap and enough FOVs\n")

cat("\n=== the PD-1 axis, odds ratio for region with depth stratified out\n")
print(RES[gene %in% c("PDCD1","CD274","PDCD1LG2") & !is.na(or)][order(refdef, stratum, contrast, gene)][
  , .(refdef, stratum, contrast, gene, or=round(or,3),
      ci=sprintf("[%.2f, %.2f]", lo, hi), q=signif(q,2))])

cat("\n=== what cleaning the reference changed (contrasts that USE IB1 background)\n")
cmpk <- c("stratum","contrast","gene")
a <- RES[refdef=="all_background", c(cmpk, "or","q"), with=FALSE]
b <- RES[refdef=="far_background", c(cmpk, "or","q"), with=FALSE]
cp <- merge(a, b, by=cmpk, suffixes=c("_all","_far"))
cp <- cp[grepl("bg|background", contrast) & is.finite(or_all) & is.finite(or_far)]
cp[, shift := or_far/or_all]
print(cp[order(-abs(log(shift)))][1:25, .(stratum, contrast, gene,
      or_all=round(or_all,3), or_far=round(or_far,3), shift=round(shift,3),
      q_all=signif(q_all,2), q_far=signif(q_far,2))])

cat("\n=== everything significant\n")
sg <- RES[!is.na(q) & q < .05][order(-abs(log(or)))]
if (nrow(sg)) print(sg[, .(refdef, stratum, contrast, gene, or=round(or,3),
      ci=sprintf("[%.2f, %.2f]", lo, hi), q=signif(q,2))]) else cat("  none\n")

## ---- cross-check the bootstrap against a mixed model ------------------------
## lme4 is installed, so the primary contrasts get a glmer with FOV as a random intercept.
## If the two disagree badly the cluster bootstrap is not doing its job.
if (requireNamespace("lme4", quietly=TRUE)) {
  cat("\n=== glmer cross-check (FOV random intercept) on the PD-1 axis, stratum 'all'\n")
  ## Braces, not a bare expression, for the outer lapply body. Written inline as
  ## `lapply(CON, function(z) rbindlist(lapply(...), fill=TRUE), fill=TRUE)` the second
  ## fill= binds to lapply, not to the outer rbindlist, and R errors with
  ## "unused argument (fill = TRUE)". This is the fourth site in this project.
  cc <- rbindlist(lapply(CON, function(z) {
    rbindlist(lapply(c("PDCD1","CD274"), function(g) {
    zz <- d[(samp==z[[2]] & bin_region==z[[3]]) | (samp==z[[4]] & bin_region==z[[5]])]
    if (!nrow(zz)) return(NULL)
    zz[, dec := cut(nCount_RNA, quantile(nCount_RNA, seq(0,1,length.out=N_DEC+1)),
                    include.lowest=TRUE, labels=FALSE)]
    a <- zz[, .(n=.N, k=sum(get(paste0("pos_",g)))),
            by=.(fov_key=paste(samp,fov_id,sep="_"), grp=bin_region, dec=factor(dec))]
    if (uniqueN(a$grp) < 2L || a[, sum(k)] < MIN_POS) return(NULL)
    a[, grp := factor(grp, levels=c(z[[5]], z[[3]]))]
    m <- tryCatch(lme4::glmer(cbind(k, n-k) ~ grp + dec + (1|fov_key),
                              family=binomial, data=a,
                              control=lme4::glmerControl(optimizer="bobyqa")),
                  error=function(e) NULL, warning=function(w) NULL)
    if (is.null(m)) return(NULL)
    s <- summary(m)$coefficients
    data.table(contrast=paste(z[[3]],"vs",z[[5]]), gene=g,
               or_glmer=exp(s[2,1]), p_glmer=s[2,4])
    }), fill=TRUE)
  }), fill=TRUE)
  if (nrow(cc)) {
    cmp <- merge(RES[stratum=="all" & refdef=="all_background", .(contrast, gene, or_boot=or, q_boot=q)], cc,
                 by=c("contrast","gene"))
    cmp[, agree := is.finite(or_boot) & sign(log(or_boot))==sign(log(or_glmer))]
    print(cmp[, .(contrast, gene, or_boot=round(or_boot,3), or_glmer=round(or_glmer,3),
                  q_boot=signif(q_boot,2), p_glmer=signif(p_glmer,2), agree)])
    sc("binomial_glmer_check.csv", cmp)
  } else cat("  glmer did not converge on any contrast; bootstrap stands unchecked\n")
}
writeLines(format(Sys.time()), file.path(OUT, "binomial.done"))
cat("============================================================\nOK\n")
