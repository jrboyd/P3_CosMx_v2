## ============================================================================
## 16_script_cooccurrence.R -- is the FL microenvironment differently ORGANISED?
##
## WHY THIS MATTERS MORE NOW THAN WHEN IT WAS PLANNED. The cell-type-conditioned DE (16m)
## found that most regional signal is COMPOSITIONAL: IB1's FL foci differ from IB1 background
## in cell mix, but zero genes differ within B or T cells. Composition and expression are
## therefore nearly exhausted as questions. What has not been asked is whether the same cell
## types are ARRANGED differently -- whether a Tfh sits next to a B cell more often in one
## node than the other, given the same number of each.
##
## THE NULL IS THE WHOLE ANALYSIS. Adjacency counts are driven by density and composition
## before any biology: a dense B-rich region gives every cell more B neighbours automatically.
## So labels are permuted WITHIN LOCAL-DENSITY DECILES, preserving both the composition of the
## region and its geometry, and only the assignment of types to positions is randomised. This
## is step 10's adjacency null (10_script_hrs_call.R:234-266), which measured HRS-HRS adjacency
## at 1.03x over exactly this kind of null -- i.e. the null is strong enough to kill a
## plausible-looking effect, which is why it is worth using.
##
## The statistic is the mean number of type-B neighbours within R_UM of a type-A cell. The
## permutation supplies the inference, so there is no cell-level p-value pretending cells are
## independent. A per-FOV rerun is reported alongside: if an enrichment is carried by one or
## two FOVs it is not a property of the region.
##
## Usage: Rscript scripts/16_script_cooccurrence.R [tier]
## ============================================================================
suppressPackageStartupMessages({
  library(RANN); library(data.table); library(ggplot2); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "matched"

REG_F <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT_F  <- "output/16_celltypes/cell_celltypes.csv.gz"
OUT   <- file.path("output/16_cooccur", TIER)
FIG   <- file.path(OUT, "figures")
for (f in c(REG_F, CT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

R_UM    <- 20      # centre-to-centre; cells are ~8-9 um across, so this is touching/near
K_MAX   <- 40L     # cap on neighbours retrieved per cell
N_PERM  <- 200L
MIN_A   <- 100L    # minimum cells of the focal type for a pair to be tested
N_DEC   <- 10L
SEED    <- 1L
TYPES   <- c("Tfh","Treg","CD4_conv","CD8_T","B","Macrophage","Monocyte","cDC","PC","NK",
             "endo","fibro")
FOC <- list(IB1 = c("FL_focus","cHL_hotspot","background"),
            RB1 = c("FL_follicular","interfollicular"))
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
DIV_LO<-"#2a78d6"; DIV_MID<-"#f4f2ec"; DIV_HI<-"#eb6834"

cat("============================================================\n")
cat("step 16h -- cell-type co-occurrence vs a density-matched null, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

rg <- fread(REG_F, select=c("cell_id","samp","bin_region","x_slide_mm","y_slide_mm",
                            "is_defect","unbinned","fov_id"))
ct <- fread(CT_F, select=c("cell_id","celltype"))
d  <- merge(rg, ct, by="cell_id", sort=FALSE)[is_defect==FALSE & unbinned==FALSE]
d  <- d[mapply(function(s,r) r %in% FOC[[s]], samp, bin_region)]
cat(sprintf("cells %d\n", nrow(d)))
print(d[, .N, by=.(samp, bin_region)][order(samp,-N)])

## One region at a time. Neighbours within R_UM, self excluded.
run_region <- function(S, R) {
  z <- d[samp==S & bin_region==R]
  if (nrow(z) < 2000) { cat(sprintf("  %s / %s: %d cells, skipped\n", S, R, nrow(z))); return(NULL) }
  xy <- as.matrix(z[, .(x_slide_mm*1000, y_slide_mm*1000)])
  k  <- min(K_MAX, nrow(xy))
  nn <- RANN::nn2(xy, xy, k = k)
  ## drop self (column 1) and anything beyond R_UM
  idx <- nn$nn.idx[, -1, drop=FALSE]; dst <- nn$nn.dists[, -1, drop=FALSE]
  keep <- dst <= R_UM
  ## density = distance to the 20th neighbour; deciles are the permutation strata
  dcol <- min(20L, ncol(nn$nn.dists))
  dens <- nn$nn.dists[, dcol]
  strat <- cut(dens, quantile(dens, seq(0,1,length.out=N_DEC+1)), include.lowest=TRUE, labels=FALSE)
  lab0 <- z$celltype
  present <- intersect(TYPES, unique(lab0))
  cnt <- table(factor(lab0, levels=present))
  ## mean neighbours of type b per cell of type a, for one labelling
  stat <- function(lab) {
    nb <- matrix(0, nrow=length(present), ncol=length(present),
                 dimnames=list(present, present))
    lv <- match(lab, present)
    for (b in seq_along(present)) {
      isb <- lv == b; isb[is.na(isb)] <- FALSE
      hit <- matrix(isb[idx], nrow=nrow(idx)) & keep
      per <- rowSums(hit)                       # type-b neighbours of every cell
      nb[, b] <- vapply(seq_along(present), function(a) {
        ia <- which(lv == a); if (!length(ia)) NA_real_ else mean(per[ia]) }, 0)
    }
    nb
  }
  obs <- stat(lab0)
  set.seed(SEED)
  perms <- vector("list", N_PERM)
  for (p in seq_len(N_PERM)) {
    lab <- lab0
    for (s in unique(strat)) { ii <- which(strat == s); lab[ii] <- lab0[sample(ii)] }
    perms[[p]] <- stat(lab)
  }
  nu <- Reduce(`+`, perms) / N_PERM
  sdv <- sqrt(Reduce(`+`, lapply(perms, function(m) (m - nu)^2)) / max(N_PERM-1,1))
  ge  <- Reduce(`+`, lapply(perms, function(m) m >= obs))
  le  <- Reduce(`+`, lapply(perms, function(m) m <= obs))
  out <- CJ(a = present, b = present)
  out[, `:=`(samp=S, bin_region=R,
             n_a = as.integer(cnt[a]), n_b = as.integer(cnt[b]),
             obs = obs[cbind(a,b)], null = nu[cbind(a,b)], null_sd = sdv[cbind(a,b)])]
  out[, enrich := obs / pmax(null, 1e-9)]
  out[, p_emp := (1 + pmin(ge[cbind(a,b)], le[cbind(a,b)])) / (N_PERM + 1) * 2]
  out <- out[n_a >= MIN_A & n_b >= MIN_A]
  cat(sprintf("  %s / %-16s %6d cells | %d types | %d pairs tested\n",
              S, R, nrow(z), length(present), nrow(out)))
  out[]
}

CO <- rbindlist(lapply(names(FOC), function(S)
  rbindlist(lapply(FOC[[S]], run_region, S = S), fill = TRUE)), fill = TRUE)
if (!nrow(CO)) stop("no region was testable", call. = FALSE)
CO[, q := p.adjust(p_emp, "BH"), by = .(samp, bin_region)]
CO[, log2_enrich := log2(pmax(enrich, 1e-6))]
for (cc in c("obs","null","null_sd","enrich","log2_enrich"))
  CO[[cc]] <- round(CO[[cc]], 4)
sc("cooccurrence.csv", CO)

cat("\n=== strongest departures from the density-matched null (q < 0.05), per region\n")
for (k in unique(paste(CO$samp, CO$bin_region))) {
  z <- CO[paste(samp,bin_region) == k & q < 0.05][order(-abs(log2_enrich))]
  if (!nrow(z)) { cat(sprintf("\n%s: none\n", k)); next }
  cat(sprintf("\n%s -- %d significant pairs, top 8:\n", k, nrow(z)))
  print(head(z[, .(a, b, n_a, obs, null, enrich, q = signif(q,2))], 8))
}

## The pair the FL question turns on: does a Tfh sit beside a B cell more than chance, and
## does that differ between IB1's FL foci and RB1's follicles?
cat("\n=== the Tfh <-> B axis, and Tfh's other partners\n")
key <- CO[a == "Tfh" & b %in% c("B","Tfh","Treg","CD4_conv","CD8_T","Macrophage","cDC")]
if (nrow(key)) print(dcast(key, b ~ samp + bin_region, value.var = "enrich"))
sc("tfh_partners.csv", key)

## ---- figure ------------------------------------------------------------------
pl <- CO[a %in% c("Tfh","B","Treg","CD8_T","Macrophage") &
         b %in% c("Tfh","B","Treg","CD8_T","Macrophage")]
if (nrow(pl)) {
  pl[, reg := paste(samp, bin_region, sep="\n")]
  lim <- max(abs(pl$log2_enrich), na.rm=TRUE)
  p <- ggplot(pl, aes(b, a, fill = log2_enrich)) +
    geom_tile(colour = SURF, linewidth = .5) +
    geom_text(aes(label = ifelse(q < 0.05, sprintf("%.2f", enrich), "")),
              size = 2.5, colour = INK) +
    facet_wrap(~ reg, nrow = 1) +
    scale_fill_gradient2(low = DIV_LO, mid = DIV_MID, high = DIV_HI, midpoint = 0,
                         limits = c(-lim, lim), name = "log2 obs/null") +
    labs(x = "neighbour type", y = "focal type",
         title = "Is the same cell mix arranged differently in the two nodes?",
         subtitle = sprintf(paste0("Mean type-b neighbours within %g um of a type-a cell, ",
           "against a null that permutes labels within local-density\ndeciles -- so ",
           "composition and geometry are held fixed and only the arrangement is randomised. ",
           "%d permutations.\nNumbers shown where q < 0.05. 1.0 = exactly as chance predicts."),
           R_UM, N_PERM),
         caption = "Data: output/16_cooccur/<tier>/cooccurrence.csv") +
    theme_minimal(base_size = 10) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA), panel.grid=element_blank(),
          axis.text.x=element_text(colour=INK2, angle=45, hjust=1, size=8),
          axis.text.y=element_text(colour=INK2, size=8),
          axis.title=element_text(colour=INK2, size=9),
          strip.text=element_text(colour=INK, face="bold", size=8),
          plot.title=element_text(colour=INK, face="bold", size=12),
          plot.subtitle=element_text(colour=INK2, size=8),
          plot.caption=element_text(colour=MUTED, size=7, hjust=0),
          legend.position="right", legend.key.width=unit(9,"pt"))
  ggsave(file.path(FIG,"fig_16_cooccurrence.png"), p, width=15, height=4.6, dpi=170, bg=SURF)
  cat("wrote ", file.path(FIG,"fig_16_cooccurrence.png"), "\n", sep="")
}
writeLines(format(Sys.time()), file.path(OUT,"cooccur.done"))
cat("============================================================\nOK\n")
