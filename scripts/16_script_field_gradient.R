## ============================================================================
## 16_script_field_gradient.R -- does the Hodgkin chemokine field change the FL beside it?
##
## THIS IS THE ONE CROSS-DISEASE QUESTION THAT SURVIVES EVERYTHING. It runs entirely inside
## IB1, so it is immune to all three problems that broke the cross-sample contrasts: the
## IB1-only staining QC, the depth asymmetry, and the composition disjointness (16c).
##
## THE HYPOTHESIS IS A FIELD, NOT CO-LOCATION, and that is not a hedge -- step 12 already
## tested co-location and found none (hotspot-to-nearest-FL-focus distance, p = 0.23). The
## chemokines are SECRETED, step 11 found ~27% of IB1 bins clear RB1's 99th percentile, and
## step 11 called that a tissue-wide field rather than a niche. A graded effect is therefore
## the shape the data actually supports.
##
## TWO EXPOSURE AXES, because they can disagree and the disagreement is informative:
##   distance  to the nearest cHL hotspot bin (geometric)
##   field     the bin's own CCL17+CCL22 cp10k (local concentration)
## Distance is the cleaner causal story; the field level is the better dose measure. A result
## that holds on both is worth more than one that holds on either.
##
## CCL17 and CCL22 SIGNAL THROUGH CCR4, which recruits Treg and Th2 cells. So the pre-declared
## readout is CCR4/FOXP3/CTLA4/PDCD1, and the genome-wide scan is reported second -- declared
## in this order so the CCR4 result cannot be a post-hoc pick from 6,175 genes.
##
## CIRCULARITY: CCL17/CCL22 define the exposure, so they are excluded and flagged. Hotspot
## bins themselves are EXCLUDED from the response set -- the question is about the tissue
## around a hotspot, not the hotspot.
##
## Usage: Rscript scripts/16_script_field_gradient.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(ggplot2); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "matched"

PB   <- file.path("output/16_pseudobulk", TIER)
REG  <- file.path("output/16_regions", TIER)
OUT  <- file.path("output/16_gradient", TIER)
FIG  <- file.path(OUT, "figures")
for (f in c(file.path(PB,"pb_sum.Rds"), file.path(PB,"bin_meta.csv.gz"),
            file.path(REG,"bins_all_scales.csv.gz")))
  if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

RINGS   <- c(0, 200, 400, 800, 1600, Inf)   # um from the nearest hotspot
MIN_BIN <- 12L; Q_CUT <- 0.05; SEED <- 1L
PRIMARY <- c("CCR4","FOXP3","CTLA4","PDCD1","IL2RA","IKZF2","CCL17","CCL22")
SEL     <- c("CCL17","CCL22")
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
DIV_HI<-"#2a78d6"; DIV_LO<-"#eb6834"

cat("============================================================\n")
cat("step 16g -- Hodgkin field gradient inside IB1, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

S  <- readRDS(file.path(PB,"pb_sum.Rds"))
BM <- fread(file.path(PB,"bin_meta.csv.gz")); setkey(BM, bin); BM <- BM[colnames(S)]
stopifnot(identical(BM$bin, colnames(S)))
B2 <- fread(file.path(REG,"bins_all_scales.csv.gz"))[bin_um == 100L]
B2[, bin := paste(samp, bx, by_, sep="_")]
BM <- merge(BM, B2[, .(bin, chemo_cp10k, x_mm, y_mm, niche, fb_a, f_assigned, arch)],
            by = "bin", sort = FALSE)
GENES <- rownames(S)
CP <- as.matrix(S) %*% Diagonal(x = 1e4/pmax(Matrix::colSums(S),1)); dimnames(CP) <- dimnames(S)

## ---- exposure ---------------------------------------------------------------
ib <- BM[samp == "IB1"]
hot <- ib[niche == TRUE]
if (nrow(hot) < 3) stop("too few hotspot bins", call. = FALSE)
ii <- match(ib$bin, colnames(CP))
ib[, d_hot_um := RANN::nn2(as.matrix(hot[, .(x_mm*1000, y_mm*1000)]),
                           as.matrix(ib[, .(x_mm*1000, y_mm*1000)]), k = 1)$nn.dists[,1]]
## Hotspot bins are the exposure source, not the response.
resp <- ib[niche == FALSE]
resp[, ring := cut(d_hot_um, RINGS, labels = c("0-200","200-400","400-800","800-1600",">1600"),
                   right = FALSE)]
resp[, field_q := cut(chemo_cp10k, quantile(chemo_cp10k, c(0,.25,.5,.75,1)),
                      labels = c("Q1 low","Q2","Q3","Q4 high"), include.lowest = TRUE)]
cat(sprintf("IB1 bins %d | hotspots %d | response bins %d\n", nrow(ib), nrow(hot), nrow(resp)))
cat("\n=== bins per distance ring, and their composition (is the gradient confounded?)\n")
rt <- resp[, .(bins = .N, fovs = uniqueN(fov), med_d = as.numeric(median(d_hot_um)),
               med_chemo = round(median(chemo_cp10k),1), med_fb_a = round(median(fb_a,na.rm=TRUE),3),
               med_depth = round(median(med_depth)), med_cells = as.numeric(median(cells)),
               pct_FL = round(100*mean(bin_region=="FL_focus"),1)), by = ring][order(ring)]
print(rt); sc("ring_summary.csv", rt)
cat("\n=== bins per field quartile\n")
ft <- resp[, .(bins = .N, fovs = uniqueN(fov), med_chemo = round(median(chemo_cp10k),1),
               med_fb_a = round(median(fb_a,na.rm=TRUE),3), med_depth = round(median(med_depth)),
               med_d_hot = round(median(d_hot_um))), by = field_q][order(field_q)]
print(ft); sc("field_quartile_summary.csv", ft)

## ---- the pre-declared CCR4 axis ---------------------------------------------
## Spearman across bins, then across FOV means as the conservative reading. Reported for the
## whole response set and for FL_focus bins alone -- the latter is the actual question ("does
## the field change the FL"), the former says whether it is FL-specific.
grad1 <- function(sub, lab) {
  if (nrow(sub) < MIN_BIN) return(NULL)
  jj <- match(sub$bin, colnames(CP))
  gg <- intersect(PRIMARY, GENES)
  rbindlist(lapply(gg, function(g) {
    v <- CP[g, jj]
    ct <- suppressWarnings(cor.test(v, sub$chemo_cp10k, method = "spearman"))
    cd <- suppressWarnings(cor.test(v, sub$d_hot_um,    method = "spearman"))
    fv <- sub[, .(f = mean(v[.I]), c = mean(chemo_cp10k), d = mean(d_hot_um)), by = fov]
    cf <- if (nrow(fv) >= MIN_BIN)
      suppressWarnings(cor.test(fv$f, fv$c, method="spearman")) else list(estimate=NA,p.value=NA)
    data.table(set = lab, gene = g, bins = nrow(sub), fovs = uniqueN(sub$fov),
               rho_field = unname(ct$estimate), p_field = ct$p.value,
               rho_dist  = unname(cd$estimate), p_dist  = cd$p.value,
               rho_field_fov = unname(cf$estimate), p_field_fov = cf$p.value,
               is_selection = g %in% SEL)
  }))
}
G1 <- rbindlist(list(grad1(resp, "all non-hotspot bins"),
                     grad1(resp[bin_region == "FL_focus"], "FL_focus bins only")), fill = TRUE)
if (nrow(G1)) {
  G1[, `:=`(q_field = p.adjust(p_field,"BH"), q_dist = p.adjust(p_dist,"BH")), by = set]
  cat("\n=== PRE-DECLARED CCR4 axis: correlation with the local chemokine field\n")
  cat("    (positive rho_field = higher where the field is higher; rho_dist should be the mirror)\n")
  print(G1[order(set, -abs(rho_field)),
           .(set, gene, bins, fovs, rho_field = round(rho_field,3), q_field = signif(q_field,2),
             rho_dist = round(rho_dist,3), rho_field_fov = round(rho_field_fov,3),
             is_selection)])
  sc("ccr4_axis.csv", G1)
  ok <- G1[set == "all non-hotspot bins" & gene == "CCR4"]
  if (nrow(ok)) cat(sprintf("\n  CCR4: rho %.3f (q %.3g) across bins, %.3f across FOV means -> %s\n",
    ok$rho_field, ok$q_field, ok$rho_field_fov,
    if (is.finite(ok$q_field) && ok$q_field < Q_CUT && sign(ok$rho_field) == sign(ok$rho_field_fov))
      "field-associated at both units" else "NOT robust across units"))
}

## ---- genome-wide, reported second -------------------------------------------
jj <- match(resp$bin, colnames(CP))
rho <- suppressWarnings(apply(CP[, jj, drop=FALSE], 1, function(v)
  cor(v, resp$chemo_cp10k, method="spearman")))
n <- nrow(resp); z <- atanh(pmax(pmin(rho, 1-1e-9), -1+1e-9)) * sqrt(n-3)
GW <- data.table(gene = GENES, rho_field = rho, p = 2*pnorm(-abs(z)))
GW[, q := p.adjust(p,"BH")][, is_selection := gene %in% SEL]
## FOV-collapsed rerun, the conservative unit.
fv <- resp[, .(rows = .I, fov)][, .(k = list(rows)), by = fov]
cm <- do.call(cbind, lapply(fv$k, function(k) rowMeans(CP[, jj[k], drop=FALSE])))
cc <- vapply(fv$k, function(k) mean(resp$chemo_cp10k[k]), 0)
if (ncol(cm) >= MIN_BIN) {
  rf <- suppressWarnings(apply(cm, 1, function(v) cor(v, cc, method="spearman")))
  zf <- atanh(pmax(pmin(rf,1-1e-9),-1+1e-9))*sqrt(ncol(cm)-3)
  GW[, `:=`(rho_field_fov = rf, p_fov = 2*pnorm(-abs(zf)))][, q_fov := p.adjust(p_fov,"BH")]
} else GW[, `:=`(rho_field_fov = NA_real_, p_fov = NA_real_, q_fov = NA_real_)]
GW[, robust := !is.na(q_fov) & q < Q_CUT & q_fov < Q_CUT & sign(rho_field) == sign(rho_field_fov)]
fwrite(GW, file.path(OUT,"genomewide_field.csv.gz"), compress="gzip")
cat(sprintf("\n=== genome-wide: %d of %d genes field-associated at both units (q<%.2f)\n",
            GW[robust & !is_selection, .N], nrow(GW), Q_CUT))
cat("\ntop POSITIVE (higher where the Hodgkin field is stronger):\n")
print(head(GW[robust & !is_selection & rho_field > 0][order(-rho_field),
      .(gene, rho_field=round(rho_field,3), rho_fov=round(rho_field_fov,3), q=signif(q,2))], 15))
cat("\ntop NEGATIVE:\n")
print(head(GW[robust & !is_selection & rho_field < 0][order(rho_field),
      .(gene, rho_field=round(rho_field,3), rho_fov=round(rho_field_fov,3), q=signif(q,2))], 15))
cat("\nwhere the pre-declared genes rank genome-wide:\n")
GW[, rank_pos := frank(-rho_field)]
print(GW[gene %in% PRIMARY][order(rank_pos),
      .(gene, rho_field=round(rho_field,3), rank_of = paste0(rank_pos,"/",nrow(GW)),
        q=signif(q,2), robust, is_selection)])

## ---- figure: the pre-declared axis across rings ------------------------------
gg <- intersect(setdiff(PRIMARY, SEL), GENES)
if (length(gg)) {
  pl <- rbindlist(lapply(gg, function(g) resp[, .(gene = g, ring = ring,
        cp10k = CP[g, match(bin, colnames(CP))])]))
  pm <- pl[!is.na(ring), .(mean = mean(cp10k),
           se = sd(cp10k)/sqrt(.N), n = .N), by = .(gene, ring)]
  ## Carry rho and q INTO the panel strips. With free y-axes -- which are necessary, because
  ## these genes differ 20-fold in absolute level -- a null panel spanning 1.05 to 1.09 fills
  ## the same height as a real gradient spanning 1.55 to 1.99, and reads as a trend. The
  ## statistic printed on the panel is what stops that misreading; so does the span label.
  st <- G1[set == "all non-hotspot bins", .(gene, rho = rho_field, q = q_field)]
  spans <- pm[, .(span = sprintf("%.2f-%.2f", min(mean), max(mean))), by = gene]
  st <- merge(st, spans, by = "gene")
  ## Label by EFFECT SIZE, not by a significance star. With 7,868 bins, q < 0.05 flags any
  ## |rho| above about 0.03, so a star next to PDCD1's rho of +0.031 is not information -- it
  ## invites the reader to conclude PD-1 tracks the field when the correlation is nil. The
  ## magnitude band is what separates the two real gradients from four negligible ones.
  st[, mag := fifelse(abs(rho) >= 0.30, "strong",
              fifelse(abs(rho) >= 0.10, "modest",
              fifelse(q < 0.05, "weak", "null")))]
  ## Two lines, gene name FIRST. A single-line strip carrying gene + rho + magnitude + span
  ## overflowed the panel width and ggplot clipped the START of the string, so the third
  ## column lost the gene name -- the one thing a reader needs -- while keeping the caveat.
  st[, lab := sprintf("%s   rho %+.3f\n%s   span %s", gene, rho, mag, span)]
  pm <- merge(pm, st[, .(gene, lab, rho, mag)], by = "gene")
  pm[, lab := factor(lab, levels = st[order(-rho), lab])]
  p <- ggplot(pm, aes(ring, mean, group = lab)) +
    geom_line(aes(colour = abs(rho) >= 0.10 & rho > 0), linewidth = .7) +
    geom_errorbar(aes(ymin = mean-1.96*se, ymax = mean+1.96*se), width = .12, colour = MUTED) +
    geom_point(aes(colour = abs(rho) >= 0.10 & rho > 0), size = 2) +
    scale_colour_manual(values = c(`TRUE` = DIV_HI, `FALSE` = MUTED), guide = "none") +
    facet_wrap(~ lab, scales = "free_y", ncol = 3) +
    labs(x = "distance from the nearest cHL chemokine hotspot (um)", y = "cp10k per 100 um bin",
         title = "CCR4 axis against distance from the Hodgkin chemokine field, inside IB1",
         subtitle = sprintf(paste0("tier '%s' | hotspot bins excluded | CCL17/CCL22 define the ",
           "exposure and are held out | rho is against the field, not the ring\n",
           "CCL17+CCL22 signal through CCR4. CCR4 and CTLA4 track the field; the Treg markers ",
           "FOXP3, IL2RA and IKZF2 do NOT --\nso this is a CCR4 gradient without a Treg ",
           "signature. Blue = |rho| >= 0.10 and positive. 'weak' = q < 0.05 but |rho| < 0.10,\n",
           "i.e. significant only because n is 7,868 bins. Y-axes are free: read the span, ",
           "not the slope."),
           TIER),
         caption = "Data: output/16_gradient/<tier>/ccr4_axis.csv | 95% CI") +
    theme_minimal(base_size = 11) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_line(colour=GRID, linewidth=.25),
          axis.text=element_text(colour=MUTED), axis.text.x=element_text(angle=45, hjust=1),
          axis.title=element_text(colour=INK2),
          plot.title=element_text(colour=INK, face="bold"),
          plot.subtitle=element_text(colour=INK2, size=9),
          plot.caption=element_text(colour=MUTED, size=8, hjust=0),
          strip.text=element_text(colour=INK2, face="bold", size=8, lineheight=1.15))
  ggsave(file.path(FIG,"fig_16_field_gradient.png"), p, width=10, height=7, dpi=150, bg=SURF)
  cat("wrote ", file.path(FIG,"fig_16_field_gradient.png"), "\n", sep="")
}
writeLines(format(Sys.time()), file.path(OUT,"gradient.done"))
cat("============================================================\nOK\n")
