## ============================================================================
## 16_script_composition.R -- composition, and a difference-in-differences ON composition.
##
## THE IDEA. The cross-sample EXPRESSION contrast is not identifiable (16c: 2 of 44 bins
## match). But composition ENRICHMENT is a dimensionless ratio -- how much more B-rich is FL
## territory than its own surroundings -- and a ratio does not care that IB1 sits at 17.6% B
## and RB1 at 50.9%. So the question "does FL organise its microenvironment differently in
## the two nodes" is still askable, on exactly the axis that broke the expression contrast.
##
## The marginal percentages hint at something sharp: B enrichment of FL over its own
## reference is 2.15x in IB1 and 2.13x in RB1 -- the same -- while T goes 1.65x in IB1 and
## 0.54x in RB1, and myeloid 0.81x vs 0.25x. That would say the FL concentrates B cells
## identically in both nodes but retains T cells and macrophages in the composite node that
## it excludes in the recurrence. This script tests it properly, because those marginals are
## pooled over cells and therefore pseudo-replicated.
##
## THE UNIT IS THE BIN, and the enrichment is computed PER BIN PAIRING via a bootstrap over
## bins, then over FOVs as the conservative reading. The reference-group asymmetry that sank
## the expression difference-in-differences still applies and is NOT solved here -- it is
## reported: IB1's reference is myeloid-rich and T-depleted, RB1's is a T-zone. That matters
## much less for the T axis (IB1 background 7.9% T vs RB1 interfollicular 9.1% -- close) than
## for the B axis (8.2% vs 23.9% -- far), so the T and myeloid results are the defensible
## ones and the B result is reported with that caveat attached.
##
## Usage: Rscript scripts/16_script_composition.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); setDTthreads(1)
})
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "strict"

REG_F <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT_F  <- "output/16_celltypes/cell_celltypes.csv.gz"
OUT   <- file.path("output/16_composition", TIER)
FIG   <- file.path(OUT, "figures")
for (f in c(REG_F, CT_F)) if (!file.exists(f)) stop("missing ", f, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

N_BOOT <- 2000L; SEED <- 1L; MIN_BINS <- 12L
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
DIV_HI<-"#2a78d6"; DIV_LO<-"#eb6834"

cat("============================================================\n")
cat("step 16k -- composition and composition DiD, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

rg <- fread(REG_F, select = c("cell_id","samp","bx","by_","bin_region","arch","is_defect",
                              "unbinned","fov_id"))
ct <- fread(CT_F, select = c("cell_id","celltype","lineage16"))
d  <- merge(rg, ct, by = "cell_id", sort = FALSE)[is_defect == FALSE & unbinned == FALSE]
d[, bin := paste(samp, bx, by_, sep = "_")]
cat(sprintf("cells: %d\n", nrow(d)))

## Per-bin fractions, coarse and fine. Fractions are of ASSIGNED cells, for the same reason
## 16a switched the architecture axis to that denominator: otherwise "fewer T cells" and
## "fewer typable cells" are the same number.
SENT <- c("ambiguous","below_background","all_markers_0")
d[, assigned := !lineage16 %in% SENT]
BIN <- d[, .(cells = .N, n_assigned = sum(assigned), fov = fov_id[1],
             samp = samp[1], bin_region = bin_region[1], arch = arch[1]), by = bin]
add_frac <- function(col, types) {
  for (k in types) {
    tmp <- d[get(col) == k, .(n = .N), by = bin]
    BIN[tmp, on = "bin", (paste0("n_", k)) := i.n]
    BIN[is.na(get(paste0("n_", k))), (paste0("n_", k)) := 0L]
    BIN[, (paste0("f_", k)) := get(paste0("n_", k)) / pmax(n_assigned, 1)]
  }
}
COARSE <- c("B","T","myeloid","PC","NK","endo","fibro")
FINE   <- setdiff(unique(d$celltype), c(SENT, COARSE))
add_frac("lineage16", COARSE); add_frac("celltype", FINE)
TYPES <- c(COARSE, FINE)
BIN <- BIN[n_assigned >= 10]
cat(sprintf("bins with >= 10 assigned cells: %d\n", nrow(BIN)))
print(BIN[, .(bins = .N), by = .(samp, bin_region)][order(samp,-bins)])

## ---- 1. composition per region, bin as the unit ------------------------------
comp <- rbindlist(lapply(TYPES, function(k) BIN[, .(
  celltype = k, bins = .N,
  mean_frac = mean(get(paste0("f_", k))), med_frac = median(get(paste0("f_", k))),
  cells = sum(get(paste0("n_", k)))), by = .(samp, bin_region)]))
comp[, `:=`(mean_frac = round(mean_frac,4), med_frac = round(med_frac,4))]
sc("composition_by_region.csv", comp)
cat("\n=== fraction of ASSIGNED cells, mean over bins\n")
print(dcast(comp[celltype %in% c("B","T","myeloid","PC")],
            samp + bin_region ~ celltype, value.var = "mean_frac"))

## ---- 2. enrichment of FL territory over its own reference --------------------
## FL vs reference, within each sample, then the ratio of those ratios across samples.
## Bootstrap resamples BINS; the FOV version resamples FOVS and is the conservative read.
## TWO DEFINITIONS OF "FL TERRITORY", run side by side, because the first run made it
## clear that one of them cannot be trusted on its own.
##
## Under the `lightchain` pairing every single cell type came out enriched in IB1's FL
## relative to RB1's -- B by 0.39 log2 but T by 1.73 and myeloid by 1.83. A uniform shift
## across every lineage is the signature of a DEFINITIONAL difference, not of biology:
## IB1's FL_focus is a 40 um halo drawn around scattered light-chain-positive cells, so it
## necessarily sweeps in whatever those cells are sitting next to, whereas RB1's
## FL_follicular is a top-quartile B-density bin and is B-dominated by construction.
##
## The `architecture` pairing removes that asymmetry: BOTH samples use B_dense vs B_sparse,
## the identical within-sample quantile rule. If the T and myeloid results survive it, they
## are about the tissue. If they collapse, they were about the region definition.
PAIRING <- list(
  lightchain   = list(col = "bin_region",
                      IB1 = list(fl = "FL_focus",      ref = "background"),
                      RB1 = list(fl = "FL_follicular", ref = "interfollicular")),
  architecture = list(col = "arch",
                      IB1 = list(fl = "B_dense", ref = "B_sparse"),
                      RB1 = list(fl = "B_dense", ref = "B_sparse")))
enrich_boot <- function(k, S, by_fov, PN) {
  P <- PAIRING[[PN]]; cc <- P$col
  fl <- BIN[samp == S & get(cc) == P[[S]]$fl]
  rf <- BIN[samp == S & get(cc) == P[[S]]$ref]
  if (nrow(fl) < MIN_BINS || nrow(rf) < MIN_BINS) return(NULL)
  nc <- paste0("n_", k)
  ## Enrichment = (cells of type k / assigned cells) in FL, over the same in the reference.
  ## Pooled within each resample, so a bin contributes in proportion to its size.
  est <- function(a, b) (sum(fl[[nc]][a])/max(sum(fl$n_assigned[a]),1)) /
                        pmax((sum(rf[[nc]][b])/max(sum(rf$n_assigned[b]),1)), 1e-9)
  obs <- est(seq_len(nrow(fl)), seq_len(nrow(rf)))
  set.seed(SEED)
  if (by_fov) {
    fa <- split(seq_len(nrow(fl)), fl$fov); fb <- split(seq_len(nrow(rf)), rf$fov)
    if (length(fa) < MIN_BINS || length(fb) < MIN_BINS) return(NULL)
    bs <- vapply(seq_len(N_BOOT), function(i)
      est(unlist(fa[sample.int(length(fa), length(fa), TRUE)]),
          unlist(fb[sample.int(length(fb), length(fb), TRUE)])), 0)
  } else {
    bs <- vapply(seq_len(N_BOOT), function(i)
      est(sample.int(nrow(fl), nrow(fl), TRUE), sample.int(nrow(rf), nrow(rf), TRUE)), 0)
  }
  data.table(pairing = PN, samp = S, celltype = k, unit = if (by_fov) "FOV" else "bin",
             enrich = obs, lo = quantile(bs,.025), hi = quantile(bs,.975),
             n_fl = nrow(fl), n_ref = nrow(rf))
}
EN <- rbindlist(lapply(names(PAIRING), function(PN)
  rbindlist(lapply(c(FALSE, TRUE), function(bf)
    rbindlist(lapply(c("IB1","RB1"), function(S)
      rbindlist(lapply(TYPES, enrich_boot, S = S, by_fov = bf, PN = PN))))))), fill = TRUE)
EN[, `:=`(enrich = round(enrich,3), lo = round(lo,3), hi = round(hi,3))]
sc("fl_enrichment.csv", EN)
cat("\n=== enrichment of FL territory over its OWN reference (>1 = enriched in FL)\n")
print(dcast(EN[celltype %in% c("B","T","myeloid","PC","NK","fibro") & unit == "bin"],
            pairing + celltype ~ samp, value.var = c("enrich","lo","hi"))[order(pairing, celltype)])

## ---- 3. the composition DiD --------------------------------------------------
## ratio of the two enrichments. CI by pairing the two samples' bootstrap draws.
did_boot <- function(k, by_fov, PN) {
  P <- PAIRING[[PN]]; cc <- P$col
  g <- function(S) {
    fl <- BIN[samp == S & get(cc) == P[[S]]$fl]
    rf <- BIN[samp == S & get(cc) == P[[S]]$ref]
    nc <- paste0("n_", k)
    est <- function(a,b) (sum(fl[[nc]][a])/max(sum(fl$n_assigned[a]),1)) /
                         pmax((sum(rf[[nc]][b])/max(sum(rf$n_assigned[b]),1)), 1e-9)
    if (nrow(fl) < MIN_BINS || nrow(rf) < MIN_BINS) return(NULL)
    set.seed(SEED)
    if (by_fov) {
      fa <- split(seq_len(nrow(fl)), fl$fov); fb <- split(seq_len(nrow(rf)), rf$fov)
      if (length(fa) < MIN_BINS || length(fb) < MIN_BINS) return(NULL)
      list(obs = est(seq_len(nrow(fl)), seq_len(nrow(rf))),
           bs = vapply(seq_len(N_BOOT), function(i)
             est(unlist(fa[sample.int(length(fa), length(fa), TRUE)]),
                 unlist(fb[sample.int(length(fb), length(fb), TRUE)])), 0))
    } else {
      list(obs = est(seq_len(nrow(fl)), seq_len(nrow(rf))),
           bs = vapply(seq_len(N_BOOT), function(i)
             est(sample.int(nrow(fl), nrow(fl), TRUE), sample.int(nrow(rf), nrow(rf), TRUE)), 0))
    }
  }
  a <- g("IB1"); b <- g("RB1"); if (is.null(a) || is.null(b)) return(NULL)
  r  <- log2(a$obs / b$obs); rb <- log2(a$bs / b$bs)
  data.table(pairing = PN, celltype = k, unit = if (by_fov) "FOV" else "bin",
             enr_IB1 = a$obs, enr_RB1 = b$obs, log2_did = r,
             lo = quantile(rb,.025), hi = quantile(rb,.975),
             p_boot = pmax(2*min(mean(rb <= 0), mean(rb >= 0)), 1/N_BOOT))
}
DD <- rbindlist(lapply(names(PAIRING), function(PN)
  rbindlist(lapply(c(FALSE,TRUE), function(bf)
    rbindlist(lapply(TYPES, did_boot, by_fov = bf, PN = PN))))), fill = TRUE)
if (nrow(DD)) {
  DD[, q_boot := p.adjust(p_boot, "BH"), by = .(pairing, unit)]
  DD[, excludes_zero := (lo > 0) | (hi < 0)]
  DD[, `:=`(enr_IB1 = round(enr_IB1,3), enr_RB1 = round(enr_RB1,3),
            log2_did = round(log2_did,3), lo = round(lo,3), hi = round(hi,3))]
  sc("composition_did.csv", DD)
  cat("\n=== composition DiD: log2( IB1 FL-enrichment / RB1 FL-enrichment )\n")
  cat("    0 = the FL concentrates this cell type equally in both nodes\n")
  print(DD[order(pairing, unit, -abs(log2_did)),
           .(pairing, unit, celltype, enr_IB1, enr_RB1, log2_did,
             ci = sprintf("[%.2f, %.2f]", lo, hi), q = signif(q_boot,2), excludes_zero)])
  ## The comparison that decides whether any of this is biology.
  ## ---- NORMALISE TO THE B-CELL ENRICHMENT ------------------------------------
  ## The raw DiD came out POSITIVE for all 16 cell types under both definitions, which is
  ## not 16 findings -- it is one scale artifact. A within-sample quantile contrast has very
  ## different dynamic range in the two samples: B_dense vs B_sparse spans a B share of
  ## 0.286 -> 0.029 in IB1 (a ~10x swing) but only 0.790 -> 0.484 in RB1 (~1.6x), because
  ## RB1 is uniformly B-rich. So IB1's B enrichment is 17.99x against RB1's 1.65x before any
  ## other cell type is considered, and everything inherits that.
  ##
  ## Dividing each cell type's enrichment by the B enrichment IN THE SAME SAMPLE removes it.
  ## The quantity that survives is interpretable on its own terms: how strongly does this
  ## population follow the B cells into FL territory, as a fraction of how strongly B does.
  ## A ratio near 1 between samples means the FL recruits that population to the same
  ## relative degree in both nodes.
  ENB <- EN[unit == "bin"]
  bref <- ENB[celltype == "B", .(pairing, samp, b_enrich = enrich)]
  ENB <- merge(ENB, bref, by = c("pairing","samp"))
  ENB[, rel_to_B := enrich / b_enrich]
  relw <- dcast(ENB, pairing + celltype ~ samp, value.var = c("enrich","rel_to_B"))
  relw[, follows_B_ratio := rel_to_B_IB1 / rel_to_B_RB1]
  relw <- relw[order(pairing, -rel_to_B_IB1)]
  sc("enrichment_relative_to_B.csv", relw)
  cat("\n=== enrichment RELATIVE TO B in the same sample (symmetric `architecture` pairing)\n")
  cat("    rel_to_B = how far this type follows B into FL territory, as a fraction of B.\n")
  cat("    follows_B_ratio near 1 = the FL recruits it to the same relative degree in both nodes.\n")
  print(relw[pairing == "architecture",
             .(celltype, enr_IB1 = round(enrich_IB1,3), enr_RB1 = round(enrich_RB1,3),
               rel_IB1 = round(rel_to_B_IB1,3), rel_RB1 = round(rel_to_B_RB1,3),
               follows_B_ratio = round(follows_B_ratio,3))])
  top <- relw[pairing == "architecture" & celltype != "B"][order(-rel_to_B_IB1)][1]
  cat(sprintf("\nmost FL-concentrated non-B population: %s (IB1 %.3f, RB1 %.3f of the B enrichment)\n",
              top$celltype, top$rel_to_B_IB1, top$rel_to_B_RB1))

  cmp <- dcast(DD[unit == "bin"], celltype ~ pairing, value.var = c("log2_did","excludes_zero"))
  if (all(c("log2_did_lightchain","log2_did_architecture") %in% names(cmp))) {
    cmp[, survives_both := excludes_zero_lightchain & excludes_zero_architecture &
                           sign(log2_did_lightchain) == sign(log2_did_architecture)]
    cat("\n=== does the DiD survive BOTH region definitions?\n")
    print(cmp[order(-abs(log2_did_architecture)),
              .(celltype, lightchain = round(log2_did_lightchain,2),
                architecture = round(log2_did_architecture,2), survives_both)])
    sc("composition_did_definition_check.csv", cmp)
    cat(sprintf("\nsurviving both definitions: %d of %d cell types\n",
                cmp[survives_both == TRUE, .N], nrow(cmp)))
  }
  cat("\nNOTE the reference groups differ (16a): IB1's is myeloid-rich and T-depleted,\n")
  cat("RB1's is a genuine T-zone. IB1 background is 7.9% T vs RB1 interfollicular 9.1%\n")
  cat("(close), but 8.2% vs 23.9% B (far). Read the T and myeloid rows as evidence and\n")
  cat("the B row as descriptive only.\n")

  pl <- melt(relw[pairing == "architecture" & celltype %in% c("Tfh","T","CD4_conv","CD8_T",
              "Treg","Macrophage","Monocyte","myeloid","fibro","NK","PC","endo")],
             id.vars = "celltype", measure.vars = c("rel_to_B_IB1","rel_to_B_RB1"),
             variable.name = "samp", value.name = "rel_to_B")
  pl[, samp := factor(sub("rel_to_B_", "", samp), levels = c("IB1","RB1"))]
  if (nrow(pl)) {
    p <- ggplot(pl, aes(rel_to_B, reorder(celltype, rel_to_B))) +
      geom_line(aes(group = celltype), colour = GRID, linewidth = 1.2) +
      geom_point(aes(colour = samp), size = 3) +
      scale_colour_manual(values = c(IB1 = DIV_HI, RB1 = DIV_LO), name = NULL) +
      labs(x = "enrichment in FL territory, as a fraction of the B-cell enrichment", y = NULL,
           title = "Does FL organise its microenvironment differently in the two nodes?",
           subtitle = paste0("Symmetric definition: B_dense vs B_sparse bins, identical rule in both samples.\n",
                             "Normalised to each sample's own B enrichment, which removes the dynamic-range\n",
                             "difference (B enriches 17.99x in IB1 but only 1.65x in RB1)."),
           caption = "Data: output/16_composition/<tier>/enrichment_relative_to_B.csv") +
      theme_minimal(base_size = 11) +
      theme(plot.background = element_rect(fill=SURF,colour=NA),
            panel.background = element_rect(fill=SURF,colour=NA),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour=GRID, linewidth=.25),
            panel.grid.minor = element_blank(),
            axis.text = element_text(colour=MUTED), axis.title = element_text(colour=INK2),
            plot.title = element_text(colour=INK, face="bold"),
            plot.subtitle = element_text(colour=INK2, size=9),
            plot.caption = element_text(colour=MUTED, size=8, hjust=0))
    ggsave(file.path(FIG,"fig_16_composition_did.png"), p, width=9, height=6, dpi=150, bg=SURF)
    cat("wrote ", file.path(FIG,"fig_16_composition_did.png"), "\n", sep="")
  }
}
writeLines(format(Sys.time()), file.path(OUT, "composition.done"))
cat("============================================================\nOK\n")
