## ============================================================================
## 16_script_matched_regions.R -- composition-matched bin pairs across the two samples.
##
## WHY THIS IS NOT OPTIONAL, established empirically by the first 16d run.
##
## C3 (IB1 FL_focus vs RB1 FL_follicular, unmatched) returned CCL19 +3.59, CCL21 +2.95,
## LYZ +3.12, APOE +3.67, PTGDS +2.67 as the top "up in IB1 FL" genes. Those are precisely
## the genes that mark RB1's INTERFOLLICULAR zone in C2 (CCL21 -2.80, CXCL12 -1.85,
## CCL19 -1.64). The contrast was not measuring FL biology at all -- it was measuring that
## IB1's FL bins are 5% B cells sitting in a myeloid/stromal node while RB1's follicular
## bins are 52% B cells.
##
## The difference-in-differences was supposed to cancel that and did not: it made it worse.
## (IB1 FL - IB1 background) - (RB1 follicular - RB1 interfollicular) returned the SAME gene
## list, because the two reference groups are not the same tissue either -- 16a showed IB1's
## B-poor tissue is myeloid-rich and T-DEPLETED (T share 0.53x, myeloid 1.86x) while RB1's is
## a real T-zone (T 2.07x, myeloid 5.09x). Subtracting a T-zone-referenced contrast from a
## myeloid-referenced one leaves the T-zone difference behind, dressed as a result.
##
## So the composition difference is removed directly, by matching, which is what step 14 did
## for depth and cell count (14:190-247) -- extended here with the composition term and
## applied ACROSS samples. BOTH the FL groups and the REFERENCE groups are matched, so the
## difference-in-differences finally has symmetric references by construction.
##
## Matched on: fb_a (B share among lineage-assigned cells -- the composition axis),
## log cell count, log median depth, and f_assigned (how typable the bin is). Greedy
## nearest-neighbour on standardised covariates, without replacement, with a caliper.
##
## Usage: Rscript scripts/16_script_matched_regions.R [tier]
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2); setDTthreads(1) })
if (!file.exists("scripts/lib/lineage.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "strict"

BIN_F <- file.path("output/16_regions", TIER, "bins_all_scales.csv.gz")
OUT   <- file.path("output/16_matched", TIER)
FIG   <- file.path(OUT, "figures")
if (!file.exists(BIN_F)) stop("missing ", BIN_F, call. = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

FOCUS_BIN <- 100L
COVS      <- c("fb_a","log_cells","log_depth","f_assigned")
CALIPER   <- 0.5     # max standardised distance; unmatched IB1 bins are REPORTED, not forced
SEED      <- 1L
SURF<-"#fcfcfb"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
DIV_HI<-"#2a78d6"; DIV_LO<-"#eb6834"

cat("============================================================\n")
cat("step 16c -- cross-sample composition-matched bins, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

B <- fread(BIN_F)[bin_um == FOCUS_BIN]
B[, bin := paste(samp, bx, by_, sep = "_")]
B[, `:=`(log_cells = log(cells), log_depth = log(med_depth))]
B <- B[!is.na(fb_a) & is.finite(log_depth)]
cat(sprintf("bins at %d um: %d\n", FOCUS_BIN, nrow(B)))
print(B[, .N, by = .(samp, bin_region)][order(samp,-N)])

## Standardise on the POOLED spread of the two groups being matched, recomputed per pair --
## using a global SD would make the caliper mean different things for different contrasts.
match_groups <- function(lab, a_sel, b_sel) {
  A <- B[a_sel]; Bp <- B[b_sel]
  if (!nrow(A) || !nrow(Bp)) { cat("\n", lab, ": one side empty; skipped\n", sep=""); return(NULL) }
  mu <- colMeans(rbind(as.matrix(A[, ..COVS]), as.matrix(Bp[, ..COVS])))
  sd <- apply(rbind(as.matrix(A[, ..COVS]), as.matrix(Bp[, ..COVS])), 2, stats::sd)
  sd[sd == 0] <- 1
  za <- scale(as.matrix(A[, ..COVS]), center = mu, scale = sd)
  zb <- scale(as.matrix(Bp[, ..COVS]), center = mu, scale = sd)
  ## Greedy nearest neighbour without replacement. IB1 bins are the scarce side (44 FL bins
  ## against 2605 RB1 follicular), so they drive the matching and each takes its own best
  ## available partner. Deterministic given the row order, which is fixed by the seed.
  set.seed(SEED)
  used <- rep(FALSE, nrow(zb)); pick <- rep(NA_integer_, nrow(za)); dist <- rep(NA_real_, nrow(za))
  ord <- sample.int(nrow(za))
  for (i in ord) {
    d <- sqrt(colSums((t(zb) - za[i, ])^2))
    d[used] <- Inf
    j <- which.min(d)
    if (is.finite(d[j]) && d[j] <= CALIPER) { pick[i] <- j; dist[i] <- d[j]; used[j] <- TRUE }
  }
  ok <- !is.na(pick)
  cat(sprintf("\n=== %s\n  %d vs %d bins | matched %d of %d (%.0f%%) within caliper %.2f\n",
              lab, nrow(A), nrow(Bp), sum(ok), nrow(A), 100*mean(ok), CALIPER))
  if (!sum(ok)) { cat("  NOTHING matched -- the two groups do not overlap on these covariates\n")
                  return(NULL) }
  ## Standardised mean difference per covariate, before and after. This is the whole point of
  ## the script: a matching step that silently failed looks identical to one that worked.
  bal <- rbindlist(lapply(COVS, function(v) data.table(
    group = lab, covariate = v,
    smd_before = (mean(A[[v]]) - mean(Bp[[v]])) / sd[[v]],
    smd_after  = (mean(A[[v]][ok]) - mean(Bp[[v]][pick[ok]])) / sd[[v]],
    a_before = mean(A[[v]]), b_before = mean(Bp[[v]]),
    a_after  = mean(A[[v]][ok]), b_after = mean(Bp[[v]][pick[ok]]))))
  bal[, `:=`(smd_before = round(smd_before,3), smd_after = round(smd_after,3),
             a_before = round(a_before,3), b_before = round(b_before,3),
             a_after = round(a_after,3), b_after = round(b_after,3))]
  print(bal[, .(covariate, a_before, b_before, smd_before, a_after, b_after, smd_after)])
  worst <- max(abs(bal$smd_after))
  cat(sprintf("  worst |SMD| after matching: %.3f -> %s\n", worst,
              if (worst < 0.1) "BALANCED (<0.1)" else if (worst < 0.25) "acceptable (<0.25)"
              else "STILL IMBALANCED -- do not read this contrast as composition-free"))
  list(pairs = data.table(group = lab, a_bin = A$bin[ok], b_bin = Bp$bin[pick[ok]],
                          a_samp = A$samp[ok], b_samp = Bp$samp[pick[ok]],
                          a_fov = A$fov[ok], b_fov = Bp$fov[pick[ok]],
                          dist = round(dist[ok],4)),
       balance = bal, unmatched = nrow(A) - sum(ok), worst_smd = worst)
}

IB <- B$samp == "IB1"; RB <- B$samp == "RB1"; R <- B$bin_region
res <- list(
  ## the FL-to-FL contrast: the question
  fl  = match_groups("FL: IB1 FL_focus ~ RB1 FL_follicular",
                     IB & R == "FL_focus",   RB & R == "FL_follicular"),
  ## the reference-to-reference contrast, so a difference-in-differences has symmetric
  ## references instead of comparing a myeloid-rich background to a T-zone
  ref = match_groups("REF: IB1 background ~ RB1 (any non-follicular)",
                     IB & R == "background", RB & R %in% c("interfollicular","follicular_mid")),
  ## and the cHL hotspots against composition-matched RB1 tissue, which is the cleanest
  ## available test of whether a hotspot is more than its cell mix
  hot = match_groups("HOT: IB1 cHL_hotspot ~ RB1 (any)",
                     IB & R == "cHL_hotspot", RB & R != "unbinned"))
res <- res[!vapply(res, is.null, TRUE)]
if (!length(res)) stop("no group matched at all", call. = FALSE)

PAIRS <- rbindlist(lapply(res, `[[`, "pairs"))
BAL   <- rbindlist(lapply(res, `[[`, "balance"))
sc("matched_pairs.csv", PAIRS)
sc("matched_balance.csv", BAL)
sm <- rbindlist(lapply(names(res), function(k) data.table(
  group = res[[k]]$pairs$group[1], pairs = nrow(res[[k]]$pairs),
  unmatched_a = res[[k]]$unmatched, worst_smd_after = round(res[[k]]$worst_smd,3),
  a_fovs = uniqueN(res[[k]]$pairs$a_fov), b_fovs = uniqueN(res[[k]]$pairs$b_fov))))
cat("\n=== summary\n"); print(sm); sc("matched_summary.csv", sm)

## Before/after balance figure -- one row per covariate, the SMD collapsing toward zero.
pl <- melt(BAL[, .(group, covariate, before = smd_before, after = smd_after)],
           id.vars = c("group","covariate"), variable.name = "stage", value.name = "smd")
p <- ggplot(pl, aes(abs(smd), covariate, colour = stage)) +
  geom_vline(xintercept = c(0.1, 0.25), linetype = "22", colour = MUTED, linewidth = .3) +
  geom_point(size = 2.6, alpha = .9) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(before = DIV_LO, after = DIV_HI), name = NULL) +
  labs(x = "|standardised mean difference|", y = NULL,
       title = "Cross-sample composition matching",
       subtitle = sprintf("tier '%s' | dashed lines at 0.1 (balanced) and 0.25 (acceptable)", TIER),
       caption = "Data: output/16_matched/<tier>/matched_balance.csv") +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = GRID, linewidth = .25),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK2),
        strip.text = element_text(colour = INK2, face = "bold", size = 9),
        legend.position = "top", legend.justification = "left")
ggsave(file.path(FIG, "fig_16_matching_balance.png"), p, width = 8, height = 7, dpi = 150, bg = SURF)
cat("wrote ", file.path(FIG, "fig_16_matching_balance.png"), "\n", sep="")
writeLines(format(Sys.time()), file.path(OUT, "matched.done"))
cat("============================================================\nOK\n")
