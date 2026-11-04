## ============================================================================
## Step 32 -- Do candidate HRS cells sit inside a WIDER CD45-poor zone? Radial CD45 to
## 500 um, with the over-fixation defect burden as the covariate.
##
## WHY. Step 31 (5 um rings to 100 um, each profile over its own 50-100 um mean) found
## candidates CD45-dim at the centre and no rosette ring, but its normalisation cancels
## any CD45-poor zone wider than ~100 um by construction. Two hints that one exists: a
## shallow -0.02 to -0.03 dip from 15-60 um in candidates only, and raw 50-100 um CD45
## ~6% dimmer around candidates than around same-FOV controls (log2 -0.09 [-0.26, 0.01]).
## The chemokine field itself runs 400-800 um (CLAUDE.md fact 3), and the defect -- whose
## cells are dropped from the tables but still sit in the image, dim -- is confounded with
## it (fact 9: hotspots 0.59 defect fraction vs 0.45 background).
##
## DESIGN (fixed before running)
##   cells      step 31's cells.csv: 500 random >= 50 candidates, 500 area-matched < 10
##              controls (sample-wide draw, so controls sample ordinary tissue)
##   image      pyramid level 3 (~0.96 um/px), mean CD45 and DNA in 25 um bands 0-500 um
##              (32_wide_cd45_profiles.py). A band with > 5% blank pixels is NA; a cell
##              whose 400-500 um reference ring is > 50% blank is dropped.
##   reference  each profile over its OWN 400-500 um mean, log2. Slide-level staining
##              variation cancels; a zone up to ~400 um across does not.
##   defect     per cell and band, fraction of ALL symmetric-tier cells (defect included)
##              within the band that are InSituType `i`, from cell_regions in slide mm.
##   contrast   per band: candidate minus control mean log2 (FOV cluster bootstrap), and
##              the group coefficient of lm(log2 rel ~ group + band defect fraction), same
##              bootstrap. If the dip survives adjustment it is not (only) the defect.
##
##   Rscript scripts/32_script_wide_cd45.R select|analyze [smoke]
## ============================================================================
suppressPackageStartupMessages(library(data.table))
args  <- commandArgs(trailingOnly = TRUE)
MODE  <- args[1]
SMOKE <- "smoke" %in% args
S31   <- if (SMOKE) "output/31_radial_cd45_smoke" else "output/31_radial_cd45"
OUT   <- if (SMOKE) "output/32_wide_cd45_smoke" else "output/32_wide_cd45"
FIG   <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
BAND_UM <- 25; R_MAX <- 500; NB <- R_MAX / BAND_UM
REF <- which(seq_len(NB) * BAND_UM > 400)              # 400-500 um
band_lab <- sprintf("%d-%d", (seq_len(NB) - 1L) * BAND_UM, seq_len(NB) * BAND_UM)
band_mid <- (seq_len(NB) - 0.5) * BAND_UM
set.seed(32)

if (MODE == "select") {
  ## Reuse step 31's cells; add slide-mm and the per-band defect fraction here, in R,
  ## so the Python stage only touches the image.
  suppressPackageStartupMessages(library(RANN))
  cells <- fread(file.path(S31, "cells.csv"))
  reg <- fread(cmd = "zcat output/16_regions/matched/cell_regions.csv.gz",
               select = c("cell_id", "samp", "x_slide_mm", "y_slide_mm", "is_defect"))[samp == "IB1"]
  ## step-16 argmax lineage, for a composition profile: a CD45-poor zone could simply be
  ## B-dense tissue (B cells are the dimmest CD45 lineage here, step 30 gate table)
  l16 <- fread(cmd = "zcat output/16_celltypes/cell_celltypes.csv.gz",
               select = c("cell_id", "samp", "lineage16"))[samp == "IB1"]
  reg <- merge(reg, l16[, .(cell_id, lineage16)], by = "cell_id", all.x = TRUE)
  LIN <- c("B", "T", "myeloid", "fibro", "endo", "PC", "NK")
  reg[!lineage16 %in% LIN | is_defect == TRUE, lineage16 := NA_character_]
  cells <- merge(cells, reg[, .(cell_id, x_slide_mm, y_slide_mm)], by = "cell_id")
  if (!SMOKE) stopifnot(nrow(cells) == 1000L)
  allxy <- as.matrix(reg[, .(x_slide_mm * 1000, y_slide_mm * 1000)])
  q <- as.matrix(cells[, .(x_slide_mm * 1000, y_slide_mm * 1000)])
  K <- 15000L
  nn <- nn2(allxy, q, k = K, searchtype = "radius", radius = R_MAX)
  stopifnot(all(nn$nn.idx[, K] == 0L))                  # K was large enough everywhere
  isdef <- reg$is_defect; lin <- reg$lineage16
  prof <- t(sapply(seq_len(nrow(q)), function(i) {
    ok <- nn$nn.idx[i, ] > 0L
    j <- nn$nn.idx[i, ok]
    b <- pmin(as.integer(nn$nn.dists[i, ok] / BAND_UM), NB - 1L) + 1L
    n <- tabulate(b, NB); k <- tabulate(b[isdef[j]], NB)
    lab <- !is.na(lin[j])
    c(n, k, tabulate(b[lab], NB),
      unlist(lapply(c("B", "T", "myeloid"), function(L) tabulate(b[lab & lin[j] == L], NB))))
  }))
  colnames(prof) <- c(sprintf("ncell_b%02d", seq_len(NB) - 1L), sprintf("ndef_b%02d", seq_len(NB) - 1L),
                      sprintf("nlab_b%02d", seq_len(NB) - 1L),
                      unlist(lapply(c("B", "T", "myeloid"), function(L) sprintf("n%s_b%02d", L, seq_len(NB) - 1L))))
  cells <- cbind(cells, as.data.table(prof))
  fwrite(cells, file.path(OUT, "cells.csv"))
  cat(sprintf("wrote %d cells with per-band defect counts (K = %d, max neighbours used %d)\n",
              nrow(cells), K, max(rowSums(nn$nn.idx > 0L))))
  quit(save = "no")
}

stopifnot(MODE == "analyze")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
N_BOOT <- if (SMOKE) 20L else 1000L
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
theme_viz <- function(base=11, xgrid=FALSE) {
  th <- theme_minimal(base_size=base) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.major.y=element_line(colour=GRID,linewidth=.25),
          panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.line.x=element_line(colour=AXIS,linewidth=.3), axis.ticks=element_blank(),
          axis.text=element_text(colour=MUTED,size=base-2),
          axis.title=element_text(colour=INK2,size=base-1),
          plot.title=element_text(colour=INK,face="bold",size=base+2),
          plot.subtitle=element_text(colour=INK2,size=base-1),
          plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
          legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2),
          legend.position="top", legend.justification="left", legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF);cat("wrote",file.path(FIG,n),"\n")}
GCOL <- c(candidate = "#0b3d7a", control = "#9a968d")

cells <- fread(file.path(OUT, "cells.csv"))
pr <- merge(cells, fread(file.path(OUT, "profiles.csv.gz")), by = c("cell_id", "group"))
cat(sprintf("profiled: %s; dropped: %s\n",
            paste(pr[ok == TRUE, .N, by = group][, sprintf("%s %d", group, N)], collapse = ", "),
            paste(pr[ok == FALSE, .N, by = .(group, reason)][, sprintf("%s %s %d", group, reason, N)],
                  collapse = ", ")))
pr <- pr[ok == TRUE]
CD <- sprintf("cd45_b%02d", seq_len(NB) - 1L); DN <- sprintf("dna_b%02d", seq_len(NB) - 1L)
BL <- sprintf("blank_b%02d", seq_len(NB) - 1L)
NC <- sprintf("ncell_b%02d", seq_len(NB) - 1L); ND <- sprintf("ndef_b%02d", seq_len(NB) - 1L)
raw <- as.matrix(pr[, ..CD]); dna <- as.matrix(pr[, ..DN]); blank <- as.matrix(pr[, ..BL])
raw[blank > 0.05] <- NA; dna[blank > 0.05] <- NA
ref <- rowMeans(raw[, REF, drop = FALSE], na.rm = TRUE)
## a cell with every reference band NA (ring partly blank) has no reference: drop it
keep <- is.finite(ref) & ref > 0
cat(sprintf("cells with a usable 400-500 um reference: %d of %d\n", sum(keep), length(keep)))
pr <- pr[keep]; raw <- raw[keep, ]; dna <- dna[keep, ]; blank <- blank[keep, ]; ref <- ref[keep]
rel <- log2(raw / ref); dna_rel <- log2(dna / rowMeans(dna[, REF, drop = FALSE], na.rm = TRUE))
NL <- sprintf("nlab_b%02d", seq_len(NB) - 1L)
nlab <- as.matrix(pr[, ..NL])
share <- lapply(c("B", "T", "myeloid"), function(L) {
  M <- as.matrix(pr[, sprintf("n%s_b%02d", L, seq_len(NB) - 1L), with = FALSE]) / pmax(nlab, 1)
  M[nlab < 5] <- NA; colnames(M) <- band_lab; M })
names(share) <- c("B", "T", "myeloid")
defrac <- as.matrix(pr[, ..ND]) / pmax(as.matrix(pr[, ..NC]), 1)
defrac[as.matrix(pr[, ..NC]) < 5] <- NA
colnames(rel) <- colnames(dna_rel) <- colnames(defrac) <- band_lab
fwrite(cbind(pr[, .(cell_id, group, fov_id, Area.um2, chemo, hrs_score, cd45_ref = round(ref, 1))],
             as.data.table(round(rel, 4)),
             setnames(as.data.table(round(defrac, 4)), paste0("defect_", band_lab))),
       file.path(OUT, "profiles_relative.csv.gz"))

## ---- 1. average profiles and the group contrast, raw and defect-adjusted ----------
gi <- list(candidate = which(pr$group == "candidate"), control = which(pr$group == "control"))
fovs <- unique(pr$fov_id); idx <- split(seq_len(nrow(pr)), pr$fov_id)
stat_of <- function(rows) {                 # rows may repeat (bootstrap draw)
  g <- pr$group[rows]
  m_c <- colMeans(rel[rows[g == "candidate"], , drop = FALSE], na.rm = TRUE)
  m_k <- colMeans(rel[rows[g == "control"], , drop = FALSE], na.rm = TRUE)
  d_c <- colMeans(defrac[rows[g == "candidate"], , drop = FALSE], na.rm = TRUE)
  d_k <- colMeans(defrac[rows[g == "control"], , drop = FALSE], na.rm = TRUE)
  adj <- sapply(seq_len(NB), function(b) {
    y <- rel[rows, b]; x <- defrac[rows, b]; gg <- as.integer(g == "candidate")
    ok <- is.finite(y) & is.finite(x)
    if (sum(ok) < 10 || sd(x[ok]) == 0) return(c(NA, NA))
    coef(lm(y[ok] ~ gg[ok] + x[ok]))[2:3]
  })
  list(m_c = m_c, m_k = m_k, diff = m_c - m_k, d_c = d_c, d_k = d_k,
       adj = adj[1, ], slope_def = adj[2, ])
}
S0 <- stat_of(seq_len(nrow(pr)))
B <- lapply(seq_len(N_BOOT), function(b)
  stat_of(unlist(idx[as.character(sample(fovs, replace = TRUE))], use.names = FALSE)))
ci <- function(what) {
  M <- do.call(rbind, lapply(B, `[[`, what))
  list(lo = apply(M, 2, quantile, .025, na.rm = TRUE), hi = apply(M, 2, quantile, .975, na.rm = TRUE))
}
res <- data.table(band = band_lab, r_um = band_mid,
                  cand = S0$m_c, cand_lo = ci("m_c")$lo, cand_hi = ci("m_c")$hi,
                  ctrl = S0$m_k, ctrl_lo = ci("m_k")$lo, ctrl_hi = ci("m_k")$hi,
                  diff = S0$diff, diff_lo = ci("diff")$lo, diff_hi = ci("diff")$hi,
                  defect_cand = S0$d_c, defect_ctrl = S0$d_k,
                  diff_adj = S0$adj, adj_lo = ci("adj")$lo, adj_hi = ci("adj")$hi,
                  slope_defect = S0$slope_def,
                  dna_cand = colMeans(dna_rel[gi$candidate, ], na.rm = TRUE),
                  dna_ctrl = colMeans(dna_rel[gi$control, ], na.rm = TRUE),
                  B_cand = colMeans(share$B[gi$candidate, ], na.rm = TRUE), B_ctrl = colMeans(share$B[gi$control, ], na.rm = TRUE),
                  T_cand = colMeans(share$T[gi$candidate, ], na.rm = TRUE), T_ctrl = colMeans(share$T[gi$control, ], na.rm = TRUE),
                  my_cand = colMeans(share$myeloid[gi$candidate, ], na.rm = TRUE), my_ctrl = colMeans(share$myeloid[gi$control, ], na.rm = TRUE),
                  n_cand = colSums(is.finite(rel[gi$candidate, ])),
                  n_ctrl = colSums(is.finite(rel[gi$control, ])))
num <- setdiff(names(res), c("band", "r_um", "n_cand", "n_ctrl"))
res[, (num) := lapply(.SD, round, 4), .SDcols = num]
fwrite(res, file.path(OUT, "wide_profile.csv"))
cat("\n== log2 CD45 vs own 400-500 um; diff = candidate - control; diff_adj = adjusted for the band's defect fraction\n")
print(res[, .(band, cand, ctrl, diff, diff_lo, diff_hi, defect_cand, defect_ctrl, diff_adj, adj_lo, adj_hi, slope_defect)])
cat("\n== lineage composition per band (share of step-16-labelled, defect-free cells)\n")
print(res[, .(band, B_cand, B_ctrl, T_cand, T_ctrl, my_cand, my_ctrl)])
## does the CD45 deficit survive adjustment for the band's B share? (cell-level, FOV bootstrap)
adjB <- function(rows) sapply(seq_len(NB), function(b) {
  y <- rel[rows, b]; x <- share$B[rows, b]; gg <- as.integer(pr$group[rows] == "candidate")
  ok <- is.finite(y) & is.finite(x); if (sum(ok) < 10 || sd(x[ok]) == 0) return(NA_real_)
  coef(lm(y[ok] ~ gg[ok] + x[ok]))[2] })
A0 <- adjB(seq_len(nrow(pr)))
AB <- do.call(rbind, lapply(seq_len(N_BOOT), function(b) adjB(unlist(idx[as.character(sample(fovs, replace = TRUE))], use.names = FALSE))))
res[, `:=`(diff_adjB = round(A0, 4), adjB_lo = round(apply(AB, 2, quantile, .025, na.rm = TRUE), 4),
           adjB_hi = round(apply(AB, 2, quantile, .975, na.rm = TRUE), 4))]
fwrite(res, file.path(OUT, "wide_profile.csv"))
cat("\n== candidate - control adjusted for the band's B share\n")
print(res[, .(band, diff, diff_adj, diff_adjB, adjB_lo, adjB_hi)])
cat("\n== DNA, same normalisation (if DNA dims like CD45, the zone is tissue density, not CD45)\n")
print(res[, .(band, cd45_cand = cand, dna_cand, cd45_ctrl = ctrl, dna_ctrl, cd45_minus_dna_cand = round(cand - dna_cand, 4))])
## is the candidates' REFERENCE ring itself dimmer? (raw, FOV-paired where possible)
both <- pr[, .(nc = sum(group == "candidate"), nk = sum(group == "control")), by = fov_id][nc > 0 & nk > 0]$fov_id
pr[, ref_ := ref]
w <- pr[fov_id %in% both, .(lr = log2(median(ref_[group == "candidate"]) / median(ref_[group == "control"]))), by = fov_id][is.finite(lr)]
if (nrow(w) >= 3) {
  bt <- replicate(2000, median(sample(w$lr, replace = TRUE)))
  rr <- data.table(n_fov = nrow(w), median_lr = median(w$lr), lo = quantile(bt, .025), hi = quantile(bt, .975))
} else rr <- data.table(n_fov = nrow(w), median_lr = NA_real_, lo = NA_real_, hi = NA_real_)
cat(sprintf("\nraw 400-500 um CD45, within-FOV log2(candidate / control), %d FOVs: %.3f [%.3f, %.3f]\n",
            rr$n_fov, rr$median_lr, rr$lo, rr$hi))
fwrite(rr, file.path(OUT, "reference_ring_within_fov.csv"))

## ---- 2. figures ----------------------------------------------------------------
L <- rbind(res[, .(group = "candidate", r_um, v = cand, lo = cand_lo, hi = cand_hi, d = defect_cand, dna = dna_cand)],
           res[, .(group = "control", r_um, v = ctrl, lo = ctrl_lo, hi = ctrl_hi, d = defect_ctrl, dna = dna_ctrl)])
p1 <- ggplot(L, aes(r_um, v, colour = group, fill = group)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) +
  geom_line(linewidth = .9) + geom_point(size = 1.3) +
  scale_colour_manual(values = GCOL, name = NULL) + scale_fill_manual(values = GCOL, name = NULL) +
  scale_x_continuous(breaks = seq(0, R_MAX, 50)) +
  labs(title = "CD45 out to 500 um around candidate HRS cells",
       subtitle = "log2 of each 25 um band's mean CD45 over the cell's own 400-500 um mean; 95% FOV cluster bootstrap",
       x = NULL, y = "log2 CD45 vs own 400-500 um") + theme_viz(10, TRUE)
p1b <- ggplot(L, aes(r_um, dna, colour = group)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_line(linewidth = .9) + geom_point(size = 1.3) +
  scale_colour_manual(values = GCOL, guide = "none") +
  scale_x_continuous(breaks = seq(0, R_MAX, 50)) +
  labs(subtitle = "DNA, same normalisation -- tissue / nuclear density", x = NULL, y = "log2 DNA vs own 400-500 um") +
  theme_viz(10, TRUE)
Lc <- rbind(res[, .(group = "candidate", r_um, B = B_cand, T = T_cand, myeloid = my_cand)],
            res[, .(group = "control", r_um, B = B_ctrl, T = T_ctrl, myeloid = my_ctrl)])
Lc <- melt(Lc, id.vars = c("group", "r_um"), variable.name = "lineage", value.name = "share")
p1c <- ggplot(Lc, aes(r_um, share, colour = group, linetype = lineage)) +
  geom_line(linewidth = .8) +
  scale_colour_manual(values = GCOL, guide = "none") + scale_linetype_manual(values = c(B = "solid", T = "22", myeloid = "12"), name = NULL) +
  scale_x_continuous(breaks = seq(0, R_MAX, 50)) + scale_y_continuous(labels = scales::percent) +
  labs(subtitle = "Lineage composition: share of step-16-labelled defect-free cells in the band", x = NULL, y = "Share") +
  theme_viz(10, TRUE) + theme(legend.position = "right")
p2 <- ggplot(L, aes(r_um, d, colour = group)) +
  geom_line(linewidth = .9) + geom_point(size = 1.3) +
  scale_colour_manual(values = GCOL, guide = "none") +
  scale_x_continuous(breaks = seq(0, R_MAX, 50)) + scale_y_continuous(labels = scales::percent) +
  labs(subtitle = "Over-fixation defect: share of all tier cells in the band that are InSituType cluster i",
       x = NULL, y = "Defect fraction") + theme_viz(10, TRUE)
D <- rbind(res[, .(kind = "unadjusted", r_um, v = diff, lo = diff_lo, hi = diff_hi)],
           res[, .(kind = "adjusted for band defect fraction", r_um, v = diff_adj, lo = adj_lo, hi = adj_hi)],
           res[, .(kind = "adjusted for band B share", r_um, v = diff_adjB, lo = adjB_lo, hi = adjB_hi)])
p3 <- ggplot(D, aes(r_um, v, colour = kind)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = .4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 6, linewidth = .4, position = position_dodge(width = 8)) +
  geom_line(linewidth = .8, position = position_dodge(width = 8)) +
  geom_point(size = 1.4, position = position_dodge(width = 8)) +
  scale_colour_manual(values = c(unadjusted = "#0b3d7a", `adjusted for band defect fraction` = "#d03b3b",
                                 `adjusted for band B share` = "#2a9d8f"), name = NULL) +
  scale_x_continuous(breaks = seq(0, R_MAX, 50)) +
  labs(subtitle = "Candidate minus control, per band; adjusted = group coefficient of log2 rel ~ group + defect fraction",
       x = "Distance from centroid (um)", y = "Difference (log2)") + theme_viz(10, TRUE)
fig <- p1 / p1b / p1c / p2 / p3 + plot_layout(heights = c(2, 1, 1.2, 1, 1.4)) +
  plot_annotation(caption = sprintf(
    "Step 31's cells: %d candidates (CCL17+CCL22 >= 50), %d area-matched controls (< 10). Pyramid level 3 (~0.96 um/px). Bands with > 5%% blank image are NA. Source: output/32_wide_cd45/wide_profile.csv.",
    length(gi$candidate), length(gi$control)),
    theme = theme(plot.caption = element_text(colour = MUTED, size = 7, hjust = 0), plot.background = element_rect(fill = SURF, colour = NA)))
sf("fig32_wide_profile.png", fig, 10, 14.5)
cat("\ndone:", OUT, "\n")
