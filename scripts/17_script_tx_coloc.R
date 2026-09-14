## ============================================================================
## 17_script_tx_coloc.R -- the PD-1 axis measured on TRANSCRIPTS, not cells.
##
## WHY THIS IS THE ONE READOUT NOT BOUNDED BY THE DETECTION FLOOR. Every cell-level analysis of
## this axis had to threshold: PDCD1/CD274 sit at 2.6-2.9x the blank-probe rate at >=1 count, so
## step 16 used >=2, which is specific but leaves only ~1-2% of cells positive (~3,000 cells).
## A transcript needs no threshold and no cell boundary. The extraction found 189,651 PDCD1 and
## 157,644 CD274 transcripts in IB1 -- two orders of magnitude more signal than the positive-cell
## count, and none of it discarded by segmentation.
##
## TWO READOUTS:
##   1. regional abundance    PDCD1 (and ligand) transcripts per 10,000 transcripts of the
##                           region. A share, so it is depth- and cellularity-controlled by
##                           construction, exactly as step 11's bin score is.
##   2. inter-cell proximity  median distance from each PDCD1 transcript to the nearest ligand
##                           transcript IN A DIFFERENT CELL, against a null.
##
## SAME-CELL PAIRS ARE EXCLUDED, and that is a biological decision rather than a technical one:
## PD-1 and PD-L1 on the same cell is not the interaction. The interaction is a PD-1-bearing
## T cell touching a PD-L1-bearing macrophage or DC, so the statistic must be between cells.
##
## THE NULL PERMUTES WHICH CELL EACH TRANSCRIPT SITS IN, within cell type and within region.
## That preserves (a) how many transcripts of each gene exist, (b) which cell TYPES carry them
## -- PDCD1 is a T-cell transcript and CD274 a myeloid one, and their patterns differ for
## reasons that have nothing to do with interaction -- and (c) the tissue geometry. It
## randomises only WHICH PARTICULAR CELL, which is the question. A label permutation across all
## transcripts would destroy (b) and manufacture a result.
##
## Transcripts without a cell (CellComp empty, 2.6% in IB1) are reported separately rather than
## spatially reassigned: assigning them would need another coordinate transform and they cannot
## change a conclusion at that share.
##
## Usage: Rscript scripts/17_script_tx_coloc.R [tier]
## ============================================================================

suppressPackageStartupMessages({
  library(RANN); library(data.table); library(ggplot2); setDTthreads(1)
})
`%||%` <- function(a,b) if (is.null(a)) b else a
if (!file.exists("scripts/lib/cosmx_tiff.R")) stop("run from the project root", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
TIER <- if (length(args) >= 1) args[1] else "matched"

TX_F  <- "output/17_pd1/tx/%s_tx_pd1.csv.gz"
REG_F <- file.path("output/16_regions", TIER, "cell_regions.csv.gz")
CT_F  <- "output/16_celltypes/cell_celltypes.csv.gz"
OUT   <- file.path("output/17_pd1", TIER)
FIG   <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

UMPX   <- 120.280945e-6 * 1000     # COSMX_MM_PER_PX in um; tx coords are canvas px
N_PERM <- 100L
## Two caps that cost nothing statistically and everything computationally. v3 timed out at 8h
## because it did an R-level loop over every carrier cell inside every permutation, and IB1
## background alone has 28,346 PDCD1 transcripts over thousands of carriers.
##   MAX_A  the statistic is a MEDIAN nearest-neighbour distance. A median is estimated to
##          within ~1% from a few thousand draws; using all 28,346 anchors buys precision the
##          permutation envelope makes irrelevant, at ~6x the nn2 cost.
##   the sampling loop is vectorised below -- one rep()/sample() over a flattened position
##          index rather than one lapply call per carrier cell.
MAX_A  <- 5000L
MIN_TX <- 200L        # transcripts of the focal gene needed before a region is tested
SEED   <- 1L
SAMPLES<- c("IB1","RB1")
FOC <- list(IB1 = c("FL_focus","cHL_hotspot","background"),
            RB1 = c("FL_follicular","interfollicular"))
## Test pairs, plus CONTROLS whose answer is known. A null that cannot be checked is not a
## null. PDCD1->CD3E must come out CLOSER than chance -- PD-1 is a T-cell transcript, so its
## neighbourhood is T-cell-rich by construction -- and PDCD1->COL1A1 should be neutral or
## further, fibroblasts being elsewhere. If the controls do not behave, the null is broken and
## the test rows mean nothing, which is exactly what happened to the first version.
PAIRS <- list(c("PDCD1","CD274"), c("PDCD1","PDCD1LG2"), c("PDCD1","CXCL13"),
              c("PDCD1","CD3E"),    # POSITIVE control: expect closer than chance
              c("PDCD1","MS4A1"))   # NEGATIVE-ish control: B cells, expect neutral/further
CTRL_POS <- "CD3E"; CTRL_NEG <- "MS4A1"
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"
SAMP_COLS <- c(IB1="#2a78d6", RB1="#eb6834")

cat("============================================================\n")
cat("step 17 -- transcript-level PD-1 axis, tier = ", TIER, "\n", sep="")
cat("============================================================\n")

rg <- fread(REG_F, select=c("cell_id","samp","bin_region","is_defect","unbinned","fov_id"))
ct <- fread(CT_F, select=c("cell_id","celltype"))
cells <- merge(rg, ct, by="cell_id", sort=FALSE)[is_defect==FALSE & unbinned==FALSE]
cells <- cells[mapply(function(s,r) r %in% FOC[[s]], samp, bin_region)]
setkey(cells, cell_id)
cat(sprintf("cells in scope: %d\n", nrow(cells)))

TX <- rbindlist(lapply(SAMPLES, function(S) {
  f <- sprintf(TX_F, S); if (!file.exists(f)) stop("missing ", f, call.=FALSE)
  x <- fread(f); x[, samp := S][]
}))
cat(sprintf("transcripts read: %d\n", nrow(TX)))
TX[, has_cell := !is.na(cell) & nzchar(cell) & CellComp != ""]
cat(sprintf("  with a cell assignment: %d (%.1f%%)\n", TX[, sum(has_cell)],
            100*TX[, mean(has_cell)]))
cat("\n=== unassigned transcripts by target (reported, not reassigned)\n")
print(TX[has_cell == FALSE, .(unassigned = .N), by = .(samp, target)][order(samp, -unassigned)][1:8])

## join to regions via the cell id; transcripts in out-of-scope cells drop out
n0 <- nrow(TX)
TX <- merge(TX[has_cell == TRUE], cells, by.x = "cell", by.y = "cell_id", sort = FALSE)
if ("samp.y" %in% names(TX)) { TX[, samp := samp.x][, c("samp.x","samp.y") := NULL] }
cat(sprintf("\ntranscripts in in-scope cells: %d of %d\n", nrow(TX), n0))
print(TX[, .(transcripts = .N, cells = uniqueN(cell)), by = .(samp, bin_region)][order(samp,-transcripts)])

## ---- 1. regional abundance, as a share of the region's transcripts ----------
tot <- TX[, .(tot_tx = .N), by = .(samp, bin_region)]
ab  <- TX[, .(n = .N), by = .(samp, bin_region, target)]
ab  <- merge(ab, tot, by = c("samp","bin_region"))
ab[, per10k := round(1e4 * n / tot_tx, 2)]
sc("tx_abundance.csv", ab)
cat("\n=== transcripts per 10,000 region transcripts (no threshold, no cell calling)\n")
print(dcast(ab[target %in% c("PDCD1","CD274","PDCD1LG2","CXCL13","CXCL9","IFNG")],
            target ~ samp + bin_region, value.var = "per10k"))

## ---- 2. inter-cell proximity against a cell-identity-permuted null ----------
## THE NULL, THIRD ATTEMPT. The controls have now caught two distinct biases, which is what
## they are for; the stopping rule is stated at the bottom of this block.
##
##   v1 compared observed transcript-to-TRANSCRIPT distances against a null of
##      transcript-to-CELL-CENTROID distances. Centroids are more regularly spaced, so the null
##      was tighter and everything read "further than chance".
##   v2 drew pseudo-transcripts from other genes' transcripts, matched on cell type. That fixed
##      the position mismatch but introduced a worse one: REAL transcripts are OVER-DISPERSED
##      across cells -- a handful of cells carry many, most carry none -- while a donor pool
##      spreads them thinly over every cell of the type. A thin uniform spread puts a
##      pseudo-transcript nearer to any given A than the real clustered pattern does, so again
##      everything read "further".
##
## v3 preserves the over-dispersion explicitly. The unit resampled is the CARRIER CELL, not the
## transcript: take the set of cells that actually carry B and their per-cell counts, draw the
## same number of donor cells matched on cell type, and hand each donor the count of one real
## carrier. Positions come from that donor cell's own transcripts, so within-cell placement
## stays realistic. What is preserved: transcript count, cell-type usage, the number of
## carrying cells, the distribution of counts per carrier, and tissue geometry. What is
## randomised: which cells are the carriers -- the question.
##
## STOPPING RULE, committed before looking: if PDCD1->CD3E still fails to come out closer than
## chance in most regions, this measurement is not working and the proximity analysis is
## reported as NOT DELIVERED. The transcript ABUNDANCE readout above needs no null and stands
## on its own regardless.
prox <- function(S, R, ga, gb) {
  z <- TX[samp == S & bin_region == R]
  A <- z[target == ga]; B <- z[target == gb]
  if (nrow(A) < MIN_TX || nrow(B) < MIN_TX) return(NULL)
  set.seed(SEED)
  if (nrow(A) > MAX_A) A <- A[sample.int(.N, MAX_A)]
  ax <- as.matrix(A[, .(x_global_px, y_global_px)])
  nn_diff <- function(bx, bcell) {
    k <- min(12L, nrow(bx)); res <- rep(NA_real_, nrow(A))
    repeat {
      q <- RANN::nn2(bx, ax, k = k)
      for (i in seq_len(k)) {
        hit <- is.na(res) & bcell[q$nn.idx[, i]] != A$cell
        res[hit] <- q$nn.dists[hit, i]
      }
      if (!anyNA(res) || k >= nrow(bx)) break
      k <- min(k * 4L, nrow(bx))
    }
    res * UMPX
  }
  obs <- median(nn_diff(as.matrix(B[, .(x_global_px, y_global_px)]), B$cell), na.rm = TRUE)

  ## carriers of B, their counts, and their cell types
  carr <- B[, .(n = .N), by = .(cell, celltype)]
  ## position candidates per cell: that cell's own transcripts of any extracted gene
  pos <- z[, .(cell, celltype, x_global_px, y_global_px)]
  poscell <- split(seq_len(nrow(pos)), pos$cell)
  pflat  <- unlist(poscell, use.names = FALSE)
  plen   <- lengths(poscell)
  pstart <- c(0L, cumsum(plen))[seq_along(plen)]   # 0-based offset into pflat
  ## donor cells per type: any cell in the region with at least one transcript
  cellsbytype <- unique(pos[, .(cell, celltype)])
  donpool <- split(cellsbytype$cell, cellsbytype$celltype)
  need <- carr[, .(ncell = .N), by = celltype]
  if (any(need$ncell > vapply(need$celltype, function(t)
        length(donpool[[as.character(t)]]) %||% 0L, 0L))) {
    ## not enough donor cells of some type; fall back to sampling with replacement
  }
  set.seed(SEED)
  nul <- vapply(seq_len(N_PERM), function(i) {
    picks <- rbindlist(lapply(seq_len(nrow(need)), function(j) {
      tt <- as.character(need$celltype[j]); pool <- donpool[[tt]]
      if (is.null(pool) || !length(pool)) return(NULL)
      cnt <- carr[celltype == need$celltype[j], n]
      dc  <- sample(pool, length(cnt), replace = length(pool) < length(cnt))
      data.table(cell = dc, n = sample(cnt))
    }), fill = TRUE)
    if (!nrow(picks)) return(NA_real_)
    ## Vectorised draw. `pstart`/`plen` are a flattened ragged index over per-cell position
    ## rows, so every carrier's transcripts are drawn in one runif() rather than one
    ## sample() call each -- the loop this replaces was the whole reason v3 timed out.
    pk <- match(picks$cell, names(poscell))
    ok <- !is.na(pk) & plen[pk] > 0L
    if (!any(ok)) return(NA_real_)
    reps <- rep.int(pk[ok], picks$n[ok])
    idx  <- pflat[pstart[reps] + floor(runif(length(reps)) * plen[reps])]
    if (length(idx) < MIN_TX) return(NA_real_)
    median(nn_diff(as.matrix(pos[idx, .(x_global_px, y_global_px)]), pos$cell[idx]), na.rm = TRUE)
  }, 0)
  nul <- nul[is.finite(nul)]
  if (!length(nul)) return(NULL)
  data.table(samp=S, bin_region=R, a=ga, b=gb, n_a=nrow(A), n_b=nrow(B),
             n_anchors_used=nrow(A),
             n_carriers=nrow(carr), med_per_carrier=as.numeric(median(carr$n)),
             obs_um=round(obs,2), null_um=round(median(nul),2), null_sd=round(sd(nul),3),
             ratio_null_obs=round(median(nul)/obs,3),
             p_emp=(1 + sum(nul <= obs))/(length(nul)+1))
}
PX <- rbindlist(lapply(SAMPLES, function(S) {
  rbindlist(lapply(FOC[[S]], function(R) {
    rbindlist(lapply(PAIRS, function(pp) prox(S, R, pp[1], pp[2])), fill=TRUE)
  }), fill=TRUE)
}), fill=TRUE)
if (nrow(PX)) {
  PX[, verdict := fifelse(p_emp < .05 & ratio_null_obs > 1, "closer than chance",
                  fifelse(p_emp > .95 & ratio_null_obs < 1, "further than chance", "no signal"))]
  sc("tx_proximity.csv", PX)
  cat("\n=== inter-cell transcript proximity vs a cell-identity-permuted null\n")
  cat("    ratio > 1 with p < 0.05 = the two transcripts sit closer than cell type and\n")
  cat("    geometry alone predict. Same-cell pairs excluded throughout.\n")
  print(PX[order(a, b, samp, bin_region)])
  ## IS THE NULL CALIBRATED? Judged on the controls, not on the test pairs.
  cat("\n=== null calibration, from the control pairs\n")
  pc <- PX[b == CTRL_POS]; nc <- PX[b == CTRL_NEG]
  cat(sprintf("  POSITIVE control PDCD1->%s: %d of %d regions closer than chance (expect most)\n",
              CTRL_POS, pc[verdict == "closer than chance", .N], nrow(pc)))
  cat(sprintf("  NEGATIVE control PDCD1->%s: %d of %d closer, %d further (expect few closer)\n",
              CTRL_NEG, nc[verdict == "closer than chance", .N], nrow(nc),
              nc[verdict == "further than chance", .N]))
  ok <- nrow(pc) && pc[verdict == "closer than chance", .N] >= ceiling(nrow(pc)/2)
  cat(if (ok)
    "  -> the positive control behaves; the test rows can be read.\n" else
    paste0("  -> *** THE POSITIVE CONTROL FAILS. The null is still miscalibrated and NO row\n",
           "  *** in the table above should be interpreted. Do not report these.\n"))
  if (PX[, uniqueN(verdict)] == 1L)
    cat("  *** every row shares one verdict -- that is a broken null, not a finding.\n")
} else cat("\n=== no region had enough transcripts of both genes\n")

## ---- figure ------------------------------------------------------------------
pl <- ab[target %in% c("PDCD1","CD274","PDCD1LG2")]
if (nrow(pl)) {
  pl[, reg := paste(samp, bin_region, sep="\n")]
  p <- ggplot(pl, aes(reg, per10k, fill = samp)) +
    geom_col(width = .68) +
    geom_text(aes(label = per10k), vjust = -0.35, size = 2.6, colour = INK2) +
    facet_wrap(~ target, scales = "free_y") +
    scale_fill_manual(values = SAMP_COLS, name = NULL) +
    labs(x = NULL, y = "transcripts per 10,000 region transcripts",
         title = "PD-1 axis measured on transcripts, not cells",
         subtitle = paste0("No positivity threshold and no cell calling, so this readout is ",
           "not bounded by the detection floor that\nlimited every cell-level test. A share of ",
           "the region's transcripts, so cellularity and depth are controlled."),
         caption = "Data: output/17_pd1/<tier>/tx_abundance.csv") +
    theme_minimal(base_size = 10) +
    theme(plot.background=element_rect(fill=SURF,colour=NA),
          panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(colour=GRID, linewidth=.25),
          panel.grid.minor=element_blank(),
          axis.text.x=element_text(colour=MUTED, size=6.5),
          axis.text.y=element_text(colour=MUTED, size=7),
          axis.title=element_text(colour=INK2, size=9),
          strip.text=element_text(colour=INK, face="bold", size=9),
          plot.title=element_text(colour=INK, face="bold", size=12),
          plot.subtitle=element_text(colour=INK2, size=8),
          plot.caption=element_text(colour=MUTED, size=7, hjust=0),
          legend.position="top", legend.justification="left")
  ggsave(file.path(FIG,"fig_17_tx_abundance.png"), p, width=10, height=5, dpi=170, bg=SURF)
  cat("wrote ", file.path(FIG,"fig_17_tx_abundance.png"), "\n", sep="")
}
writeLines(format(Sys.time()), file.path(OUT, "tx_coloc.done"))
cat("============================================================\nOK\n")
