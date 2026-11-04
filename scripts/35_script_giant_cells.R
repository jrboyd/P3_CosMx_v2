## ============================================================================
## Step 35 -- Giant-cell reconstruction: are chemokine-bright polygons fragments of
## larger cells? Composites of touching chemokine-positive polygons, under two
## partitions, against a within-FOV permutation null.
##
## WHY. Step 34: candidate polygons are pale-filled or DNA-dim at the centroid, bright
## nuclei press against them from outside, and the HRS-leaning RNA sits in dense-nucleus
## cells while the pale large nuclei are RNA-quiet. A 20-60 um HRS cell would be carved
## into several ~100 um^2 vendor polygons -- one holding the nucleus, the rest cytoplasm
## carrying the mRNA. If so, the chemokine-bright polygons should be CONTIGUOUS beyond
## chance, the contiguous sets should hold few nuclei for their size, and pooling them
## should reunite the morphology with the RNA.
##
## STAGES (scripts/35_submit_giant_cells.sh)
##   select   (this, R)   cells of every FOV holding a seed, with counts + vendor metrics;
##                        FOV canvas origins; per-FOV polygon files split by awk
##   adjacency (35_fov_adjacency.py, SLURM array, one FOV per task)  polygons rasterised to
##                        a label image; touching pairs with contact length; per polygon
##                        dense / pale DNA fractions from the full-resolution DNA channel
##   analyze  (this, R)   partitions, null, tests, figures
##
## DESIGN (fixed before running)
##   seeds        IB1: CCL17+CCL22 >= 50; RB1 (no HRS -- what activated clusters look
##                like): >= 20, its own brightest. Positivity thresholds 10 / 20 / 50.
##   touching     contact >= MIN_CONTACT_PX pixels at a 1-3 px shift (corners excluded)
##   partition 1  CONNECTED COMPONENTS of touching positive polygons. Agnostic; can chain
##                two objects into one. Chains show as >= 2 dense nuclei.
##   partition 2  NUCLEUS-ANCHORED: positive polygons with a dense nucleus (DNA >= FOV
##                Otsu over >= 8 um^2) are anchors; nucleus-free positive polygons are
##                assigned by label propagation to the touching anchor-composite with the
##                most contact. One nucleus per composite by construction. A fragment with
##                >= 2 candidate parents at assignment is AMBIGUOUS; counted, and the RNA
##                test is run three ways (best parent / counts split by contact / dropped).
##   null         N_PERM permutations of the EXPRESSION ROWS among a FOV's polygons --
##                geometry, nuclei and per-FOV abundance kept, spatial arrangement broken;
##                both partitions rebuilt each time. p = (1 + #null >= obs) / (1 + N_PERM).
##   tests        (1) contiguity: members and area of seed composites vs null;
##                (2) dense nuclei per member in >= 2-member seed composites vs null
##                    composites of the same member count -- the giant-cell vs cluster
##                    discriminator; (3) RNA: nucleus member vs nucleus-free members of the
##                    same composite, excess per gene vs < 10 cells of the same FOV set and
##                    depth decile; (4) agreement of the two partitions; chain rate.
##   unit         FOV: the null is within-FOV; intervals are FOV cluster bootstraps.
##   seams        composites with a member < SEAM_UM from the FOV edge are flagged and
##                every statistic is reported with and without them.
##
##   scripts/35_submit_giant_cells.sh [smoke]
## ============================================================================
suppressPackageStartupMessages(library(data.table))
args  <- commandArgs(trailingOnly = TRUE)
MODE  <- args[1]
SMOKE <- "smoke" %in% args
OUT   <- if (SMOKE) "output/35_giant_cells_smoke" else "output/35_giant_cells"
FIG   <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
FLAT  <- "input/CosMx/export1/flatFiles"
SEED  <- c(IB1 = 50L, RB1 = 20L)
THR   <- c(10L, 20L, 50L)
SLIDE <- c(IB1 = 2L, RB1 = 1L)        # c_2_ = IB1 (CLAUDE.md); checked against samp
MIN_CONTACT_PX <- 3L
SEAM_UM <- 5
NUC_MIN_UM2 <- 8
UM2_PER_PX2 <- 0.120280945^2
N_PERM <- if (SMOKE) 10L else 200L
N_BOOT <- if (SMOKE) 20L else 1000L
if (SMOKE) options(warn = 1)
set.seed(35)

SETS <- list(
  HRS       = c("TNFRSF8", "FUT4", "IL13", "LTA", "CSF2"),
  ambiguous = c("CD274", "PDCD1LG2", "CD40", "BATF3", "IRF4", "CD80", "CD86", "CD83", "JAK2"),
  DC        = c("LAMP3", "CCR7", "CCL19", "IDO1", "CD1C", "CLEC9A", "XCR1", "CLEC10A",
                "FCER1A", "CD1A", "CD207", "MARCKSL1"),
  macro     = c("CD68", "CD163", "C1QA", "CD14", "LYZ"),
  Bprog     = c("CD19", "CD79A", "MS4A1"),
  T         = c("CD3E", "CD3D"))
TEST <- unlist(SETS, use.names = FALSE); GENES <- c("CCL17", "CCL22", TEST)
SET_OF <- setNames(rep(names(SETS), lengths(SETS)), TEST)

## ============================================================================ select
if (MODE == "select") {
  source("scripts/lib/cosmx_tiff.R")
  src <- fread("output/29_bright_cells/source_cells.csv.gz")
  reg <- fread(cmd = "zcat output/16_regions/matched/cell_regions.csv.gz",
               select = c("cell_id", "samp", "is_defect", "tier_matched", "nCount_RNA"))
  fl <- rbindlist(lapply(names(SEED), function(s) {
    f <- sort(unique(src[samp == s & chemo >= SEED[[s]]]$fov))
    if (SMOKE) f <- head(f, 6)
    pl <- cosmx_fov_placement(s)
    data.table(samp = s, fov = f)[as.data.table(pl), on = "fov", nomatch = 0][, .(samp, fov, x0, y0)]
  }))
  fl[, task := seq_len(.N) - 1L]
  fwrite(fl, file.path(OUT, "fov_list.csv"))
  cat(sprintf("FOVs: %s\n", paste(fl[, .N, by = samp][, sprintf("%s %d", samp, N)], collapse = ", ")))
  for (s in names(SEED)) {
    fovs <- fl[samp == s]$fov
    ff <- file.path(OUT, sprintf("fovs_%s.txt", s)); writeLines(as.character(fovs), ff)
    m <- fread(file.path(FLAT, s, paste0(s, "_metadata_file.csv.gz")),
               select = c("fov", "cell_ID", "Area", "NucArea", "Mean.DAPI", "CenterX_local_px", "CenterY_local_px"))[fov %in% fovs]
    f <- file.path(FLAT, s, paste0(s, "_exprMat_file.csv.gz"))
    hdr <- names(fread(cmd = sprintf("zcat %s | head -1", f)))
    stopifnot(all(GENES %in% hdr))
    ex <- fread(cmd = sprintf("zcat %s | awk -F, 'NR==FNR{k[$1]; next} FNR==1 || ($1 in k)' %s - | cut -d, -f%s",
                              f, ff, paste(match(c("fov", "cell_ID", GENES), hdr), collapse = ",")))
    d <- merge(m, ex, by = c("fov", "cell_ID"))
    d[, cell_id := sprintf("c_%d_%d_%d", SLIDE[[s]], fov, cell_ID)]
    d <- merge(d, reg[samp == s, .(cell_id, is_defect, tier_matched, nCount_RNA)], by = "cell_id", all.x = TRUE)
    d[, `:=`(samp = s, chemo = CCL17 + CCL22, in_tier = !is.na(tier_matched) & tier_matched == TRUE)]
    d[, depth0 := nCount_RNA - chemo]
    fwrite(d, file.path(OUT, sprintf("cells_%s.csv.gz", s)))
    cat(sprintf("%s: %d polygons in %d FOVs; seeds (>= %d, defect-free tier) %d\n", s, nrow(d), length(fovs), SEED[[s]],
                d[in_tier & is_defect == FALSE & chemo >= SEED[[s]], .N]))
  }
  quit(save = "no")
}

## ============================================================================ analyze
stopifnot(MODE == "analyze")
suppressPackageStartupMessages({ library(igraph); library(ggplot2); library(patchwork) })
SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
theme_viz <- function(base=11, xgrid=FALSE) {
  th <- theme_minimal(base_size=base) +
    theme(plot.background=element_rect(fill=SURF,colour=NA), panel.background=element_rect(fill=SURF,colour=NA),
          panel.grid.major.y=element_line(colour=GRID,linewidth=.25), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
          axis.line.x=element_line(colour=AXIS,linewidth=.3), axis.ticks=element_blank(),
          axis.text=element_text(colour=MUTED,size=base-2), axis.title=element_text(colour=INK2,size=base-1),
          plot.title=element_text(colour=INK,face="bold",size=base+2), plot.subtitle=element_text(colour=INK2,size=base-1),
          plot.caption=element_text(colour=MUTED,size=base-3,hjust=0), legend.text=element_text(colour=INK2,size=base-2),
          legend.title=element_text(colour=INK2,size=base-2), legend.position="top", legend.justification="left", legend.key.height=unit(9,"pt"))
  if (xgrid) th <- th + theme(panel.grid.major.x=element_line(colour=GRID,linewidth=.25))
  th
}
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF);cat("wrote",file.path(FIG,n),"\n")}
SAMP_COLS <- c(IB1 = "#2a78d6", RB1 = "#eb6834")
pboot <- function(x) { x <- x[is.finite(x)]; n <- length(x); min(1, 2 * min((1 + sum(x <= 0)) / (1 + n), (1 + sum(x >= 0)) / (1 + n))) }
pperm <- function(obs, nul, greater = TRUE) { nul <- nul[is.finite(nul)]; (1 + sum(if (greater) nul >= obs else nul <= obs)) / (1 + length(nul)) }

fl <- fread(file.path(OUT, "fov_list.csv"))
load_sample <- function(s) {
  d <- fread(file.path(OUT, sprintf("cells_%s.csv.gz", s)))
  fovs <- fl[samp == s]$fov
  pc <- rbindlist(lapply(fovs, function(f) { p <- file.path(OUT, "adj", sprintf("%s_%d_cells.csv", s, f)); if (file.exists(p)) fread(p) }))
  ed <- rbindlist(lapply(fovs, function(f) { p <- file.path(OUT, "adj", sprintf("%s_%d_edges.csv", s, f)); if (file.exists(p)) fread(p) }))
  d <- merge(d, pc[, .(fov, cell_ID, poly_px, dense_px, pale_px, dna_rel, edge_um, fov_blank)], by = c("fov", "cell_ID"))
  d[, `:=`(poly_um2 = poly_px * UM2_PER_PX2, dense_um2 = dense_px * UM2_PER_PX2,
           has_nuc = !is.na(dense_px) & dense_px * UM2_PER_PX2 >= NUC_MIN_UM2)]   # degenerate polygons: no nucleus
  d <- d[!is.na(chemo)]
  d[, id := seq_len(.N)]
  ed <- ed[contact_px >= MIN_CONTACT_PX]
  ed <- merge(ed, d[, .(fov, a = cell_ID, ia = id)], by = c("fov", "a"))
  ed <- merge(ed, d[, .(fov, b = cell_ID, ib = id)], by = c("fov", "b"))
  list(d = d, ed = ed[, .(fov, ia, ib, contact_px)])
}

## partitions -------------------------------------------------------------------
## pos: logical per polygon (positive at a threshold). Returns per-polygon composite ids
## for partition 1 (comp) and 2 (anch), plus per-fragment parent ambiguity.
partition <- function(d, ed, pos) {
  n <- nrow(d)
  e <- ed[pos[ia] & pos[ib]]
  comp <- rep(NA_integer_, n)
  if (nrow(e)) {
    g <- graph_from_edgelist(as.matrix(e[, .(ia, ib)]), directed = FALSE)
    cl <- components(g)$membership
    comp[seq_along(cl)] <- cl
  }
  comp[pos & is.na(comp)] <- max(c(0L, comp), na.rm = TRUE) + seq_len(sum(pos & is.na(comp)))   # singletons
  comp[!pos] <- NA_integer_
  ## anchored: anchors = pos & has_nuc; propagate to pos & !has_nuc by max contact
  anch <- rep(NA_integer_, n); is_anchor <- pos & d$has_nuc
  anch[is_anchor] <- seq_len(sum(is_anchor))
  n_parent <- rep(NA_integer_, n)
  if (nrow(e)) {
    e2 <- rbind(e[, .(from = ia, to = ib, contact_px)], e[, .(from = ib, to = ia, contact_px)])
    repeat {
      cand <- e2[is.na(anch[from]) & pos[from] & !is.na(anch[to])]
      if (!nrow(cand)) break
      cand[, parent := anch[to]]
      agg <- cand[, .(contact = sum(contact_px)), by = .(from, parent)]
      best <- agg[order(from, -contact)][, .(parent = parent[1], n_parent = .N), by = from]
      anch[best$from] <- best$parent; n_parent[best$from] <- best$n_parent
    }
  }
  orphan <- pos & is.na(anch)                       # nucleus-free with no path to an anchor
  anch[orphan] <- max(c(0L, anch), na.rm = TRUE) + seq_len(sum(orphan))
  list(comp = comp, anch = anch, n_parent = n_parent, orphan = orphan)
}

## composite summaries: one row per composite containing >= 1 seed ----------------
summarise_comp <- function(d, cid, seed) {
  k <- !is.na(cid)
  x <- data.table(cid = cid[k], fov = d$fov[k], seed = seed[k], area = d$poly_um2[k], nuc = d$has_nuc[k],
                  chemo = d$chemo[k], edge = d$edge_um[k], tier = d$in_tier[k] & !is.na(d$is_defect[k]) & d$is_defect[k] == FALSE)
  stopifnot(nrow(x) == sum(k))
  x[, .(fov = fov[1], members = .N, area_um2 = sum(area), n_nuc = sum(nuc), chemo = sum(chemo), has_seed = any(seed),
        seam = any(edge < SEAM_UM), n_tier = sum(tier)), by = cid][has_seed == TRUE]
}
stat_of <- function(cs) {
  c(n = nrow(cs), med_members = median(cs$members), frac_multi = mean(cs$members >= 2),
    frac_ge300 = mean(cs$area_um2 >= 300), med_area = median(cs$area_um2),
    nuc_per_member_multi = cs[members >= 2, sum(n_nuc) / sum(members)],
    frac_multi_le1nuc = cs[members >= 2, mean(n_nuc <= 1)], frac_chain = mean(cs$n_nuc >= 2))
}

results <- list(); nulls <- list(); rna_rows <- list(); amb_rows <- list(); comp_tabs <- list(); assign_rows <- list()
for (s in names(SEED)) {
  cat(sprintf("\n==================== %s\n", s))
  L <- load_sample(s)
  stopifnot(all(L$d$samp == s))
  ## drop polygons of FOVs with no DNA image (no nuclei can be scored there), re-index, remap edges
  keep <- L$d$fov_blank == FALSE
  remap <- cumsum(keep) * keep
  d <- L$d[keep]; d[, id := seq_len(.N)]
  ed <- L$ed[keep[ia] & keep[ib]]; ed[, `:=`(ia = remap[ia], ib = remap[ib])]
  d[, seed := in_tier & !is.na(is_defect) & is_defect == FALSE & chemo >= SEED[[s]]]
  cat(sprintf("%d polygons in %d FOVs (%d dropped in blank-DNA FOVs), %d touching pairs (contact >= %d px), %d seeds\n",
              nrow(d), uniqueN(d$fov), sum(!keep), nrow(ed), MIN_CONTACT_PX, sum(d$seed)))
  ## expression columns to permute within FOV (rows travel together)
  XCOL <- c("chemo", GENES, "depth0", "seed", "in_tier", "is_defect")
  fov_idx <- split(seq_len(nrow(d)), d$fov)
  for (t in THR) {
    pos <- !is.na(d$chemo) & d$chemo >= t
    P <- partition(d, ed, pos)
    cs1 <- summarise_comp(d, P$comp, d$seed); cs2 <- summarise_comp(d, P$anch, d$seed)
    cs1[, `:=`(samp = s, thr = t, partition = "components")]; cs2[, `:=`(samp = s, thr = t, partition = "anchored")]
    comp_tabs[[paste(s, t)]] <- rbind(cs1, cs2)
    frag <- pos & !d$has_nuc
    ## per-polygon assignment under both partitions, for the review panels (observed only)
    assign_rows[[paste(s, t)]] <- d[pos, .(samp = s, thr = t, fov, cell_ID, cell_id, chemo, has_nuc, seed, poly_um2, edge_um,
                                           comp = P$comp[id], anch = P$anch[id], n_parent = P$n_parent[id], orphan = P$orphan[id])]
    amb_rows[[paste(s, t)]] <- data.table(samp = s, thr = t, n_pos = sum(pos), n_anchor = sum(pos & d$has_nuc), n_fragment = sum(frag),
      n_orphan = sum(P$orphan), n_ambiguous = sum(frag & !is.na(P$n_parent) & P$n_parent >= 2),
      frac_ambiguous_of_assigned = mean(P$n_parent[frag & !is.na(P$n_parent)] >= 2))
    S1 <- stat_of(cs1); S2 <- stat_of(cs2)
    S1ns <- stat_of(cs1[seam == FALSE]); S2ns <- stat_of(cs2[seam == FALSE])
    ## null: permute expression rows within FOV, rebuild both partitions
    NUL <- lapply(seq_len(N_PERM), function(b) {
      dp <- copy(d)
      for (rows in fov_idx) dp[rows, (XCOL) := d[rows[sample(length(rows))], XCOL, with = FALSE]]
      posp <- !is.na(dp$chemo) & dp$chemo >= t
      Pp <- partition(dp, ed, posp)
      c1 <- summarise_comp(dp, Pp$comp, dp$seed); c2 <- summarise_comp(dp, Pp$anch, dp$seed)
      rbind(data.table(partition = "components", t(stat_of(c1))), data.table(partition = "anchored", t(stat_of(c2))))[, b := b]
    })
    NUL <- rbindlist(NUL)
    nulls[[paste(s, t)]] <- NUL[, `:=`(samp = s, thr = t)]
    for (pn in c("components", "anchored")) {
      S <- if (pn == "components") S1 else S2; Sns <- if (pn == "components") S1ns else S2ns
      nn <- NUL[partition == pn]
      results[[paste(s, t, pn)]] <- data.table(samp = s, thr = t, partition = pn, stat = names(S), obs = S, obs_noseam = Sns,
        null_med = sapply(names(S), function(k) median(nn[[k]], na.rm = TRUE)),
        null_p95 = sapply(names(S), function(k) quantile(nn[[k]], .95, na.rm = TRUE)),
        p_greater = sapply(names(S), function(k) pperm(S[[k]], nn[[k]], TRUE)),
        p_less = sapply(names(S), function(k) pperm(S[[k]], nn[[k]], FALSE)))
    }
    cat(sprintf("\n-- threshold >= %d: positive %d, anchors %d, fragments %d (ambiguous %d, orphan %d)\n", t, sum(pos), sum(pos & d$has_nuc),
                sum(frag), amb_rows[[paste(s, t)]]$n_ambiguous, sum(P$orphan)))
    print(rbindlist(results[grep(paste0("^", s, " ", t, " "), names(results))])[stat %in% c("n", "med_members", "frac_multi", "frac_ge300", "nuc_per_member_multi", "frac_multi_le1nuc", "frac_chain"),
          .(partition, stat, obs = round(obs, 3), obs_noseam = round(obs_noseam, 3), null_med = round(null_med, 3), null_p95 = round(null_p95, 3), p_greater = signif(p_greater, 3), p_less = signif(p_less, 3))])

    ## RNA: nucleus member vs nucleus-free members, anchored partition, seed composites with >= 2 members
    if (t == SEED[[s]] || (s == "IB1" && t == 20L)) {
      ref <- d[in_tier & !is.na(is_defect) & is_defect == FALSE & chemo < 10 & is.finite(depth0) & depth0 > 0]
      ref[, depth_q := cut(depth0, unique(quantile(depth0, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
      brk <- unique(quantile(ref$depth0, 0:10 / 10))
      rates <- ref[, c(list(dep = sum(as.numeric(depth0))), lapply(.SD, function(z) sum(as.numeric(z)))), by = depth_q, .SDcols = TEST]
      for (g in TEST) rates[, (paste0("r_", g)) := get(g) / dep]
      multi <- cs2[members >= 2]$cid
      m <- d[!is.na(P$anch) & P$anch %in% multi & in_tier & !is.na(is_defect) & is_defect == FALSE & is.finite(depth0) & depth0 > 0]
      m[, cid := P$anch[id]]; m[, role := fifelse(has_nuc, "nucleus", "fragment")]
      m[, n_parent := P$n_parent[id]]
      m[, depth_q := cut(depth0, brk, include.lowest = TRUE, labels = FALSE)]
      m <- merge(m, rates[, c("depth_q", paste0("r_", TEST)), with = FALSE], by = "depth_q")
      for (g in TEST) m[, (paste0("e_", g)) := depth0 * get(paste0("r_", g))]
      cat(sprintf("RNA set: %d composites with >= 2 members; %d nucleus members, %d fragments (%d ambiguous)\n",
                  uniqueN(m$cid), sum(m$role == "nucleus"), sum(m$role == "fragment"), sum(m$role == "fragment" & !is.na(m$n_parent) & m$n_parent >= 2)))
      if (uniqueN(m$cid) < 5 || sum(m$role == "fragment") < 3 || uniqueN(m$fov) < 3) {
        cat("too few multi-member composites for the RNA contrast; skipped\n"); next
      }
      ex_of <- function(mm, w = rep(1, nrow(mm))) {
        rbindlist(lapply(c("nucleus", "fragment"), function(r) { i <- mm$role == r
          o <- colSums(as.matrix(mm[i, ..TEST]) * w[i]); e <- colSums(as.matrix(mm[i, paste0("e_", TEST), with = FALSE]) * w[i])
          data.table(role = r, gene = TEST, lr = log((o + 0.5) / (e + 0.5)), obs = o) })) }
      variants <- list(best = m, drop = m[!(role == "fragment" & !is.na(n_parent) & n_parent >= 2)])
      ## fractional: an ambiguous fragment's counts are shared across its candidate parents -> weight 1/n_parent on the composite axis;
      ## for the role contrast that is a down-weighting of ambiguous fragments
      wfrac <- fifelse(m$role == "fragment" & !is.na(m$n_parent) & m$n_parent >= 2, 1 / m$n_parent, 1)
      fovs_m <- unique(m$fov); idx_m <- split(seq_len(nrow(m)), m$fov)
      for (v in c("best", "fractional", "drop")) {
        mm <- if (v == "fractional") m else variants[[v]]
        ww <- if (v == "fractional") wfrac else rep(1, nrow(mm))
        E0 <- ex_of(mm, ww)
        idx_v <- split(seq_len(nrow(mm)), mm$fov); fv <- unique(mm$fov)
        B <- lapply(seq_len(N_BOOT), function(b) { rows <- unlist(idx_v[as.character(sample(fv, replace = TRUE))], use.names = FALSE)
          if (length(rows) < 2) return(NULL); ex_of(mm[rows], ww[rows]) })
        B <- Filter(Negate(is.null), B)
        Bn <- sapply(B, function(x) x[role == "nucleus"]$lr); Bf <- sapply(B, function(x) x[role == "fragment"]$lr)
        rna_rows[[paste(s, t, v)]] <- data.table(samp = s, thr = t, variant = v, gene = TEST, set = SET_OF[TEST],
          obs_nucleus = E0[role == "nucleus"]$obs, obs_fragment = E0[role == "fragment"]$obs,
          excess_nucleus = exp(E0[role == "nucleus"]$lr), excess_fragment = exp(E0[role == "fragment"]$lr),
          ratio_frag_over_nuc = exp(E0[role == "fragment"]$lr - E0[role == "nucleus"]$lr),
          lo = exp(apply(Bf - Bn, 1, quantile, .025, na.rm = TRUE)), hi = exp(apply(Bf - Bn, 1, quantile, .975, na.rm = TRUE)),
          p_boot = apply(Bf - Bn, 1, pboot))
      }
      r <- rna_rows[[paste(s, t, "best")]]
      cat("fragment / nucleus excess ratio (best-parent variant), HRS + DC + gauge:\n")
      print(r[gene %in% c("IL13", "TNFRSF8", "LTA", "BATF3", "IRF4", "CCR7", "LAMP3", "XCR1", "CD3E", "CD3D", "LYZ"),
              .(gene, set, excess_nucleus = round(excess_nucleus, 2), excess_fragment = round(excess_fragment, 2),
                ratio = round(ratio_frag_over_nuc, 2), lo = round(lo, 2), hi = round(hi, 2), p = signif(p_boot, 2))])
    }
  }
}
res <- rbindlist(results); nul <- rbindlist(nulls); rna <- rbindlist(rna_rows); amb <- rbindlist(amb_rows); ct <- rbindlist(comp_tabs)
num <- c("obs", "obs_noseam", "null_med", "null_p95"); res[, (num) := lapply(.SD, round, 4), .SDcols = num]
fwrite(res, file.path(OUT, "composite_tests.csv")); fwrite(nul, file.path(OUT, "composite_null.csv.gz"))
fwrite(rna[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 4) else z)], file.path(OUT, "rna_nucleus_vs_fragment.csv"))
fwrite(amb, file.path(OUT, "ambiguity.csv")); fwrite(ct, file.path(OUT, "seed_composites.csv.gz"))
fwrite(rbindlist(assign_rows), file.path(OUT, "polygon_assignments.csv.gz"))
cat("\n== ambiguity\n"); print(amb)

## ---- figures -----------------------------------------------------------------------
## (a) seed composite size, observed vs null, both samples, three thresholds
sz <- res[stat %in% c("frac_multi", "frac_ge300", "med_members")]
sz[, stat_lab := factor(c(frac_multi = "share with >= 2 members", frac_ge300 = "share >= 300 um^2", med_members = "median members")[stat],
                        levels = c("share with >= 2 members", "share >= 300 um^2", "median members"))]
nl <- nul[, .(samp, thr, partition, frac_multi, frac_ge300, med_members)]
nl <- melt(nl, id.vars = c("samp", "thr", "partition"), variable.name = "stat", value.name = "v")
nl[, stat_lab := factor(c(frac_multi = "share with >= 2 members", frac_ge300 = "share >= 300 um^2", med_members = "median members")[as.character(stat)],
                        levels = levels(sz$stat_lab))]
nl[, x := factor(paste0(samp, " >= ", thr), levels = unique(paste0(rep(names(SEED), each = 3), " >= ", THR)))]
sz[, x := factor(paste0(samp, " >= ", thr), levels = levels(nl$x))]
f1 <- ggplot() +
  geom_violin(data = nl, aes(x, v), fill = "#e9e7e1", colour = NA, scale = "width", width = .7) +
  geom_point(data = sz, aes(x, obs, colour = samp, shape = partition), size = 2.6, position = position_dodge(width = .5)) +
  geom_point(data = sz, aes(x, obs_noseam, colour = samp, shape = partition), size = 1.4, alpha = .5, position = position_dodge(width = .5)) +
  facet_grid(stat_lab ~ partition, scales = "free_y") +
  scale_colour_manual(values = SAMP_COLS, name = NULL) + scale_shape_manual(values = c(components = 16, anchored = 17), name = NULL) +
  labs(title = "Are chemokine-bright polygons contiguous beyond chance?",
       subtitle = sprintf("Seed composites (IB1 seeds >= 50, RB1 >= 20) under two partitions; grey = %d within-FOV permutations of the expression rows. Small faint point = seam-free composites only", N_PERM),
       x = NULL, y = NULL, caption = "Source: output/35_giant_cells/composite_tests.csv, composite_null.csv.gz.") +
  theme_viz(10) + theme(strip.text = element_text(colour = INK2, face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1))
sf("fig35_contiguity.png", f1, 11, 8)

## (b) nuclei per member in multi-member composites, components partition
nu <- res[partition == "components" & stat %in% c("nuc_per_member_multi", "frac_multi_le1nuc", "frac_chain")]
nu[, stat_lab := c(nuc_per_member_multi = "dense nuclei per member (>= 2-member composites)", frac_multi_le1nuc = "share of >= 2-member composites with <= 1 nucleus",
                   frac_chain = "share of composites with >= 2 nuclei (chains)")[stat]]
nn <- melt(nul[partition == "components", .(samp, thr, nuc_per_member_multi, frac_multi_le1nuc, frac_chain)], id.vars = c("samp", "thr"), variable.name = "stat", value.name = "v")
nn[, stat_lab := c(nuc_per_member_multi = "dense nuclei per member (>= 2-member composites)", frac_multi_le1nuc = "share of >= 2-member composites with <= 1 nucleus",
                   frac_chain = "share of composites with >= 2 nuclei (chains)")[as.character(stat)]]
nn[, x := factor(paste0(samp, " >= ", thr), levels = levels(nl$x))]; nu[, x := factor(paste0(samp, " >= ", thr), levels = levels(nl$x))]
f2 <- ggplot() +
  geom_violin(data = nn, aes(x, v), fill = "#e9e7e1", colour = NA, scale = "width", width = .7) +
  geom_point(data = nu, aes(x, obs, colour = samp), size = 2.6) +
  facet_wrap(~ stat_lab, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = SAMP_COLS, name = NULL) +
  labs(title = "One nucleus and several fragments, or one nucleus per member?",
       subtitle = "Connected-component composites containing a seed; grey = permutation null (same geometry and nuclei, expression rows shuffled within FOV)",
       x = NULL, y = NULL, caption = "A giant cell carved into polygons has few nuclei for its members; a cluster of activated cells has about one per member. Source: composite_tests.csv.") +
  theme_viz(10) + theme(strip.text = element_text(colour = INK2, face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1))
sf("fig35_nuclei_per_member.png", f2, 13, 4.8)

## (c) RNA: fragment over nucleus member, by gene, three ambiguity variants
if (nrow(rna) == 0) { cat("no RNA contrast (too few multi-member composites); figure skipped\n"); quit(save = "no") }
rr <- rna[variant == "best" & thr == SEED[samp]]
rr[, gene := factor(gene, levels = rev(TEST))]
GC <- c(HRS = "#0b3d7a", ambiguous = "#7a5aa6", DC = "#2a9d8f", macro = "#9a968d", Bprog = "#9a968d", T = "#d03b3b")
f3 <- ggplot(rr, aes(ratio_frag_over_nuc, gene, colour = set)) +
  geom_vline(xintercept = 1, colour = AXIS, linewidth = .4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = .5, orientation = "y") + geom_point(size = 1.8) +
  facet_wrap(~ samp, nrow = 1) +
  scale_x_continuous(trans = "log2", breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = function(x) paste0(x, "x")) +
  scale_colour_manual(values = GC, name = NULL) +
  labs(title = "Within a composite: do the nucleus-free fragments carry the programme the nucleus member lacks?",
       subtitle = "Ratio of depth-matched gene excess, fragments / nucleus member, anchored partition, >= 2-member seed composites; 95% FOV cluster bootstrap",
       x = "Fragment / nucleus excess ratio (log scale)", y = NULL,
       caption = "Ambiguous fragments assigned to the best parent here; the split and dropped variants are in rna_nucleus_vs_fragment.csv. T genes = leakage gauge.") +
  theme_viz(9, xgrid = TRUE)
sf("fig35_rna_fragment_vs_nucleus.png", f3, 11, 7)
cat("\ndone:", OUT, "\n")
