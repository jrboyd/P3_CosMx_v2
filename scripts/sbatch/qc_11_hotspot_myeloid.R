## Why are macrophages DEPLETED in the CCL17/CCL22 hotspots?
##
## THE ANOMALY. Step 11 found hotspots enriched for T cells (11.7% vs 6.4%) but
## DEPLETED for myeloid content: CD68 0.82x, CD163 0.60x, lineage-myeloid 10.1% vs
## 15.7%. Tumour-associated macrophages are a defining cHL feature, so this contradicts
## the textbook picture and has been carried as unexplained since it was measured.
##
## THE HYPOTHESIS THIS TESTS. CCL17/TARC and CCL22/MDC are not only HRS products --
## they are the signature secretions of MATURE / MIGRATORY DENDRITIC CELLS
## (LAMP3+ CCR7+). If the hotspots are DC-rich rather than macrophage-rich, then:
##   - LAMP3/CCR7 and the cDC1 axis (CLEC9A/XCR1/BATF3/THBD) should be UP
##   - the macrophage axis (CD163/MRC1/MSR1/MARCO/SIGLEC1) should be DOWN
##   - and step-01's `lineage` caller would score those cells as something other than
##     myeloid, because its myeloid set (CD68/CD163/LYZ/ITGAX/C1QA/C1QB/CD14) is
##     macrophage-weighted and contains no DC-specific gene at all.
## That would make "macrophage depletion" a DC-enrichment result, not a paradox.
##
## THREE READOUTS, because they can disagree and the disagreement is informative:
##   A  bin-level DEPTH-WEIGHTED share   sum(marker counts) / sum(nCount)  [step-11's]
##   B  bin-level MEAN OF PER-CELL share mean(marker cp10k per cell)       [unweighted]
##   C  fraction of cells by CALLED lineage
##   D  within called-myeloid cells only: is the per-cell marker level itself lower,
##      or are there simply fewer of them?
## A is dominated by deep cells; B weights every cell equally; C ignores expression
## level entirely. A real compositional shift moves all three the same way.
##
##   sbatch --mem=96G --cpus-per-task=2 --time=2:00:00 \
##     scripts/sbatch/run_rscript.sh scripts/sbatch/qc_11_hotspot_myeloid.R
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); setDTthreads(1)
})

IN_RDS <- "output/04_merge/merged_processed.Rds"
BIN_F  <- "output/11_hodgkin_niche/niche_bins_all_scales.csv.gz"
OUT    <- "output/11_hodgkin_niche"; FIG <- file.path(OUT, "figures")
SAMPLE <- "IB1"; DEFECT_CL <- "i"; BIN_UM <- 100; NICHE_Q <- 0.99

SETS <- list(
  chemokine  = c("CCL17","CCL22"),
  DC_mature  = c("LAMP3","CCR7"),                       # the CCL17/CCL22 source
  DC_cDC1    = c("CLEC9A","XCR1","BATF3","THBD"),
  DC_cDC2    = c("CD1C","FCER1A"),
  pDC        = c("LILRA4","IL3RA"),
  macrophage = c("CD163","MRC1","MSR1","MARCO","SIGLEC1"),
  mac_core   = c("CD68","C1QA","C1QB","CD14"),
  mono_myelo = c("LYZ","ITGAX","FLT3","IRF8"),
  Tcell      = c("CD3E","CD4"))
GENES <- unique(unlist(SETS, use.names = FALSE))

SURF<-"#fcfcfb"; INK<-"#0b0b0b"; INK2<-"#52514e"; MUTED<-"#898781"; GRID<-"#e1e0d9"; AXIS<-"#c3c2b7"
GCOLS <- c(hotspot = "#2a78d6", background = "#a9a7a1")
theme_viz <- function(base=11) theme_minimal(base_size=base) +
  theme(plot.background=element_rect(fill=SURF,colour=NA),
        panel.background=element_rect(fill=SURF,colour=NA),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour=GRID,linewidth=.25),
        plot.title=element_text(colour=INK,face="bold",size=base+2),
        plot.subtitle=element_text(colour=INK2,size=base-1),
        plot.caption=element_text(colour=MUTED,size=base-3,hjust=0),
        legend.position="top", legend.justification="left")
sf <- function(n,p,w,h){ggsave(file.path(FIG,n),p,width=w,height=h,dpi=150,bg=SURF)
                       cat("wrote",file.path(FIG,n),"\n")}
sc <- function(n,d){fwrite(d,file.path(OUT,n));cat("wrote",file.path(OUT,n),"\n")}

cat("============================================================\n")
cat("step 11 QC -- macrophage depletion in the hotspots: real, or DCs?\n")
cat("============================================================\n")

obj <- readRDS(IN_RDS); rn <- rownames(obj@meta.data)
md <- as.data.table(obj@meta.data[, c("merge_sample","x_slide_mm","y_slide_mm",
                                      "nCount_RNA","inSitu_clusters","lineage"),
                                  drop = FALSE])
md[, cell_id := rn]
cnt <- LayerData(obj[["RNA"]], layer = "counts")
G <- intersect(GENES, rownames(cnt))
if (length(G) < length(GENES))
  message("absent and dropped: ", paste(setdiff(GENES, G), collapse = ", "))
M <- as(cnt[G, , drop = FALSE], "dgCMatrix"); colnames(M) <- rn
rm(obj, cnt); invisible(gc())
SETS <- lapply(SETS, intersect, G); SETS <- SETS[lengths(SETS) > 0]
## Persist the marker sets AS SCORED, after the panel intersection above has had its say.
## The results CSV carries only set names, so without this the exact gene membership
## behind every ratio lives solely in this source file and cannot be checked against a
## finished run.
sc("hotspot_myeloid_marker_genes.csv",
   rbindlist(lapply(names(SETS), function(s)
     data.table(set = s, gene = SETS[[s]], n_in_set = length(SETS[[s]]),
                dropped_absent_from_panel = FALSE, sample = SAMPLE))))

md <- md[merge_sample == SAMPLE & nCount_RNA > 0 & inSitu_clusters != DEFECT_CL]
md <- md[!is.na(x_slide_mm) & !is.na(y_slide_mm)]
gi <- match(md$cell_id, colnames(M)); M <- M[, gi, drop = FALSE]
for (s in names(SETS))
  set(md, j = paste0("ct_", s), value = as.vector(Matrix::colSums(M[SETS[[s]], , drop=FALSE])))
for (s in names(SETS))
  set(md, j = paste0("cp_", s), value = 1e4 * md[[paste0("ct_", s)]] / md$nCount_RNA)
md[, `:=`(bx = floor(x_slide_mm*1000/BIN_UM), by_ = floor(y_slide_mm*1000/BIN_UM))]

## Join the bin CALLS from step 11 rather than recomputing them, so this is testing the
## same hotspots that were reported and not a lookalike set.
B <- fread(BIN_F)[bin_um == BIN_UM & samp == SAMPLE]
B[, grp := fifelse(chemo_cp10k >= quantile(chemo_cp10k, NICHE_Q), "hotspot", "background")]
md <- merge(md, B[, .(bx, by_, grp, seam_flag)], by = c("bx","by_"))
cat(sprintf("\ncells inside called bins: %d (%s)\n", nrow(md),
            paste(md[, .N, by=grp][, paste0(grp,"=",N)], collapse=" ")))

## ---- A / B: two ways of averaging the same markers ---------------------------
cpn <- paste0("cp_", names(SETS)); ctn <- paste0("ct_", names(SETS))
A <- md[, c(.(readout = "A depth-weighted share"),
            lapply(setNames(names(SETS), cpn), function(s)
              1e4 * sum(get(paste0("ct_", s))) / sum(nCount_RNA))), by = grp]
B2 <- md[, c(.(readout = "B mean of per-cell share"),
             lapply(.SD, mean)), by = grp, .SDcols = cpn]
AB <- rbind(A, B2, fill = TRUE)
long <- melt(AB, id.vars = c("grp","readout"), variable.name = "set", value.name = "v")
long[, set := sub("^cp_", "", set)]
w <- dcast(long, readout + set ~ grp, value.var = "v")
w[, ratio := round(hotspot / pmax(background, 1e-9), 2)]
w[, `:=`(hotspot = round(hotspot, 3), background = round(background, 3))]
cat("\n=== A and B: marker sets, hotspot vs background\n")
cat("    DC hypothesis -> DC_mature and DC_cDC1 UP, macrophage DOWN, in BOTH readouts\n")
print(dcast(w, set ~ readout, value.var = "ratio")[order(-`A depth-weighted share`)])
sc("hotspot_myeloid_marker_sets.csv", w)

## ---- C: called cell types ------------------------------------------------------
lin <- md[, .N, by = .(grp, lineage)][, pct := 100*N/sum(N), by = grp]
lw <- dcast(lin, lineage ~ grp, value.var = "pct")
lw[, ratio := round(hotspot / pmax(background, 1e-9), 2)]
cat("\n=== C: fraction of cells by CALLED lineage (%)\n")
cat("    step-01's myeloid set is CD68/CD163/LYZ/ITGAX/C1QA/C1QB/CD14 -- macrophage-\n")
cat("    weighted, with no DC-specific gene, so DCs cannot be called as such here\n")
print(lw[order(-hotspot)]); sc("hotspot_myeloid_lineage.csv", lin)

## ---- D: within called-myeloid cells only ---------------------------------------
## Separates "fewer macrophages" from "the same macrophages expressing less".
myl <- md[lineage == "myeloid"]
if (nrow(myl) > 50) {
  D <- myl[, c(.(cells = .N), lapply(.SD, mean)), by = grp, .SDcols = cpn]
  dl <- melt(D, id.vars = c("grp","cells"), variable.name = "set", value.name = "v")
  dl[, set := sub("^cp_", "", set)]
  dw <- dcast(dl, set ~ grp, value.var = "v")
  dw[, ratio := round(hotspot / pmax(background, 1e-9), 2)]
  dw[, `:=`(hotspot = round(hotspot,3), background = round(background,3))]
  cat(sprintf("\n=== D: WITHIN called-myeloid cells (%d hotspot, %d background)\n",
              D[grp=="hotspot", cells], D[grp=="background", cells]))
  cat("    ratio ~1 -> the macrophages present are normal, there are just fewer\n")
  cat("    ratio <1 -> the macrophages themselves are different in the hotspots\n")
  print(dw[order(-hotspot)]); sc("hotspot_myeloid_within_myeloid.csv", dw)
} else message("too few myeloid cells for readout D")

## ---- verdict --------------------------------------------------------------------
g <- function(st, ro) {
  v <- w[set == st & readout == ro, ratio]; if (length(v) == 1) v else NA_real_
}
dcm_A <- g("DC_mature","A depth-weighted share"); dcm_B <- g("DC_mature","B mean of per-cell share")
mac_A <- g("macrophage","A depth-weighted share"); mac_B <- g("macrophage","B mean of per-cell share")
dc1_A <- g("DC_cDC1","A depth-weighted share")
cat("\n============================================================\n")
cat(sprintf("DC_mature (LAMP3/CCR7)  A %.2fx  B %.2fx\n", dcm_A, dcm_B))
cat(sprintf("DC_cDC1                 A %.2fx\n", dc1_A))
cat(sprintf("macrophage              A %.2fx  B %.2fx\n", mac_A, mac_B))
if (is.finite(dcm_A) && is.finite(mac_A) && dcm_A > 1.3 && mac_A < 0.85)
  cat("VERDICT: DC-mature markers UP and macrophage markers DOWN in the same regions.\n",
      "         The 'macrophage depletion' is DC enrichment -- the hotspots look like\n",
      "         mature/migratory DC territory, which is exactly what secretes CCL17 and\n",
      "         CCL22. The step-01 lineage caller has no DC set, so those cells were\n",
      "         never going to be counted as myeloid.\n", sep = "")
if (is.finite(dcm_A) && dcm_A <= 1.3)
  cat("VERDICT: DC markers are NOT elevated. The DC hypothesis fails and the macrophage\n",
      "         depletion stands as unexplained -- do not paper over it.\n", sep = "")

pl <- w[readout == "A depth-weighted share"]
p1 <- ggplot(pl, aes(reorder(set, ratio), ratio, fill = ratio > 1)) +
  geom_col(width = .7) + coord_flip() +
  geom_hline(yintercept = 1, colour = INK, linetype = "22", linewidth = .4) +
  scale_fill_manual(values = c(`TRUE` = GCOLS[["hotspot"]], `FALSE` = "#eb6834"),
                    guide = "none") +
  labs(title = "Are the hotspots macrophage-poor, or dendritic-cell-rich?",
       subtitle = paste0("Marker sets in hotspot bins relative to IB1 background. ",
                         "CCL17/CCL22 are secreted by mature\n(LAMP3+CCR7+) DCs as well ",
                         "as by HRS cells, so DC enrichment would explain the apparent\n",
                         "macrophage loss without contradicting cHL biology."),
       x = NULL, y = "hotspot : background ratio",
       caption = sprintf("Data: %s/hotspot_myeloid_marker_sets.csv", OUT)) + theme_viz()
sf("fig_11_hotspot_myeloid_vs_dc.png", p1, 8.5, 5.5)
cat("============================================================\nOK\n")
