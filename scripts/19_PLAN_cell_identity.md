# Step 19 — cell identity, done properly

## Context

Cell identity has to be solid for this project, and right now it is not. Four different
names for a cell are in circulation (`inSitu_clusters`, `pearson_clusters_low/high`, the
step-09b audited labels, the step-16 argmax `celltype`), they disagree, and the project has
no single place that says which one a given claim rests on.

The design here keeps **two** concepts of cell type deliberately, because they answer
different questions and fail in different ways:

- **`celltype_argmax`** — direct per-cell marker evidence. Independent of clustering and of
  QC tier. Silent on most cells.
- **`celltype_cluster`** — guilt by association. A cell with no marker transcripts meeting
  the argmax cutoffs still co-clusters with cells that do, and inherits their label.

Both are needed, and the reason is measurable: **43.4% of the high-confidence set gets no
argmax call at all** (`ambiguous` 20.8%, `all_markers_0` 15.2%, `below_background` 7.4%).
Marker evidence alone leaves nearly half the tissue unnamed. That is the gap clustering
fills, and it is why neither concept can be dropped.

**The clustering is gated.** No cell-type call is made until the Seurat clusters come back
clean — no large cluster of low-signal or defect-residue tissue. Dropping `is_defect` is
necessary but may not be sufficient (see 19b), and a surviving low-signal cluster is a
finding about the filter, not something to label around.

Two tiers of data, per the requirement that statements be clean but maps be dense:

| | cells | use |
|---|---|---|
| **high-confidence** | 1,993,349 (IB1 725,076 / RB1 1,268,273) | every cell-type **statement** |
| **rescued / low-confidence** | 2,413,591 | spatial visualisation only |

(The high-confidence figure is 2,053,618 on the symmetric tier with the defect dropped, less
the 60,269 cells removed by the per-cell low-information floor added in 19a — see 19b.)

---

## What step 18 established, and what it forbids

Step 18 measured transcript leakage between adjacent cells. Four results constrain this step.

1. **Both callers are neighbourhood-contaminated, and the clustering is worse.** Holding a
   cell's own PC-marker count at zero, the *cluster* calls it PC 30.1% of the time within
   10 µm of a plasma cell and 3.6% beyond 80 µm — an 8.4× gradient on proximity alone.
   Across the same range the cluster's PC rate moves 6.8× against argmax's 2.5×. The likely
   mechanism is specific: the leaked mass is light chains, which argmax's PC set
   (MZB1/XBP1/PRDM1/SDC1/DERL3/TNFRSF17) excludes and whole-transcriptome clustering
   includes. **Argmax is partly protected because its marker set is narrow.**

2. **Leakage is real and directional.** Marker transcripts are displaced toward the
   marker-rich neighbour, more so where argmax and the cluster disagree: myeloid +0.127 log2
   [0.070, 0.189], T +0.114, fibro +0.103. These are **floors** — cells where leakage
   captured both callers sort into the "agrees" column and read as confirmations.

3. **Compartment restriction does not fix it. Do not build it in.** Nuclear-only argmax
   scored −3.53 pp against a marker-matched random control. `CellComp` is assigned relative
   to the segmentation mask, so a bled-in transcript landing over the wrong cell's nucleus
   is labelled `Nuclear`; compartment tracks position inside a possibly-wrong boundary, not
   provenance. Markers are also the cytoplasm-enriched fraction (64.0% nuclear vs 78.4% for
   transcripts overall), so nuclei are the marker-poor compartment.

4. **Never transfer labels by spatial proximity.** Doing so manufactures exactly the halo
   artefact step 18 measured. The rescue in 19f is expression-space only.

---

## Can leakage be corrected at this stage? Only at the cluster level.

**Leakage is real at the population level and unmeasurable at the single-cell level.** This
was tested directly, and it settles the design.

As per-cell classifiers of argmax/cluster disagreement, every candidate score is at chance:

| per-cell score | AUC range across lineages |
|---|---|
| positional displacement `z_A` (the *validated* leakage signal) | 0.446 – 0.526 |
| donor marker load `donorA_k` | 0.354 – 0.605 |
| nucleus-vs-cytoplasm discrepancy | 0.383 – 0.554 |

The myeloid leakage effect is +0.127 log2 measured across 5,913 cells; per cell it is buried
under shot noise. **No per-cell leakage flag is going to work — not compartment, not
position, not donor load.** Do not build one, and do not let a future step try again.

### Nucleus-vs-cytoplasm discrepancy specifically: tested, and it measures sparsity

The intuition is that leaked transcripts should be peripheral, so a cell whose lineage
markers are anomalously non-nuclear is a leakage victim. It does not hold.

Mean compartment discrepancy by positional-`z_A` decile is a **symmetric U** — myeloid
0.329 at z = −2.95, 0.104 at z = −0.07, 0.343 at z = +3.15; fibro 0.293 / 0.070 / 0.284.
Cells with markers pushed out of the nucleus are displaced extremely in *either* direction,
equally. Leakage would give a monotone rise on the positive side only. Both quantities blow
up when a cell has few, scattered marker transcripts: this is a variance signature.

Where the linear slope of `z_A` on discrepancy is significant it runs the **wrong way** —
B −0.237 (p = 1e-6), T −0.178 (p = 0.0004). Markers are also 63.0% nuclear against 74.9%
for the same cells' other transcripts, so non-nuclear markers are the *norm*, not an anomaly.
(That is the same fact as the 64.0%/78.4% figure above, measured on the narrower subset this
test runs on — cells with an adjacent donor and at least two markers — not a second result.)

The mechanism explains it: `CellComp` is assigned relative to the segmentation mask, so a
transcript that bled into the wrong cell and landed over *that* cell's nucleus is called
`Nuclear`. Compartment tracks position inside a possibly-wrong boundary, never provenance.

*(Untested variant, if anyone wants to reopen it: a gene-matched discrepancy, comparing each
marker gene's nuclear fraction in a cell against that gene's global nuclear fraction. Given
the U-shape is a variance signature and the direction is wrong where significant, it is
unlikely to rescue the idea.)*

### What the step does instead

1. **Correct at the cluster level, where the signal lives.** The PC contamination gradient is
   **8.4×** — enormous, and only visible by aggregating. That is 19e, and it is the only
   mechanism in this plan that can work.
2. **Carry donor load and donor distance as per-cell COVARIATES, never as flags.** The honest
   use of a weak signal: downstream analyses condition on it rather than filter by it.
   Available for every cell from `18_script_select_cells.R`.
3. **Prefer the narrow-marker call where a promiscuous neighbour exists.** Argmax is 2.7×
   *less* neighbourhood-contaminated than whole-transcriptome clustering for PC, precisely
   because its marker set excludes light chains. For plasma cells specifically,
   `celltype_argmax` should be the primary call and the cluster label the secondary one.

---

## Why this reverses the merged-label ban

`CLAUDE.md` fact 1 says *"Do not use the merged labels for cell identity."* The documented
reason (`09c_script_cluster_consistency.R`) is dilution, not merging as such: merged
`pearson_clusters_low` cluster 1 was 93.6% IB1 cells that IB1's own audit calls
`defect_dominated` at 70.1% defect. Pooling diluted the defect fraction to 42%, under the
`DEFECT_MAX = 0.50` guard, so it escaped the guard and acquired the label "Neutrophil".

**Dropping defect cells before clustering removes the mechanism.** With no defect cells in
the object there is nothing left to dilute. That is the whole justification for the design,
and it is why the clustering must be re-run rather than subsetted from step 04.

When this step lands, `CLAUDE.md` fact 1 must be rewritten — not deleted — to say that
merged labels are valid *only* on a defect-free object and why.

## The other trap, which this step must not walk into

`09b_script_marker_relabel.R:199` scores clusters by z **across clusters within a sample**
(`Z <- t(scale(t(log1p(CP))))`). That collapses when one lineage dominates: RB1 is ~60% B
cells and its largest cluster (836,909 cells) audits as `low-information` at only 3% defect,
even though its own table records `best_supported = B cell` on BANK1/FCER2/CD19/CD79A/CD79B
with a 0.75 margin (documented at `16_script_umap_panels.R:23-30`).

This **will** recur here. IB1 loses 59% of its cells to the defect while RB1 loses 8%, so the
high-confidence object is **62.6% RB1**, and B is **50.5% of all argmax-resolved cells in
it**. A relative, across-cluster score is the wrong instrument for this object.

→ **19d labels clusters by absolute argmax composition, not by across-cluster z.**

---

## Design

New step 19. Logically sits after step 04 (merge and process); numbered 19 because
renumbering would break every existing path. Outputs to `output/19_celltypes/`.

Inputs, all existing: `output/04_merge/merged_processed.Rds`,
`output/16_regions/cells_meta_all.csv.gz`, `output/16_celltypes/cell_celltypes.csv.gz`,
`output/18_leakage/{IB1,RB1}_centroids.csv.gz`, `scripts/lib/lineage.R`,
`scripts/lib/cluster_labels.R`.

### 19a — `19_script_hc_object.R`: build and cluster the high-confidence object

Subset merged to `tier_matched == TRUE & is_defect == FALSE` → 2,053,618 cells.

`tier_matched` is the symmetric tier (`qc_spatial_fovpx == "ok" & qcCellsPassed`). The
alternative tiers inherit the staining rule that cost IB1 2,013,225 cells and was never
evaluated for RB1, which would contaminate every IB1-vs-RB1 statement (`CLAUDE.md` fact 2).

Re-run step 03's block verbatim — `gene_frequency` → `FindVariableFeatures(nfeatures = 900)`
→ `sparse_quasipoisson_pca_seurat` → `make_umap` → `FindClusters` at resolution 0.4 and 0.8
(`03_script_filter_and_process.R:244-270`). `16_script_umap_panels.R` already replicates
this block for the recomputed-embedding case and is the code to copy; the one difference is
that **this step DOES re-run `FindClusters`**, which that script deliberately does not.

Name the results `hc_clusters_low` / `hc_clusters_high` so they can never be mistaken for
`pearson_clusters_*` from step 03, whose ids mean something else.

Feasibility is already proven: step 03 ran this pipeline on RB1's 1.87M cells and step 04 on
all 4.4M. Budget `--mem=200G --cpus-per-task=4 --time=12:00:00`.

### 19b — `19_script_cluster_health.R`: THE GATE. Nothing proceeds until clusters come back clean

**No cell-type call is made until the clustering returns clean — no large cluster of
low-signal tissue surviving.** This is a stop-and-look checkpoint, not a report at the end.

Dropping `is_defect` removes InSituType cluster `i`, but **low-information tissue is not the
same set as cluster `i`**. IB1's `pearson_clusters_low` cluster 0 audits as
`defect_dominated` at 70.1% defect — which means ~30% of that mass is non-defect cells that
**survive this filter** and can re-form a low-signal cluster in the new object. If they do,
labelling it would repeat the 09c laundering failure by a different route.

Per-cluster diagnostics, all from evidence rather than from cluster membership so the check
is not circular:

| metric | source | why |
|---|---|---|
| `frac_resolved` | share of members argmax calls a real lineage | object-wide baseline is **56.5%** (1,161,000 / 2,053,618) |
| median `nCount_RNA`, `nFeature_RNA` | `cells_meta_all.csv.gz` | depth |
| `propNegative` / negprobe ratio | vendor metadata | the actual signal-to-noise measure — depth alone predicts detection at AUC ~0.70 for *any* sparse gene, so negprobe rate is the honest one |
| `complexity` | vendor metadata | library diversity |
| marker share (lineage markers as % of cluster transcripts) | `panel_counts_all.Rds` | 09b's `MIN_MK_SHARE`, reused |
| sample skew | `merge_sample` | the 09c cluster-1 signature |

**Pass criteria.** No cluster larger than 5% of the object may have `frac_resolved` below
half the object baseline (< ~28%) together with depressed depth/complexity and an elevated
negprobe rate. No cluster may draw > 90% of its cells from one sample while that sample's
own per-sample audit declines it.

**If the gate fails**, do not label and do not loosen the criteria. Work the remedy ladder:

1. Add **per-cell** low-information criteria — negprobe ratio, `nFeature_RNA` floor,
   `complexity` — defined on the cell's own evidence, never on its cluster.
2. Re-run 19a and re-test.
3. **Do not** reach for the FOV-level `qc_spatial_i_fraction` flag: `CLAUDE.md` records that
   it is far too aggressive and cuts IB1 to 196,316 cells. The defect is dropped per cell.

Report the gate result before proceeding. A surviving low-signal cluster is a finding about
the filter, not a nuisance to be labelled around.

### 19c — per-cell argmax: **reuse, do not recompute**

Take `output/16_celltypes/cell_celltypes.csv.gz` as it stands. It is tier-independent by
design — called once over all 4,406,940 pre-QC cells so the fold-over-mean gene scaling
comes from whole tissue rather than from whichever subset is in play.

Step 18 measured the cost of recomputing on a subset: **12.7% of calls change** from the
scaling basis shifting alone. Re-running it on the high-confidence cells would silently
produce a third, incompatible vocabulary.

Extend only by adding stage-2 sets that do not exist yet (pDC, GC B). Keep the existing
discipline: subsets are called only within a resolved parent, and every checkpoint gene
stays held out so "Tfh express PD-1" cannot become a restatement of the definition.

### 19d — `19_script_cluster_labels.R`: label clusters by argmax composition

The composition-robust replacement for 09b's z-scoring, and the literal implementation of
guilt-by-association.

For each `hc_clusters_low` cluster, over its members that argmax **resolved**:

- `label` = modal argmax family
- `purity` = share of resolved cells taking the modal call
- `margin` = purity − runner-up
- `n_resolved`, `frac_resolved` — a cluster where almost nothing resolves gets no label

Absolute, so it cannot collapse when one lineage dominates the object.

Carry over 09b's guards, which are sound and independent of the scoring method:
`MIN_FAM_GENES = 2` (so **Neutrophil, which has one gene on this panel — CSF3R; FCGR3B is
absent — can never reach `confident`**), the marker-share floor, and the `PARENT` fallback at
`09b:85` that resolves to parent lineage when a subset call is not supported.

### 19e — `19_script_halo_qc.R`: the leakage QC, new and only possible after step 18

For every cluster, run step 18's contamination diagnostic: P(cell belongs to cluster *k*),
holding the cell's **own** marker count for *k*'s family fixed, as a function of the nearest
marker-rich neighbour's load. A cluster whose membership is driven by the neighbour rather
than by the cell's own evidence is a **halo cluster**, not a cell type.

Reuse the donor machinery from `18_script_select_cells.R` — nearest cell carrying
≥ `MIN_DONOR` markers of family F, searched over all non-defect cells, geometry in canvas px
from `*_centroids.csv.gz` so the fitted mm→px affine never enters.

Emit a per-cluster `halo_suspect` flag. Expect it to fire on plasma-cell halos. Any
downstream statement must respect it.

### 19f — `19_script_rescue.R`: label the other 2,353,322 cells, in two columns

Never merged into one column, so no figure can silently mix evidence types.

- **`celltype_argmax`** — already available for every one of the 4.4M cells from 19c.
  Unaffected by clustering, QC tier, or defect status.
- **`celltype_projected`** — project low-confidence cells onto the 19a PCA and kNN-transfer
  the `hc_clusters_low` label. **Expression space only.**
- **`confidence`** — `high` (in the HC set) / `rescued` / `unresolved`.
- **`donorA_k`, `donorA_dist`** — the leakage covariates, carried per cell from
  `18_script_select_cells.R`. Covariates to condition on, **not** filters. Per-cell leakage
  is unmeasurable (AUC ~0.5, see above); these exist so a downstream analysis can adjust for
  neighbourhood load, not so anyone can drop cells with them.

Gate the transfer the way the HCA tonsil transfer failed (disabled 2026-08-06 because its
labels anti-track marker evidence). Hold out a random 10% of **high-confidence** cells,
transfer to them, and require the transferred label to agree with (a) their own `hc_clusters_low`
label and (b) `celltype_argmax` where argmax resolves. **If it anti-tracks marker evidence,
`celltype_projected` ships as NA rather than as a label.** Pre-commit that stopping rule
before looking at the result.

### 19g — `19_summary_report.Rmd`: characterisation

- Composition per sample per label — **high-confidence cells only**.
- Marker dotplot per label, reported twice, with and without the defining genes (only the
  second is evidence).
- UMAP panels: `hc_clusters_low`, `celltype_argmax`, and their agreement.
- Spatial maps over all cells with `confidence` encoded, so the reader sees which parts of
  the map are carried by rescued labels.
- Self-contained HTML in `output/19_celltypes/`, `scripts/lib/methods_notes.Rmd` as a knitr
  **child**, `mathjax: null`, plain markup for formulae.

---

## Rules the step must state in its own outputs

- Cell-type **statements** use `confidence == "high"` only.
- Spatial figures may use every cell, but must show the confidence split.
- Region is the unit of every p-value (154 neighbourhoods are not 400k observations).
- Neutrophil rests on one gene. Always. Every call, permanently caveated.
- `halo_suspect` clusters are not cell types.
- **No per-cell leakage filtering.** Every candidate score is at chance per cell (AUC
  0.35–0.61). Leakage is corrected at the cluster level or not at all.
- For plasma cells, `celltype_argmax` is primary and the cluster label secondary — argmax's
  marker set excludes light chains and is 2.7× less contaminated there.

## Verification

1. **The cluster-health gate (19b) passes.** This is the precondition for every other
   check; nothing downstream is meaningful until the clustering comes back with no large
   low-signal cluster.
2. **Regression against the one independent audit.** Compare `hc_clusters_low` labels to
   step 09b's per-sample IB1 labels — IB1 is the only sample with HierAnnot to have been
   audited against. Expect agreement on the `confident` clusters; investigate every
   disagreement rather than averaging over it.
3. **The laundering check, re-run.** Confirm no HC cluster exceeds `DEFECT_MAX` (should be
   0% by construction) and that no cluster reproduces the merged cluster-1 failure — one
   sample contributing >90% of a cluster that its own per-sample audit declines.
4. **The domination check.** Confirm the argmax-composition labelling assigns **B cell** to
   RB1's dominant cluster, where 09b's z-scoring returned `low-information`. This is the
   specific failure 19d exists to fix; if it recurs, the design is wrong.
5. **The transfer gate** (19f), with its pre-committed stopping rule.
6. **Arithmetic.** Cells by confidence class must sum to 4,406,940.
7. **Sample provenance.** Join on `merge_sample`, never on the cell-id prefix — `c_2_` is
   IB1 and `c_1_` is RB1, and guessing it the obvious way round swaps the samples in every
   row without erroring.

## Results as built (19a-19e complete)

**High-confidence object**: 1,993,349 cells (IB1 725,076 / RB1 1,268,273) after the symmetric
tier, defect drop, and a per-cell low-information floor (`nCount >= 50 & nFeature >= 30`,
60,269 cells / 2.93%) applied BEFORE clustering so the junk never enters the .Rds. Gate
passes at both resolutions.

**Vocabulary -- 100% of high-confidence cells labelled:**

| label | cells | % | how |
|---|---|---|---|
| B cell | 1,076,141 | 54.0 | argmax composition, 15 clusters merged |
| Myeloid | 269,628 | 13.5 | argmax composition |
| Stromal/FDC | 229,006 | 11.5 | argmax + the FDC rule |
| T cell | 142,361 | 7.1 | argmax composition |
| **Hodgkin chemokine niche** | 82,228 | 4.1 | module rule; **`is_niche`, not a cell type** |
| Endothelial | 67,246 | 3.4 | argmax composition |
| **pDC** | 50,829 | 2.5 | module rule |
| NK cell | 48,953 | 2.5 | enrichment route |
| **Plasma cell** | 26,957 | 1.4 | enrichment route; **`halo_suspect` at 12.7x** |

### Six things learned in the build that the design above did not anticipate

1. **The two largest unnamed clusters were real populations, not QC failures.** Cluster 9
   (80.5% IB1) is CCL17 36.7x / CCL22 30.4x / BATF3 15.3x -- the Hodgkin compartment.
   Cluster 15 is CLEC4C 11.4x / LILRA4 7.6x / IL3RA 4.8x / IRF7 11.0x -- pDC, which argmax
   split evenly between the PC and B sets (purity 0.227 vs 0.225) because pDCs express MZB1
   and IRF4, the same confusion CLAUDE.md records for HierAnnot. Neither has a family in
   `LINEAGE_SETS_COARSE`, so both are named by cluster-level MODULE rules, applied only
   where argmax composition left the cluster unnamed. `lib/lineage.R` is not modified.
2. **Purity is the wrong yardstick for a minor lineage.** The object is ~50% B among resolved
   cells, so NK could not exceed ~0.34 purity and was losing its label entirely. Enrichment
   over the family baseline is capped at 1/baseline and so cannot certify the dominant
   family; purity cannot certify a minor one. A cluster now qualifies on **either** route.
   This is what recovered NK (48,953) and first produced a Plasma cell label.
3. **The 19b gate shipped broken and reported a vacuous PASS.** Its negprobe term was a
   median of per-cell ratios; `nCount_negprobes` is nonzero in only 21.8% of cells, so that
   median is 0 everywhere and `negrate > baseline` read `0 > 0` -- FALSE always, inside an
   AND. Now pooled, and the gate is two conditions rather than a four-way AND.
4. **The FDC rule needed three conditions, not one.** CR2/FCER2 over B_core inflates
   automatically wherever B signal is thin; on the first pass it relabelled a 102,783-cell
   cluster with argmax B purity 0.806 as FDC. It now also requires B_core depressed and
   argmax not strongly B -- per-cell marker evidence must not be overridden by a ratio.
5. **Cluster labels are not stable.** A clean 107,894-cell naive/mantle B cluster (IGHD
   6.05, SELL 4.60) present in the first build vanished into the B mass after 2.9% of cells
   were removed. It is not in the shipped vocabulary. This is the strongest argument in the
   project for carrying argmax and cluster labels as separate columns.
6. **The halo QC fired on the label that had just been added.** Plasma cell scores **12.69x**
   -- worse than the 8.4x step 18 measured for InSituType -- so the cluster-level PC label is
   neighbourhood-assembled and is barred from statements. Use `celltype_argmax` for plasma
   cells. Every other label is clean (Myeloid 2.87, Endothelial 2.18, Stromal/FDC 2.17,
   T 1.47, B 1.40, NK 0.95).

### Statement rule, as implemented

`cluster_status` in {confident, weak_support, module_call} AND `is_niche == FALSE` AND
`halo_suspect == FALSE`. All three columns travel in
`hc_cells_labelled_high.csv.gz`, not only in a log.

## Follow-on

`CLAUDE.md` needs three edits once this lands: fact 1 (merged labels valid on a defect-free
object, with the reason), a new fact recording that both callers are neighbourhood-contaminated
and that clustering is the more contaminated of the two, and the closed-avenue list gaining
compartment-restricted calling.
