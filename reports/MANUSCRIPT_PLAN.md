# P3_CosMx_v2 — manuscript readiness review

**Revised 2026-09-08.** Supersedes the 2026-09-01 version. Every number is quoted from an
existing report or result table and the source is named inline; nothing here is newly
computed. The HTML companion is `manuscript_plan.html` in this directory.

Sources: `output/13_summary/IB1_summary.html` (rev. 2026-09-01),
`output/16_summary/step16_FL_comparison.html`, `output/17_summary/step17_PD1_PDL1.html`,
`output/19_summary/step19_cell_identity.html`, `output/21_summary/step21_interface.html`,
`output/24_summary/step22_24_gap_closure.html`, and the tables cited per line.

---

## Revision notice — what changed since 2026-09-01

**Four of the eight gaps are closed and are no longer listed below.** Gap numbering is
preserved, so the remaining gaps keep their original identifiers: **G2, G3, G6, G8**.

| gap | closed by | verdict |
|---|---|---|
| **G1** — step-19 labels terminal | steps 24 + 24b | **R9 CORRECTED.** Callers agree on 70.2% of cells; step 16's "RB1 sits at chance" was the caller's silence, not RB1's architecture. |
| **G4** — no ligand-receptor analysis | step 22 | Co-expression in proximity is measurable and **small** (1.01–1.18x). Only *CXCL13*–*CXCR5* is distance-robust. The largest value in the table is leakage-shaped. |
| **G5** — no copy-number inference | step 23 | **9p24.1 dosage NOT obtained.** The no-Hodgkin control fires and the signal is not confined to one lineage. |
| **G7** — FL/Hodgkin interface | step 21 | A **bounded null** in expression and composition; lesion arrangement **unresolved**, not negative. |

**One established result changed, and it propagates.** R9 below is rewritten, and the draft
abstract and Figure 6 are revised with it: the claim that the recurrence had *lost* spatial
organisation does not survive step 24b.

---

## 1. What the project has established

### 1.1 Results that survive their own scrutiny

| # | Result | Key numbers | Source |
|---|---|---|---|
| R1 | A Hodgkin-specific chemokine field with a **measured extent** | IB1 decays 2.71x [1.95, 3.91] from 100–400 um to >1600 um from a peak; RB1 1.20x [1.14, 1.26]; ratio **2.27x [1.62, 3.25]**. Profile 75.7 -> 47.6 -> 29.4 -> 12.1 cp10k, then flat | `17_pd1/matched/field_extent_verdict.csv` |
| R2 | The hotspots have **no RB1 equivalent** at a symmetric QC tier | 80 hotspots in IB1 (66 seam-free), **zero** RB1 bins clear the same cut | step 16, `matched` tier |
| R3 | The elevation is **not composition** | holds within every lineage; composition standardisation raises the ratio | step 11 |
| R4 | Hotspots are **DC territory, not macrophage** | cDC1 1.83x, mature DC 1.59x, macrophage 0.77x, `mac_core` 0.57x | step 11 |
| R5 | **BATF3 / CD40 / IL6** mark the field, and it is not spillover | BATF3 +2.07 log2 (q = 3e-37); BATF3 ratio **highest farthest from any cDC** (5.14x at >120 um, detection 3.08x) | `16_checkpoint/batf3_spillover.csv` |
| R6 | A **CCR4 gradient without a Treg signature** (pre-declared readout) | CCR4 rho 0.415 in FL foci (q = 1.4e-4), 0.431 at FOV level, vs 0.120 tissue-wide; FOXP3 rho -0.007 (rank 2534/6175), IL2RA negative | `16_gradient/matched/ccr4_axis.csv` |
| R7 | **Ligand-side-only** checkpoint elevation in IB1's hotspots | CD274 OR 1.83, PDCD1LG2 OR 1.92 (clean reference), PDCD1 1.32 ns. PD-1 axis module +0.472 (q = 2.1e-5), 9p24.1 +0.391 (q = 2.6e-5), **IFN-gamma flat** (-0.114 ns) | `binomial_model.csv`, `module_contrasts.csv` |
| R8 | RB1 recovers the **canonical follicular PD-1 pattern of FL**, unsought | PDCD1 in T, follicular vs interfollicular, OR 1.82 (q = 4.6e-14) on 982 positive cells; CTLA4 1.67, CXCL13 1.67, TOX 1.44, TIGIT 1.38. Tfh definition excludes PDCD1 | `binomial_model.csv` |
| **R9** | **REVISED.** Both nodes are non-randomly arranged; **IB1 about twice as strongly**. The earlier "RB1 is spatially homogeneous" is withdrawn | At matched cell counts, median \|log2 enrichment\| **IB1 0.265 vs RB1 0.128**. RB1's fraction-at-chance is 0.818 under step-16 labels but **0.535** under step-19's, and matching *n* does not close the gap — so step 16 measured the argmax caller's silence, not RB1's architecture | `24_arch_relabel/matched/control_effect_size.csv` |
| R10 | **Tfh tracks FL B cells in both nodes** | 0.465 (IB1) / 0.416 (RB1) of the same-sample B enrichment; leads the next population 3–4x in each independently | `16_composition/matched/` |
| R11 | Regional signal is **compositional, not cell-state** | Conditioned on cell type: **zero** significant genes within B or within T for FL focus vs background; the BATF3 hotspot signal also vanishes within B and T | step 16, contrast CT1 |
| R12 | IB1's B-dense tissue carries **the Hodgkin field, not the FL programme** | IB1: CCL17 21.8 -> 7.3, CCL22 11.4 -> 3.3; BCL2/BCL6/MME flat. RB1: chemokines flat, BCL2 -0.70, BCL6 -0.70, MME -0.74 | `16_rings/matched/ring_trends.csv` |
| R13 | The **over-fixation lesion is characterised with a mechanism** | 37.1% of IB1 vs 8.1% of RB1; 1.11x counts at 1.03x area; max log2FC 0.98 vs 3.73; 0.61% vs 1.55% of budget to markers. 9 of 9 criteria fired | `07_overfixation/ib1_overfixation_summary.md` |
| R14 | A **two-column cell-identity framework** with held-out validation | 2.05M high-confidence cells, 10 labels; 66.2% argmax/cluster agreement; every label tops its own held-out marker set; rescue to 61.5% coverage | `19_summary/step19_cell_identity.html` |
| **R15** | **NEW.** The FL does not respond to the Hodgkin field distinguishably from ordinary IB1 tissue — a **bounded null** | 0 of 44 pre-declared genes in each of B/T/myeloid; 0 of 6 cell types on composition. Median resolvable effect **1.09–1.14 fold** on 20–22 FL FOVs | `21_summary/step21_interface.html` |
| **R16** | **NEW.** The two lesions' arrangement is **UNRESOLVED**, not negative | Observed hotspot-to-FL 949 um vs null median 1105 (p = 0.53); under a defect-matched null, 1503 um (p = 0.152; 0.085 the other direction). Step 12's own statistic with FOV as the unit gives p = 0.332 vs the published 0.233 | `21_interface/matched/arrangement_verdict.csv` |
| **R17** | **NEW.** Ligand-receptor co-expression in proximity is real and **small** | 1.01–1.18 fold over 2,315 curated pairs; only *CXCL13*–*CXCR5* enriched across **all four** distance bands (1.14/1.16/1.15/1.13). *CCL17*/*CCL22*–*CCR4* is **below 1 inside the FL foci** — the CCR4 gradient is regional, not pairwise | `22_lr/matched/lr_proximity_tier1.csv` |
| **R18** | **NEW.** IB1's hotspots carry an elevated 9p24.1 neighbourhood, but it is **not dosage** | Rank 5/439 (T) and 5/443 (fibro) of 440+ cytobands, non-immune genes rising more than immune. But RB1's fibroblasts also fire (q = 0.039) and the signal spans 3 of 4 lineages when a clonal event must be confined to one | `23_dosage/matched/dosage_9p24.csv` |

### 1.2 Methodological results — arguably the most transferable content

- **Depth-vs-positivity AUC ~0.70 is generic** to every sparse gene on this platform
  (IGKC 0.693, CD3E 0.695, CD68 0.731), so it is not evidence a gene is noise.
- **cp10k is compositional** and fakes depression in the panel's most abundant genes; gate
  every call on detection rate agreeing in direction.
- **Spearman's rho cannot separate a decaying field from a slide offset.** The **amplitude**
  of the decay discriminates; its monotonicity does not.
- **Pseudoreplication with a new name.** Five T compartments have ICC 0.875 -> 1.11 effective
  observations per gene; a sign-consistency *p* of 2.7e-12 becomes **p = 0.143**.
- **The between-slide offset bounds the design** (reference-zone ORs to 5.33) — but **gene by
  gene**: CD274's reference contrast is null, so PD-L1 is underpowered, not foreclosed.
- **Both cell-type callers are neighbourhood-contaminated**, clustering the worse (6.8x vs
  2.5x proximity gradient); leakage is real at the population level and **unmeasurable per
  cell** (AUC 0.35–0.61).
- **A neighbourhood gradient is not itself evidence of leakage.** Rising lines are the
  expected null.
- **Merged clustering launders defect tissue into lineage labels** unless defect cells are
  dropped per cell *before* clustering.
- **NEW — the over-fixation defect is spatially confounded with almost everything in IB1**,
  and per-cell removal does not make the surviving area uniform: per-bin defect fraction is
  0.33 in FL bins, 0.45 in background, 0.59 in hotspots; it is a rim (rho -0.578 with edge
  distance) and therefore rises with distance to FL (rho +0.379) and with the chemokine axis
  (rho 0.412). Carry it as a covariate. IB1's "tissue pieces" are 2 pre-drop and 3 post-drop —
  coverage, not anatomy.
- **NEW — a percentile bootstrap p needs the `(1 + count)/(1 + n)` correction.** Unfloored it
  returns exactly 0 and BH reads that as infinite significance; in step 21 that alone
  manufactured **143 genome-wide "hits"**.
- **NEW — composition is a SHARE, never a count.** Endothelium is the diagnostic: rho -0.378
  on counts, **+0.236** on shares.
- **NEW — a permutation null often has a closed form.** Step 22's first design was ~24 million
  permutations; the same quantity is two sparse matrix products.
- **NEW — a verdict must consult its own negative control.** Step 23's first verdict did not
  look at RB1 and labelled RB1's own fibroblasts "DOSAGE-CONSISTENT".
- **NEW — a 1.0x-scale enrichment must be looked at, not just tabulated.** Step 22's snapshots
  show the receptor side sits at the 2-count positivity floor in 8 of 9 representative panels
  even for the strongest pair.

### 1.3 Closed negatives that constrain the paper

Per-cell HRS/CD30 calling. Light-chain regionalisation of RB1. Cross-sample neighbourhood
matching (2 of 89 bins). The difference-in-differences (0 of 30, MDE 1.8–2.5 fold).
Transcript-level PD-1/PD-L1 proximity (**NOT DELIVERED**). Segmentation fragmentation of HRS.
Cluster-based T-cell subsetting. Compartment-restricted cell typing. HCA tonsil label
transfer. **NEW:** 9p24.1 copy-number inference; distance-to-hotspot as a field-exposure axis;
the FL-response-to-field interaction (bounded null).

---

## 2. Gaps and missing analyses

Four gaps are closed (see the revision notice) and are not repeated. **Numbering is preserved.**

### G2 — There is no morphological or protein confirmation of anything  **[BLOCKING — top priority]**

Five morphology channels are on the slide and `output/15_imagej/` holds a reviewable export
(overview, 34 full-resolution crops, ROIs for 29 hotspots and the FL neighbourhoods).
**No reviewed verdict is recorded anywhere in `output/`.**

The central claim is a chemokine field radiating from focal peaks that have never been shown
to contain Reed–Sternberg cells. As it stands the hotspots are *defined* by a chemokine
measurement and *validated* by the same chemokine measurement.

- **No pathologist read of the hotspot ROIs.** The export exists for this and was never used.
- **No RB1 ImageJ export** — an image-side positive with no image-side negative control.
- **No orthogonal assay.** CD30 (and PD-L1) IHC on the block is outside this dataset but is
  the obvious cheap anchor, and it does not have the panel's dynamic-range problem.

**What is new since 2026-09-01:** step 22's snapshots
(`output/22_lr/matched/snapshots/figures/`) prove the rendering pipeline works end to end on
the real OME-TIFFs — tight crops with segmentation boundaries and transcripts on 5-channel
morphology. The same machinery would render the hotspot ROIs for review, so the remaining
cost of G2 is a pathologist's time, not engineering.

### G3 — The 4.07x far-field IB1:RB1 residual is unattributed  **[BLOCKING for R1's framing]**

Beyond 1600 um from any peak IB1 is still 4.07x RB1 and flat in distance — something a
diffusible field cannot produce. Two live explanations, a between-slide offset or **diffuse
sub-threshold HRS cells**, and the far-field coherence test built to separate them is
uninformative (the no-Hodgkin control is *more* clustered). Cleaning the reference did **not**
shrink it, which excludes field contamination as the cause.

Untried: a negprobe-anchored technical-offset estimate; the morphology channels (diffuse HRS
cells are large and CD45-dim and might be visible without transcript calling); an FOV-level
offset model fit on housekeeping genes.

### G6 — RB1 is analytically under-served  **[MEDIUM]**

1.4M cells that mostly serve as a null — and step 24b makes this sharper, because RB1 turns
out to have real architecture that the old labels could not see.

- **Step 14 never ran on RB1** — `output/14_lightchain_in_out/` contains `IB1/` only.
- The **12,394 light-chain-evaluable RB1 cells with no neighbourhood** have never been
  examined. If they are spatially clustered they are RB1's only internal reference. (Note the
  corrected denominator: 466,764 RB1 cells were light-chain-evaluable, not 1.29M — the older
  "the clone *is* the node" phrasing overstated it ~3x.)
- No ImageJ export (see G2).

### G8 — Composition is reported without uncertainty  **[CLOSED 2026-09-08 — step 25a]**

`25_script_composition_fov.R` computes the per-FOV share of every label among high-confidence
cells; `output/25_figures/composition_fov_summary.csv` and Figure 1d carry the spread.

**And it qualifies a phrasing.** IB1's composition is far less stable across the section than
RB1's: its two largest fractions have almost completely overlapping IQRs — myeloid 17.7–34.1%
against B cell 17.5–35.2% — so **"IB1 is myeloid-dominant" is not supportable**, though
"IB1 23.6% B vs RB1 71.4%" is (RB1's B IQR 66.1–77.7%, disjoint). The pooled myeloid estimate
also sits 4.4 points above the FOV median, so a few dense fields carry it.

### Lower priority

- **No external reference anchor.** A published cHL or FL spatial/scRNA dataset would
  independently support the DC (R4) and Tfh (R10) calls.
- **Neutrophil rests on one gene** (CSF3R). Fold it into `unresolved` rather than shipping a
  permanently caveated label.
- **Treg has no independent validator on this panel.** The paper should avoid the word; R6
  already does.
- **No point-pattern statistics** (cross-K, pair correlation). The density-stratified
  permutation is arguably better, but a reviewer may ask.
- **Bin-size sensitivity** is swept for hotspots but not for the architecture/region axis that
  R9–R12 rest on.
- **The `arch` axis means different biology in the two samples.** Cleaner to define IB1's FL by
  light chain throughout and use the arch axis only for RB1.
- **NEW — minor-lineage identity is a coin flip between callers.** NK 32.1% and PC 37.5%
  agreement (step 24). Any claim resting on NK or PC identity needs both callers to agree, or
  needs dropping.

### Out of reach — state as limitations, do not chase

n = 1 patient, two sections; the between-section offset is not estimable from two sections.
Per-cell HRS calling on this panel. Cross-sample regional contrasts for the offset-affected
genes. Spatial engagement of the PD-1 axis. Per-cell leakage correction. 9p24.1 copy number
from a 6K targeted panel.

---

## 3. Proposed abstract

**Framing decision unchanged.** **(A)** the biology — a composite node where the Hodgkin
component projects a bounded chemokine and checkpoint-ligand field over adjacent tissue — is
the stronger paper, with the methodological discipline visible in its limitations and
supplement.

### Abstract (framing A, ~305 words)

> **Composite lymphoma places two diseases in one microenvironment.** We profiled a lymph
> node containing classic Hodgkin lymphoma (cHL) together with follicular lymphoma (FL), and
> the same patient's subsequent FL recurrence, by single-molecule spatial transcriptomics
> (CosMx 6K; 4.4 million segmented cells, two sections, one patient).
>
> Reed–Sternberg cells were not identifiable per cell on this panel — *TNFRSF8* reaches six
> counts with 78% singletons, and the Hodgkin-free recurrence was 6.27% CD30-positive against
> the composite node's 7.01% — so we replaced per-cell calling with a regional readout.
> *CCL17*/*CCL22* defined 80 focal hotspots in the composite node and none in the recurrence
> at an identical threshold and a symmetric quality tier. The elevation is a spatially
> decaying field rather than a section-level offset: chemokine falls 2.71-fold [1.95–3.91]
> between 100–400 µm and beyond 1600 µm from a peak, against 1.20-fold [1.14–1.26] in the
> recurrence, with most of the decay spent by 400–800 µm, and it persists within every
> lineage. Hotspots are dendritic-cell rather than macrophage territory and carry *BATF3*,
> *CD40* and *IL6*, with *BATF3* enrichment highest farthest from any dendritic cell,
> excluding transcript spillover. The field is accompanied by a *CCR4* gradient roughly
> 3.5-fold steeper inside FL foci than tissue-wide, with *FOXP3* flat — chemoattraction
> without a regulatory T-cell signature. Hotspots carry elevated *CD274* and *PDCD1LG2* with
> *PDCD1* and interferon-γ response flat: ligand-side activation, the pattern 9p24.1 copy
> gain predicts, though a direct dosage test on this panel could not confirm gain. The
> recurrence independently reproduced the canonical follicular PD-1 pattern of FL. Under a
> composition- and density-preserving permutation both nodes departed from random cellular
> arrangement, about twice as strongly in the composite node.
>
> Conditioned on cell type, no gene distinguished FL foci from background, and FL tissue's
> response to the Hodgkin field was indistinguishable from ordinary tissue's. The Hodgkin
> component therefore remodels the *composition and milieu* of neighbouring tissue over
> several hundred micrometres without altering the transcriptional state of individual
> bystander cells. Findings are descriptive (n = 1).

**Two sentences changed from the 2026-09-01 draft**, both because later work corrected them:
the architecture sentence no longer says the recurrence had *lost* organisation (step 24b),
and the 9p24.1 clause no longer asserts the mechanism (step 23).

### Title options

1. *A Hodgkin chemokine field with a measured extent remodels adjacent follicular lymphoma in
   a composite lymph node*
2. *Ligand-side checkpoint activation within a spatially bounded Hodgkin chemokine field*
3. *Composition, not cell state: what a Hodgkin lesion does to the tissue around it*

---

## 4. Figure list with legends

Status — **[have]** exists as a PNG; **[make]** data exists, figure does not; **[blocked]**
needs work not yet done, with the gap named.

### Figure 1. The case, the tissue, and the analytical substrate
- **(a)** Whole-section 5-channel morphology, both samples, with ROI overlays **[blocked · G2]** — IB1 exists as `15_imagej/IB1_overview_L5.tif`; RB1 never exported
- **(b)** `19_celltypes/figures/celltype_spatial.png` **[have]**
- **(c)** `19_celltypes/figures/qc_spatial.png` **[have]**
- **(d)** `output/25_figures/fig_25_1d_composition_fov.png` **[have]** — G8 closed by `25_script_composition_fov.R`

> **Figure 1. A composite Hodgkin/follicular lymphoma node and its follicular recurrence.**
> (a) Whole-section morphology; cyan outlines mark *CCL17*+*CCL22* hotspots, orange the
> light-chain-defined FL neighbourhoods. (b) Cell-type assignment across both sections, from a
> merged clustering of quality-passing, non-defect cells with per-cell marker argmax carried as
> a separate column; colour scales are per panel and not comparable between types. (c) Two
> distinct technical geographies: IB1's over-fixation defect forms a peripheral rim while
> quality-tier failures form an FOV grid. Cells failing either are excluded per cell before
> clustering. Note the defect is spatially confounded with the regions analysed later — 0.33
> defect fraction in FL bins, 0.45 in background, 0.59 in hotspots — and is carried as a
> covariate throughout. (d) Composition by sample among high-confidence cells; bars show the
> FOV-level interquartile range. RB1 is 71.4% B cells against IB1's 23.6%.

### Figure 2. Cell identity: two callers, held-out validation
- **(a)** `19_celltypes/figures/celltype_umap_high.png` **[have]**
- **(b)** `19_celltypes/figures/marker_heatmap.png` **[have]**
- **(c)** `19_celltypes/figures/composition_figure.png` **[have]**
- **(d)** `output/25_figures/fig_25_2d_caller_agreement.png` **[have]**

> **Figure 2. Cell identity is carried in two columns, and the disagreement is informative.**
> (a) UMAP of 2.05M high-confidence cells. (b) Marker profiles: left, the genes that named the
> labels, circular by construction; right, canonical markers no labelling rule consulted —
> every one of the ten labels tops its own held-out set, and Plasmablast ranks second on *both*
> parent sets. (c) Per-label agreement between cluster label and per-cell argmax: 66.2%
> overall, with 43.4% of cells receiving no argmax call at all, which is why both columns are
> kept. (d) Agreement between the two callers by lineage, on the 998,544 cells both name: 70.2%
> overall, but **NK 32.1% and plasma cell 37.5%**. Minor-lineage identity is a coin flip
> between callers and no claim should rest on it alone.

### Figure 3. A Hodgkin-specific chemokine field with a measured extent
- **(a)** `11_hodgkin_niche/figures/fig_11_niche_map_IB1.png` + `16_regions/matched/figures/fig_16_regions_*.png` **[have]**
- **(b)** `17_pd1/matched/figures/fig_17_field_extent.png` **[have]**
- **(c)** `11_hodgkin_niche/figures/fig_11_confound_within_lineage.png` **[have]**
- **(d)** `11_hodgkin_niche/figures/fig_11_hotspot_myeloid_vs_dc.png` **[have]**
- **(e)** Annotated morphology crop of a hotspot ROI **[blocked · G2]** — crops exist, annotation does not

> **Figure 3. The Hodgkin component projects a chemokine field that decays over 400–800 µm.**
> (a) *CCL17*+*CCL22* per 100 µm bin; 80 hotspots in the composite node, **none** in the
> recurrence at the same threshold. Sequential scale clamped to the 5th–99th percentile.
> (b) Decay from each sample's own top-1% bins, so the construction is identical in both and
> the recurrence is a no-Hodgkin control: IB1 falls 2.71-fold [1.95, 3.91], RB1 1.20-fold
> [1.14, 1.26], ratio 2.27-fold [1.62, 3.25] on an FOV cluster bootstrap. Spearman's rho does
> *not* discriminate (−0.180 vs −0.095) — the amplitude of the decay is the informative
> quantity, not its monotonicity. Beyond 1600 µm a flat 4.07-fold residual remains,
> unattributed. (c) The elevation holds within every lineage. (d) Hotspots are dendritic-cell
> rather than macrophage territory: cDC1 1.83-fold, mature DC 1.59-fold, macrophage 0.77-fold.
> (e) Morphology at a representative hotspot.

### Figure 4. What the field does to the tissue around it
- **(a)** `16_gradient/matched/figures/fig_16_field_gradient.png` **[have]**
- **(b)** `output/25_figures/fig_25_4b_batf3_vs_cdc_distance.png` **[have]**
- **(c)** `output/25_figures/fig_25_4c_ccr4_axis.png` **[have]**
- **(d)** `16_rings/matched/figures/fig_16_fl_edge_rings.png` **[have]**

> **Figure 4. The field is transcriptionally coherent and reaches into the FL tissue.**
> (a) Genome-wide association with the field axis independently recovers the Hodgkin
> programme: *BATF3* 0.555, *CD40* 0.389, *IL6* 0.335 lead the positive list, the negative list
> being normal T-zone architecture. (b) *BATF3* enrichment is *highest* farthest from any cDC
> (5.14-fold beyond 120 µm), the opposite of what spillover predicts. (c) *CCR4* was declared
> as the readout before any genome-wide scan; its association with the field is ~3.5-fold
> stronger inside FL foci (rho 0.415, q = 1.4e-4) than tissue-wide, while *FOXP3* is flat and
> *IL2RA* negative — a *CCR4* gradient without a regulatory T-cell signature. Note this is a
> **regional** association: at the level of individual cell pairs, *CCL17*/*CCL22*–*CCR4*
> co-expression in proximity is *below* chance inside the FL foci (Figure 7). (d) Across the FL
> boundary the two nodes diverge: the recurrence grades its germinal-centre programme, the
> composite node grades the chemokines instead.

### Figure 5. Checkpoint: ligand-side in the Hodgkin component, receptor-side in the FL
- **(a)** `output/25_figures/fig_25_5a_module_contrasts.png` **[have]** — cross-sample facets marked as bounded by the design
- **(b)** `output/25_figures/fig_25_5b_forest_ib1_hotspot.png` **[have]**
- **(c)** `output/25_figures/fig_25_5c_forest_rb1_follicular.png` **[have]**
- **(d)** `output/25_figures/fig_25_5d_did_vs_reference.png` **[have]**
- **(e)** `output/25_figures/fig_25_5e_9p24_dosage.png` **[have]**

> **Figure 5. The two lesions engage the PD-1 axis from opposite ends, and the mechanism is
> not established.** Positivity is ≥2 counts, modelled per cell with region and depth decile as
> covariates; intervals from resampling fields of view. IB1's reference is `background`
> restricted to >800 µm from any hotspot. (a) The hotspots carry the PD-1 axis (+0.472 log2,
> q = 2.1e-5) and the 9p24.1 module (+0.391, q = 2.6e-5) while the interferon-γ response does
> *not* rise. (b) Per cell the elevation is ligand-side only: *CD274* OR 1.83, *PDCD1LG2* 1.92,
> *PDCD1* 1.32 (ns). (c) The recurrence shows the reciprocal pattern in the expected cells:
> *PDCD1* in T cells 1.82-fold more likely inside follicles (q = 4.6e-14), recovered without
> being sought. (d) Cross-sample contrasts shown beside the two sections' *reference* contrast:
> where the reference is large the apparent difference is the section (*CXCL13* 5.21 -> 0.98).
> *CD274* is the exception — reference null, residual 1.76-fold, unresolved rather than refuted.
> (e) **The dosage test does not support 9p24.1 gain.** The band scores in the 99th percentile
> of 440+ cytobands in IB1's hotspots with non-immune 9p24.1 genes rising too — but the same
> statistic fires in the no-Hodgkin recurrence (q = 0.039) and the signal spans three of four
> lineages, when a clonal copy event must be confined to one. The module stays exploratory.

### Figure 6. Architecture: both nodes are organised, the composite node more strongly
- **(a)** `16_cooccur/matched/figures/fig_16_cooccurrence.png` **[have]**
- **(b)** `output/25_figures/fig_25_6b_architecture_both_callers.png` **[have]**
- **(c)** `16_composition/matched/figures/fig_16_composition_did.png` **[have]**
- **(d)** `16_tfh_views/matched/tfh_*.png` **[have]**

> **Figure 6. Composition held fixed, the two nodes are arranged differently — but neither is
> random.** The null permutes cell-type labels within local-density deciles, so composition and
> geometry are fixed and only the assignment of types to positions is randomised. (a) Pairwise
> co-occurrence enrichment by sample and region: in IB1's background Tfh self-aggregate at
> 2.74-fold chance, preferentially neighbour B cells at 1.62-fold and avoid macrophages at
> 0.47-fold. (b) **A correction.** An earlier analysis reported the recurrence as sitting at
> chance throughout. Repeating it under an improved cell-type caller, with every lineage
> down-sampled to matched cell counts so power is equalised, the recurrence's enrichments sit
> **2.15-fold further from chance** than the old labels showed, and the same pair is further
> from chance under the better labels in 80.3% of comparisons (p = 2.9e-39). Both nodes are
> non-randomly arranged; the composite node about twice as strongly (median |log2 enrichment|
> 0.265 vs 0.128). The earlier reading measured the caller's silence, not the tissue.
> (c) Tfh is nonetheless the only population that follows the FL B cells at a comparable rate
> in both nodes. (d) Registered morphology at representative Tfh-rich regions.

### Figure 7. NEW — What a small proximity enrichment actually looks like
- **(a–b)** `22_lr/matched/snapshots/figures/fig_22_snap_CXCL13CXCR5_{cherry,repr}.png` **[have]**
- **(c)** `22_lr/matched/snapshots/figures/fig_22_snap_PDCD1LG2PDCD1_repr.png` **[have]**
- **(d)** `output/25_figures/fig_25_7d_lr_distance_profile.png` **[have]** — built on `enrich_fov`, the quantity step 22 reports

> **Figure 7. Ligand–receptor co-expression in proximity, and why the largest number in it is
> an artefact.** Co-expression in proximity is measured, **not** signalling or engagement:
> spatial engagement of a receptor–ligand axis remains unmeasured in this project after three
> null designs failed their own positive controls. (a, b) *CXCL13*–*CXCR5*, the strongest pair,
> shown as the clearest available examples and as a uniform draw from the same 4,882-pair
> candidate pool. Even here the receptor side sits at the 2-count positivity floor in 8 of 9
> representative panels — a "pair" is typically a well-expressed ligand cell beside a barely
> positive receptor cell. (c) *PDCD1LG2*–*PDCD1*, the largest enrichment in the whole table
> (2.28-fold): its entire candidate pool is 10 pairs in 6 fields of view, every one 2-versus-2
> counts, with cell-type combinations (B–B, endothelial–stromal, myeloid–NK) that are not the
> myeloid/DC-to-T axis PD-L2 signals on. (d) The discriminator: enrichment across five distance
> bands. *CXCL13*–*CXCR5* persists to 40–80 µm; *PDCD1LG2*–*PDCD1* runs 2.28 -> 1.54 -> 1.14 ->
> 1.01, the monotone decay that transcript leakage produces. Never quote a single band.

### Supplementary figures

| # | Content | Status |
|---|---|---|
| S1 | Over-fixation diagnostic: nine criteria, channel dissociation, radial gradient | **[have]** `07_overfixation/figures/*` |
| S2 | QC attrition and the two tiers | **[have]** `output/25_figures/fig_25_S2_qc_attrition.png` |
| S3 | The measurement floor: depth AUC, module ceilings, MDE per contrast | **[have]** `output/25_figures/fig_25_S3_measurement_floor.png` |
| S4 | Non-identifiability of the cross-sample FL comparison | **[have]** `output/25_figures/fig_25_S4_composition_disjoint.png` |
| S5 | Neighbourhood contamination: leakage profiles, the halo test, plasma-cell conditioning | **[have]** `18_leakage/figures/*`, `halo_qc_high.png` |
| S6 | Retracted and closed analyses: sign-consistency pseudoreplication, CD30 dynamic range, transcript proximity not delivered | **[have]** `fig_17_sign_consistency.png` |
| S7 | Light-chain neighbourhoods, IB1 | **[have]** IB1 only — **RB1 half missing · G6** |
| S8 | The IB1 B-marker deficit is the label, not the biology | **[have]** `b_discrepancy.png` |
| S9 | Seam, scale and threshold robustness of the hotspot definition | **[have]** `fig_11_seam_check.png` |
| **S10** | the over-fixation defect is spatially confounded | **[have]** `output/25_figures/fig_25_S10_defect_confounding.png` |
| **S11** | **NEW** — axis selection for the FL/Hodgkin interface, and the two arrangement nulls | **[have]** `21_interface/matched/figures/fig_21_{axis_selection,arrangement}.png` |
| **S12** | **NEW** — the FL field response is a bounded null: 44 pre-declared genes × 3 cell types with the detectable-effect band drawn | **[have]** `fig_21_interface_response.png` |
| **S13** | **NEW** — the remaining 13 ligand–receptor snapshot figures | **[have]** `22_lr/matched/snapshots/figures/*` |

---

## 5. Recommended order of work before submission

1. **G2** — get the hotspot ROIs in front of a pathologist, and export RB1. The rendering
   pipeline is proven (step 22's snapshots); what remains is a pathologist's time. This is what
   a reviewer will demand first, and it is the one gap that would convert the central claim
   from "a chemokine field with focal peaks" into "a chemokine field centred on HRS cells".
2. ~~**The [make] panels.**~~ **DONE** — all 15 built by `scripts/25_script_figure_panels.R`
   into `output/25_figures/`. Plotting only: each reads a saved result table and recomputes
   nothing, so no panel can drift from the result it depicts.
3. **G3** — the far-field residual. A reviewer will ask what fraction of "the field" is a field.
5. **G6** — RB1's missing step 14 and the 12,394 unexamined evaluable cells.
6. **Lower-priority items**, of which the most substantive is dropping or double-gating any
   claim that rests on NK or plasma-cell identity.
