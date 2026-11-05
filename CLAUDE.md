# P3_CosMx_v2

CosMx 6K spatial transcriptomics, **n = 1 patient, two slides**. Every result is
descriptive; nothing here is population-level inference.

| | IB1 | RB1 |
|---|---|---|
| What it is | **composite lymph node: classic Hodgkin (cHL) + follicular lymphoma (FL)** | **the FL recurrence — same clone, no Hodgkin component** |
| Role | the case | negative control for Hodgkin signal; positive control for full FL involvement |

## Read these first

These are generated documents and are the project's institutional memory. Read them
before proposing analysis — most questions are already answered, and several avenues
are closed with evidence.

- `output/13_summary/IB1_summary.html` (source `scripts/13_summary_report.Rmd`) — the
  project report for steps 07-12: what was asked, what held, what failed, and the retracted
  results. **REVISED 2026-09-01** and it carries a dated "Revision notice" changelog at the
  top, because an earlier copy was already circulated to the primary researchers. Four
  findings were superseded by steps 16-17: the PRE-QC negative control (FAILED -> resolved),
  the CD30 depth-AUC leg (withdrawn, conclusion unchanged), "tissue-wide" field (now measured),
  and the FL-foci null (power -> compositional). **Keep the filename stable and keep the
  changelog** — collaborators need to tell their copy from this one.
- `output/07_overfixation/ib1_overfixation_summary.md` — the IB1 over-fixation
  diagnostic (InSituType cluster `i`, 37.1% of IB1 vs 8.1% of RB1).
- `output/16_summary/step16_FL_comparison.html` (source `scripts/16_summary_report.Rmd`) —
  the IB1-vs-RB1 FL comparison. Read its "What this step cannot support" section before
  proposing any cross-sample regional analysis: the direct contrast is **not identifiable**
  (IB1's most B-dense quartile is less B-rich than RB1's least B-dense quartile, so 2 of 89
  FL bins find a composition-matched partner), and the difference-in-differences that was
  meant to rescue it is withdrawn.
- `output/24_summary/step22_24_gap_closure.html` (source `scripts/24_summary_report.Rmd`) —
  steps 22, 23 and 24 together. Read before proposing a ligand-receptor analysis, any 9p24.1
  copy-number claim, or any re-run of the step-16 architecture result on step-19 labels.
- `output/21_summary/step21_interface.html` (source `scripts/21_summary_report.Rmd`) — the
  FL/Hodgkin interface, the design step 17 asked for. Read fact 9 before proposing any
  within-IB1 regional analysis: the over-fixation defect is spatially confounded with the
  regions, with the chemokine axis and with distance-to-FL, so per-cell removal is not enough.
- `output/17_summary/step17_PD1_PDL1.html` (source `scripts/17_summary_report.Rmd`) — the
  checkpoint question, exhausted. Read facts 3 and 5 below before proposing any cross-sample
  contrast at all; step 17 bounds the two-section design, not just the FL question.

## Ten things that are not derivable from the code

1. **Cluster labels come through `scripts/lib/cluster_labels.R`, never from HierAnnot.**
   HierAnnot's `annot_path` is wrong for a substantial minority of clusters on this data
   (stroma called macrophage, endothelium called neutrophil, pDC called plasma cell).
   `scripts/09b_script_marker_relabel.R` holds the audit. The lookup `stop()`s rather
   than falling back, on purpose.

   **The audit covers IB1, RB1 and merged, at three clusterings each** —
   `pearson_clusters_low`, `pearson_clusters_high` and `inSitu_clusters` (the vendor
   InSituType call). Only IB1 has HierAnnot labels to compare against, so RB1's and
   merged's `agrees` column is empty: marker calls with no independent comparator.

   **Do not use the merged labels for cell identity.** `09c_script_cluster_consistency.R`
   shows why: merged `pearson_clusters_low` cluster 1 is 93.6% IB1 cells that IB1's own
   audit calls `low-information, defect_dominated` at 70.1% defect — pooling dilutes the
   defect fraction to 42%, under the 50% `DEFECT_MAX` guard, so it escapes the guard and
   acquires the label "Neutrophil". Merged cluster 0 does the same for RB1's unresolved
   mass and gets "B cell". **Pooling launders defect-dominated tissue into a lineage
   label.** Per-sample audits decline both. Cross-sample cell type comes from the argmax
   marker caller (`output/16_celltypes/`), which is now a checked result and not an
   assumption.

   **The ban is about DILUTION, not about merging, and step 19 lifts it under one
   condition.** The failure needs defect cells present to dilute. `output/19_celltypes/`
   drops them per cell BEFORE clustering, so there is nothing left to launder, and its
   merged labels are valid — the gate in `19_script_cluster_health.R` confirms no large
   low-signal cluster survives. Merged labels built any other way remain banned.

   **`inSitu_clusters` is the only clustering whose ids mean the same thing in every
   object** (the vendor assigns them once; `pearson_clusters_*` are re-derived per object
   and so IB1's cluster 3 and RB1's cluster 3 are unrelated). Both samples' audits agree on
   7 of 10 vendor clusters. The three that differ each have a mechanical cause visible in
   the audit's own columns: `f` is Mast cell in all three but capped out in IB1 (14.97% vs
   the 6% cap); `d` is `ambiguous` everywhere with margins of 0.009-0.034; `c` is T cell on
   clean CD3D/CD3E/CD3G in IB1 and merged, and "Neutrophil" in RB1 on a **single gene**.

   **The Neutrophil family has one gene on this panel** (CSF3R; FCGR3B is absent), so it
   cannot reach `confident` — `MIN_FAM_GENES` enforces that. Every Neutrophil call in the
   project rests on one marker. Treat them all as unresolved.

2. **The spatial QC is asymmetric between samples.** IB1 lost 2,013,225 cells to a
   staining/section rule that RB1 was never evaluated for (`qc_spatial_staining ==
   "not_evaluated"`, `fail_staining = 0`). Any cross-sample claim on the `qc_passing_all`
   pools (IB1 400,111 vs RB1 1,399,189) inherits that. A symmetric tier is
   `qc_spatial_fovpx == "ok" & qcCellsPassed`, which gives IB1 1,877,697 / RB1 1,399,189.
   See `output/02_spatial_QC/qc_spatial_attrition.csv`.

3. **The Hodgkin chemokine field is real, has a measured extent, and contaminates
   "background".** Step 13 called the field "tissue-wide"; step 17 read the same IB1-vs-RB1
   background difference as a slide offset. `field_extent_verdict.csv` settles the part that
   is settleable, by asking whether the elevation *decays with distance* from the hotspots —
   a field must, a slide effect cannot. Each sample's OWN top-1% bins are the peaks, so RB1
   is a construction-identical control.

   - **The field is real and Hodgkin-specific.** IB1 falls 2.71x [1.95, 3.91] from near field
     (100-400um) to far field (>1600um); RB1 falls only 1.20x [1.14, 1.26]. Ratio **2.27x
     [1.62, 3.25]**, FOV cluster bootstrap, excludes 1. Most of the decay is spent by
     **400-800um** (75.7 -> 47.6 -> 29.4 -> 12.1 cp10k), then it plateaus.
   - **The field does NOT explain the whole background elevation.** At >1600um IB1 is still
     **4.07x** RB1, flat in distance. Whether that residual is a slide offset or diffuse
     sub-threshold HRS cells is **NOT RESOLVED** — a far-field spatial-coherence test was
     built and is uninformative (IB1 ratio 0.909 p=0.030, RB1 0.944 p=0.007: the no-Hodgkin
     control is *more* clustered, so clustering here is generic).
   - **Do not use all of IB1 `background` as disease-free reference.** 40% of it (2,360 of
     5,898 non-peak bins) sits within 800um of a hotspot and carries field signal. The clean
     reference is `background` **>800um from any hotspot**: 3,538 bins, 60% of it. The binomial
     model is fitted under BOTH definitions (`refdef` column in `binomial_model.csv`); use
     `far_background` for anything referencing IB1 background.
   - **What cleaning it actually did, both directions tested:** every within-IB1 ligand result
     strengthened (CD274 pooled 1.743 -> 1.833, myeloid 1.930 -> 2.014; PDCD1LG2 pooled
     1.810 -> 1.915, B 3.132 -> 3.494), confirming the contaminated reference was diluting
     them. But it did **not** shrink the between-slide offset (CXCL13 in T 5.33 -> 5.55, LAG3
     3.20 -> 2.98), so **the offset is NOT field contamination — that explanation is excluded**.
     And the DiD stayed **0 of 30**: CD274 myeloid reached nominal p = 0.047 and CD274 T
     p = 0.075, i.e. 1 nominal hit against 1.5 expected by chance, nothing surviving BH.
     Cleaning the reference is closed as a route to rescuing the cross-sample comparison.
   - Spearman rho does NOT separate the samples (IB1 -0.180 vs RB1 -0.095 at bin level;
     -0.287 vs -0.280 at FOV level). Both decline monotonically because the field is
     spatially autocorrelated. **The amplitude of the decay is the discriminating quantity,
     not its monotonicity** — a general lesson about rho.

4. **`B_dense` means different biology in the two samples.** The shared architecture axis
   (per-bin B share among typed cells, within-sample quantiles) is the only way to regionalise
   RB1 — its light-chain readout is degenerate, 178/178 neighbourhoods kappa-restricted. But
   RB1's `B_dense` carries the FL germinal-centre programme graded across its boundary
   (BCL2 -0.70, BCL6 -0.70, MME -0.74) while IB1's does not (all flat). IB1's B-dense tissue is
   where the Hodgkin chemokines concentrate instead (CCL17 21.8 -> 7.3 across the boundary,
   flat in RB1). Any `arch`-based result is about B-dense tissue, not about FL.

5. **The between-slide offset is larger than the within-slide regional effects.** Measured
   on the two samples' *reference* zones (IB1 background vs RB1 interfollicular, within T
   cells, depth stratified out), the checkpoint/exhaustion genes differ by odds ratios up to
   **5.33** (CXCL13 5.33, LAG3 3.20, CTLA4 1.88). The largest *within*-sample regional effect
   anywhere in the project is 1.93 (CD274, IB1 hotspot) or 1.82 (PDCD1, RB1 follicles).
   **Check the reference contrast per gene before concluding anything cross-sample** —
   `did_within_celltype.csv` carries `or_reference`/`ref_class` for exactly this. Where the
   reference contrast is large the comparison is foreclosed by the design, not by the
   analysis, and no depth or method fixes it (CXCL13 5.21 disease vs 5.33 reference; LAG3
   4.64 vs 3.20).

   **Do not over-read the 0-of-30 null.** The DiD variance is the sum of two bootstrap
   variances, so its median detectable effect is only **1.8-2.5 fold** — the weakest test in
   the project. Where the reference contrast is *null*, a DiD null is a power failure, not a
   retraction. **CD274 is the live exception**: reference OR 1.11/1.07 (ns, so no offset), DiD
   1.76 in T cells (p = 0.054) and 1.53 in myeloid (p = 0.069). PD-L1 higher in IB1's FL than
   RB1's is unresolved and is the best remaining cross-sample candidate; it wants more FL
   tissue in IB1 (20 FL FOVs vs RB1's 346), not more genes.

   Prefer within-IB1 designs (the CCR4 field gradient is the model) over anything that
   subtracts RB1 from IB1.

6. **Both cell-type callers are neighbourhood-contaminated, and the CLUSTERING is the worse
   of the two.** Holding a cell's own PC-marker count at zero, the cluster calls it PC 30.1%
   of the time within 10um of a plasma cell and 3.6% beyond 80um — an 8.4x gradient on
   proximity alone. Across that range the cluster's PC rate moves **6.8x against argmax's
   2.5x**. The mechanism is specific: the leaked mass is light chains, which
   whole-transcriptome clustering sees and argmax's narrow PC set excludes. **Argmax is
   partly protected because its marker set is narrow** — prefer it wherever a promiscuous
   high-abundance neighbour exists. Step 18 also measured leakage directly: marker
   transcripts are displaced toward the marker-rich neighbour, more where the two callers
   disagree (myeloid +0.127 log2 [0.070, 0.189], T +0.114, fibro +0.103), and those are
   floors because cells where leakage captured BOTH callers sort into the "agrees" column.

7. **A neighbourhood gradient is NOT evidence of leakage on its own.** The step-19e halo test
   (P(label) against a neighbour's marker load, own markers held fixed) cannot separate
   leakage from ordinary tissue architecture, and rising lines are the EXPECTED null: cell
   types cluster, so a cell with no markers beside a marker-rich neighbour genuinely is more
   likely to be that type with dropout. The Plasma cell label was flagged at 12.7x and the
   flag was WRONG — conditioning on the cell's own light chain collapsed it to 1.6x, and
   PC-labelled cells carry half the B markers and 5-10x the light chain of non-PC-labelled
   cells in the SAME donor-load bins. A high fold marks a label to investigate
   (`needs_conditioning`), never one to distrust. `output/19_celltypes/pc_halo_test.csv`.

8. **`cp10k` means two different things.** Step 11 and the bin tables use the *pooled*
   share of a bin's transcriptome; step 10 per-gene uses the *mean of per-cell* shares.
   They are not comparable. cp10k is also compositional and fakes depression in
   top-abundance genes, so results are gated on cp10k **and** detection rate agreeing in
   direction (`scripts/14_script_lightchain_in_out.R:450-460`).

9. **The over-fixation defect is spatially confounded with almost everything in IB1, and
   dropping it per cell does NOT make the surviving area uniform.** Removing the cells is
   necessary and is done everywhere (`16_script_regions.R:166`,
   `16_script_celltype_pseudobulk.R:51`), but it leaves a defect-shaped *sampling* of the
   section. Measured in `output/21_interface/matched/defect_qc_verdict.csv`:

   - **Burden differs by region.** Per-bin defect fraction is **0.33 in FL bins, 0.45 in
     background, 0.59 in cHL hotspots** (IB1 loses 59.1% of its cells overall). Any two-region
     contrast is therefore between differently depleted tissue — SMD −0.60 for FL vs background.
   - **It is a rim, and the rim is the confounder.** rho(defect fraction, distance to tissue
     edge) = **−0.578**; median fraction 0.445 at 100–400 um from the edge and 0.036 beyond
     400 um. Because of that it also rises with distance to the nearest FL focus
     (rho **+0.379**, p = 1.4e-145).
   - **It tracks the chemokine axis.** Inside FL tissue, defect fraction against CCL17+CCL22
     is rho **0.412** at bin level (0.379 at FOV level, SMD 0.85 across a low/high split).
     Consistent with `qc_11_flattening_vs_chemokine.R`, which found the discarded cells
     chemokine-**enriched** (CCL17 1.20x).
   - **Consequence: carry per-bin defect fraction as a covariate**, not just the per-cell drop.
     In step 21 the adjustment moved individual gene deltas by a median 40–50% of their own
     magnitude and flipped 9–20% of signs, without changing the overall verdict.
   - **IB1's "tissue pieces" are coverage, not anatomy.** The same single-linkage decomposition
     gives **2** pieces pre-drop and **3** post-drop — the defect split the lower tissue mass.
     Do not read a piece boundary as a tissue boundary.

10. **The two cell-type callers agree on only 70.2% of cells, and step 16's "RB1 is
   spatially homogeneous" is WRONG — corrected, with the control run.** Step 24 re-ran step
   16's co-occurrence permutation under both label sources; step 24b then ran the control
   that separates better labels from more statistical power
   (`output/24_arch_relabel/matched/`).

   - **Agreement is 70.2%** over 998,544 cells named by both, and collapses for minor
     lineages: B 77.7%, myeloid 63.5%, endo 63.3%, T 59.4%, fibro 57.8%, **PC 37.5%, NK
     32.1%**. Any claim resting on NK or PC identity rests on a coin flip between callers.
   - **The step-19 labels see more spatial structure, and it is NOT a power artefact.**
     Matching every cell type in every region to `min(n16, n19)` — so n_a and n_b are
     identical for every pair by construction — leaves RB1's fraction-at-chance at
     **0.818 (step-16 labels) vs 0.535 (step-19 labels)**, i.e. the gap does not close, it
     slightly widens against the unmatched 0.776/0.503. Threshold-free, the median
     |log2 enrichment| is **0.0594 vs 0.1276** — step-19 labels put RB1's enrichments
     **2.15x further from chance at matched n**. Paired per (region x pair x seed), the same
     pair is further from chance under step-19 labels in **80.3% of 441** comparisons
     (median +0.081 log2, sign test *p* = 2.9e-39).
   - **So step 16's "RB1 sits at chance, 0.93-1.09" is a CORRECTION, not a flag.** RB1 is
     less structured than IB1, but it is not unstructured; step 16 measured the argmax
     caller's silence, not RB1's architecture. Note the matching also reveals why: at matched
     n, step-16's number gets WORSE (0.776 -> 0.818) while step-19's barely moves
     (0.503 -> 0.535), so the step-19 labels are more informative per cell, not merely more
     numerous.
   - **IB1 is where the labels matter least.** Matched, its fraction-at-chance is 0.705 vs
     0.693 — a 0.012 difference — though threshold-free it still moves 0.1751 -> 0.2646.
     IB1's segregation was real under either caller.
   - **Two methodological notes carried from the control.** Step 24 rebuilt the neighbour
     graph separately for each label source, since it subset to labelled cells before
     `nn2()`; step 24b builds it once per region from ALL cells and lets unnamed cells carry
     NA, which moved IB1's own number by ~0.06. And `frac_at_chance` is a poor primary
     statistic: its q comes from a permutation p floored at 1/(N_PERM+1), and BH only passes
     pairs when many pile up at that floor together, so prefer the threshold-free effect size.
   - **Step 16's Tfh findings are NOT re-testable and never will be from clusters.** Step 19
     has no Tfh label by design, so Tfh self-aggregation (2.74x) and the Tfh composition
     result can only ever come from per-cell argmax subsets.

## Closed avenues — do not redo

- **Per-cell HRS/CD30 calling.** TNFRSF8 maxes out at 6 counts in IB1 with 78%
  singletons and depth predicts detection at AUC 0.716. Step 10 established this; step 11
  replaced it with a regional chemokine readout.
- **Light-chain regionalisation of RB1.** Closed on **178/178 neighbourhoods being
  kappa-restricted** (`output/07_lightchain/RB1/lightchain_neighborhoods_r40um.csv`): there is
  no lambda-restricted or polyclonal B region to contrast against, so clonality cannot
  partition the sample. That is the whole of the evidence and it is sufficient.

  **Two claims previously attached to this are WRONG — do not repeat them.** (a) "1,276,791 of
  1,286,223 cells sit inside a kappa-restricted neighbourhood" is a category error: it puts a
  B-cell-subset numerator over the whole-sample denominator. Only **466,764** RB1 cells were
  ever light-chain-evaluable, of which **454,370** have a neighbourhood and **12,394** do not.
  So "the clone *is* the node" is overstated ~3x. (b) Step 14 **never ran on RB1** —
  `output/14_lightchain_in_out/` contains `IB1/` only, and the `coverage_verdict.csv` this file
  used to cite **does not exist**. The 12,394 unassigned evaluable cells have never been
  examined; if they are spatially clustered they would be RB1's only internal reference.
- **Compartment-restricted cell typing.** Nuclear-only argmax scores **-3.53 pp** against a
  marker-matched random control; non-nuclear scores -1.57. `CellComp` is assigned relative to
  the segmentation mask, so a transcript that bled into the wrong cell and landed over that
  cell's nucleus is called `Nuclear` — compartment tracks position inside a possibly-wrong
  boundary, never provenance. Markers are also the cytoplasm-enriched fraction (64.0% nuclear
  vs 78.4% overall), so nuclei are the marker-POOR compartment.
- **Per-cell leakage flags of any kind.** Every candidate score is at chance for flagging a
  miscalled cell: positional displacement AUC 0.446-0.526, donor marker load 0.354-0.605,
  nucleus-vs-cytoplasm discrepancy 0.383-0.554. The myeloid effect is +0.127 log2 across
  5,913 cells and is invisible in any one of them. **Leakage is real at the population level
  and unmeasurable at the single-cell level** — correct it per cluster or not at all. The
  compartment version specifically measures marker SPARSITY: its relationship to positional
  displacement is a symmetric U (myeloid 0.329 / 0.104 / 0.343 across the z range), and where
  the slope is significant it runs the wrong way (B -0.237, p = 1e-6).
- **Imputation before projecting cells into a reference embedding.** 68% of the cells outside
  step 19's high-confidence set are over-fixation defect, and defect cells are NOT sparse —
  median nCount 623, deeper than the non-defect cells that project successfully at 475 — with
  label agreement **flat at 0.26-0.29 across a 7x depth span**. That is corruption, not
  sparsity; imputation fabricates plausible values from a biased observation. The two standard
  forms are also blocked here: expression-space imputation is circular with the projection it
  would feed, and spatial imputation manufactures the halo artefact by construction.
- **9p24.1 copy-number inference from this panel.** Attempted properly in step 23 and it
  FAILS ITS OWN NEGATIVE CONTROL. The panel carries **nine** 9p24.1 genes, five not
  immune-functional (`CDC37L1`, `ERMP1`, `KDM4C`, `MLANA`, `PLPP6`), so gain (raises all) is
  separable from interferon induction (raises only `CD274`/`PDCD1LG2`/`JAK2`); the null is the
  other 440+ cytobands on the same panel. IB1's hotspots DO score in the **99th percentile**
  (rank 5/439 in T, 5/443 in fibro) with the non-immune subset rising more than the immune
  one. But **RB1's fibroblasts also fire** (0.195, q = 0.039) and RB1 has no Hodgkin
  component; and the signal is spread across three of four lineages when a clonal event must
  be **confined to one**. Both say regional expression gradient, not dosage. **Step 17's
  9p24.1 module stays EXPLORATORY — step 23 is the reason not to upgrade it.**
  `output/23_dosage/matched/dosage_9p24.csv` (its `verdict` column is control-aware; an
  earlier version did not consult RB1 and labelled RB1's own fibroblasts "DOSAGE-CONSISTENT").
- **Reading a ligand-receptor proximity enrichment without looking at the cells.** Step 22's
  snapshots (`output/22_lr/matched/snapshots/figures/`, 16 figures, 8 pairs x
  {cherry, repr}) crop real cell pairs onto real morphology with every segmentation boundary
  and both genes' transcripts drawn. Two lessons visible only there: the receptor side is at
  the 2-count positivity floor in 8 of 9 representative panels even for the strongest pair
  (CXCL13-CXCR5), so a "pair" is typically a well-expressed ligand cell beside a barely
  positive receptor cell; and PDCD1LG2-PDCD1's entire candidate pool is **10 pairs in 6
  FOVs, every one 2-vs-2 counts**, with representative cell-type combinations (B-B,
  Endothelial-Stromal, Myeloid-NK) that are not the myeloid/DC-to-T axis PD-L2 signals on.
  **`Hodgkin chemokine niche` must be excluded from every cell-type role** — it is a region
  label; leaving it in put a region in 16 panels as a lineage and supplied 60% of
  CCL22-CCR4's candidate pool (50 -> 20 pairs once removed).
- **Reading a ligand-receptor proximity enrichment without its distance profile.** Step 22
  measured co-expression in proximity for 2,315 curated pairs (`scDiffCom` LRI_human,
  both partners on the panel) across five distance bands. Effects are **1.01–1.18 fold** and
  only reach significance where FOV counts run to the hundreds — IB1's hotspots (16 FOVs) and
  FL foci (25) yield nothing. What survives: `CXCL13`-`CXCR5` enriched across **all four**
  bands in IB1 background (1.14/1.16/1.15/1.13), then `CCL19`/`CCL21`-`CCR7` at 1.01–1.07;
  `CCL17`/`CCL22`-`CCR4` only 1.055–1.075 and **below 1 inside the FL foci**, so the CCR4
  gradient is regional, not pairwise. **The largest value in the whole table is an artefact
  and the design caught it**: `PDCD1LG2`-`PDCD1` in IB1 FL foci runs
  **2.28 → 1.54 → 1.14 → 1.01** across the bands (q = 0.167) — a monotone decay to nothing by
  40–80 um, which is the leakage signature step 18 measured. Never quote a single band.
- **Distance-to-hotspot as the Hodgkin field-exposure axis.** REJECTED on its own positive
  control: CCL17+CCL22 inside FL tissue is 23.3 cp10k in the near half against 26.5 in the far
  half, and exposure is uncorrelated with distance across the 22 FL FOVs (rho −0.126, p = 0.575).
  The hotspots are IB1's **top 1%** of chemokine bins and do not enumerate the field's sources —
  step 17 flagged this — and only **47 of the 80** survive the low-information filter. FL FOVs
  693/694 sit ~3,000 um from any hotspot at 37 and 77 cp10k. **Use the bin's own CCL17+CCL22
  level as the exposure axis**, which is what `16_script_field_gradient.R` already did.
  `output/21_interface/matched/gate_verdict.csv`.
- **The FL response to the Hodgkin field, as a within-IB1 interaction.** Built properly (step 21:
  FL slope minus background slope on the field axis, within cell type, FOV cluster bootstrap,
  depth and defect burden as covariates) and it is a **bounded null** — 0 of 44 pre-declared
  genes in each of B/T/myeloid, and 0 of 6 cell types on composition. Median resolvable effect
  **1.09–1.14 fold** on 20–22 FL FOVs. This is a power limit, not a demonstration that the FL is
  indifferent: it wants **more FL tissue**, the same remedy CD274 needs. The genome-wide tier
  returns 15 genes which are **not interpreted** — on the detection floor, no coherent programme,
  and *CD3D* called down in B cells while *CD3E* and *TRAC* do not move.
- **"No spatial relationship between the two diseases" as a CLEAN null.** The verdict survives
  but the claim of cleanliness does not, and the reason is fact 9. Step 13 rested on one
  bin-level Wilcoxon over 11 foci (p = 0.233); with the FOV as the unit it is p = 0.332, so
  pseudoreplication was **not** the defect. Under a null that relocates the Hodgkin foci
  preserving their number, size and clustering, observed 949 um against a null median 1105
  (p = 0.53). Under the same null **matched on defect burden** — which the real hotspots
  have and average tissue does not — the null median is 1503 um and p = 0.152 (0.085 in the
  reverse direction). Nothing reaches significance either way, so **arrangement is UNRESOLVED**,
  not negative. A label-permutation null gives a nominal p = 0.045 and is rejected: it inflates
  within-pattern spacing 1.37x, i.e. it destroys the clustering it should hold fixed.
  `output/21_interface/matched/arrangement_verdict.csv`.
- **Cross-sample neighbourhood matching for B cells.** Not identifiable: IB1's most B-dense
  quartile is less B-rich than RB1's least B-dense quartile, so 2 of 89 bins find a partner.

  **The IB1-vs-RB1 B-marker difference it was meant to test is RESOLVED, and it is mostly an
  artefact of the label.** IB1's B-cell label is only 14.8% argmax-confirmed at a median 221
  counts, against RB1's 51.3% at 658 — so a pooled cp10k over B-labelled cells compares two
  different populations, which is what makes IB1's column pale in the by-sample heatmap.
  Restricted to argmax-confirmed cells IB1's share is HIGHER (73.6 vs 68.3); depth-matched
  within nCount deciles it settles at IB1 ~10-15% lower, consistently. **Dilution was tested
  within IB1 and REFUTED**: the absolute B-marker count falls with a less-B neighbourhood
  (slope -1.11 [-1.41, -0.77]) while the share does not move (-2.80 [-8.12, 2.12]), the
  opposite of the dilution signature, and the count decline just tracks total depth. Any
  cross-sample marker comparison must use argmax-confirmed cells within depth strata; see
  `output/19_celltypes/figures/README_marker_heatmap_versions.md`.
- **Cluster-based T-cell subsetting (Tfh / Treg / CD8 / CD4).** Tested three ways in step 19
  and it fails every time. Sub-clustering the 142,361 T cells jointly separated SAMPLE, not
  subset (4.6%-82.7% IB1 against a 26.2% baseline) while every sub-cluster had CD4_conv as
  its modal argmax subset at enrichment 0.75-1.26. Removing that confound by sub-clustering
  each sample on its own changed nothing: 12 sub-clusters across two independent runs, all
  CD4_conv, enrichment 0.89-1.14 (IB1) and 0.67-1.23 (RB1). The input is too thin -- only
  **10-23% of T-labelled cells carry an argmax subset call** and CD4_conv is ~50% of those.
  At ~500 counts/cell the separating genes (FOXP3 has 3, CD8A/GZMK/GZMB, CXCL13/ICOS/CXCR5)
  cannot drive a neighbour graph that whole-transcriptome variation dominates. **T stays a
  single label; Tfh comes from the per-cell argmax subsets, which is what step 16 already
  does.** `output/19_celltypes/t_subcluster_labels_{IB1,RB1}.csv`.
- **Per-cell CD45 protein intensity as a CD45-negative (HRS) marker.** Step 30
  (`scripts/30_script_cd45_rosette.R`) built a "CD45-dim cell among CD45-bright neighbours"
  score from the vendor's `Mean.CD45`, behind a pre-declared gate, and the GATE FAILED: the
  per-cell intensity does not separate CD45+ from CD45- lineages -- T/B over fibro/endo AUC
  **0.52 in IB1 and 0.47 in RB1**, `Max.CD45` no better (0.52 / 0.47), B cells dimmer than
  fibroblasts, and only Spearman **0.08** with PTPRC RNA. RB1 fails as badly, so it is the
  per-cell measurement, not IB1's staining defect: in packed lymphoid tissue CD45 sits on
  shared membranes and a segmented cell's intensity is mostly its neighbours'. The gene
  results it produced are NOT interpretable and must not be quoted (the score marks cells
  with FEWER lymphocyte transcripts -- CD3D 0.41x, PTPRC 0.71x top vs bottom quintile).
  `output/30_cd45_rosette/gate_verdict.csv`, `gate_diagnostics.csv`. The same caution
  applies to `Mean.Membrane` and every other per-cell stain mean here; a CD45-negative cell
  can be judged by eye on the step-29 `_gray.png` panels, not by these numbers.
- **HCA tonsil label transfer.** Disabled 2026-08-06; its labels anti-track marker
  evidence. `celltype_coarse`/`celltype_fine` are not carried past step 03.
- **Transcript-level proximity of PDCD1 to CD274.** Three null designs built, all three
  failed their own positive control (centroid/transcript mismatch; lost over-dispersion;
  carrier-cell resampling). Withdrawn per a pre-committed stopping rule.
  `output/17_pd1/matched/tx_proximity.csv` exists and **is not a result** — see the
  `tx_proximity_NOT_DELIVERED.md` beside it. Spatial engagement of the axis is *unmeasured*,
  not measured-and-null.
- **Sign consistency across T compartments.** RETRACTED, and the primary reason is
  pseudoreplication, not the slide offset: the 5 T compartments have **ICC 0.875**, i.e.
  **1.11 effective observations per gene**, so the null p = 2*0.5^5 = 0.0625 is wrong. All 14
  "consistent" genes have a compartment mean exceeding their compartment SD, so agreement is
  arithmetic. Under the correct null (p = 0.5) the same 14/22 gives **p = 0.143, not
  significant** (`sign_pseudoreplication.csv`). Secondary reason: 7 of 9 checkable genes point
  the same way between the *reference* zones (`sign_offset_check.csv`). Only PDCD1 and CXCL9
  survive either check. **This is the region-as-unit rule broken with "compartment"
  substituted for "cell"** — when inventing a new consistency statistic, first ask what a
  uniform group-level shift does to it, and count the effective replicates.

## Reports

**Cross-project reviews live in `reports/`, not in the numbered `output/NN_*` namespace.**
`reports/manuscript_plan.html` (+ `MANUSCRIPT_PLAN.md`) is the manuscript-readiness review:
what is established, the remaining gaps, a draft abstract and a figure list. **Revised
2026-09-08.** Four of the original eight gaps are closed (G1 by steps 24+24b, G4 by step 22,
G5 by step 23, G7 by step 21) and are **no longer listed** -- only a summary table of them
survives. **Gap numbering is preserved**, so the open gaps keep their original identifiers:
**G2, G3, G6, G8**, with G2 (no morphological or protein confirmation) the top priority. Do not
renumber them; earlier drafts and this file cross-reference the originals. The review also
carries R9 in its CORRECTED form and a draft abstract revised to match -- the claim that RB1
had *lost* spatial organisation is withdrawn. Its figure list carries **embedded previews** of every panel:
`scripts/26_script_embed_review_figures.py` downscales each referenced PNG and inlines it as a
`data:` URI, so the page stays self-contained (10.5 MB, no external refs). It is **re-runnable**
-- it strips the previous embed before inserting the current one, so run it again after
rebuilding any panel. Two things it has to get right and does: it excludes non-primary variants
(`preqc`, `_prev`, `SUPERSEDED`) so a pre-QC map cannot appear unlabelled beside the matched-tier
one it must not be compared with, and it prints the full-resolution path under every thumbnail
so the original stays findable. Per-step reports stay with their step, in
`output/NN_summary/`. Note the global rule says `report/` singular; this project uses
`reports/` -- do not create both.


`output/17_summary/composite_node_lexicon.html` is the companion glossary — the immunology
and the statistics, terms tagged field-standard vs project-specific. **This file is the only
copy; edit it here.** It was briefly also published as a hosted artifact, which is abandoned
and stale: something republished it from outside the session and reverted a correction, so the
hosted page still teaches the retracted sign-consistency result as sound method. Do not link
to it and do not re-publish this file — see the house rule below.

**Reports are local HTML files in this repo, never hosted artifacts.** That rule exists partly
because of the above: a hosted copy acquires its own edit history that the repo cannot see, and
the two silently diverge.

Every summary report includes `scripts/lib/methods_notes.Rmd` as a knitr **child**, not a
link — collaborators get a single HTML file with no repo, so the methods text must be
compiled in. One source of truth, inlined into each report. Reports also set
`mathjax: null` and use plain markup for formulae: `self_contained: true` still fetches
MathJax from a CDN, which breaks every equation for a reader without internet.

## Manuscript figure panels

`scripts/25_script_figure_panels.R` builds every panel the review tagged `[make]` into
`output/25_figures/` (15 PNGs). It is **plotting only** -- each panel reads a saved result
table and recomputes nothing, so a panel cannot drift from the result it depicts.
`scripts/25_script_composition_fov.R` is the one piece of new analysis, and it closed **G8**.

Three things this step established that are easy to walk back into:

- **`lr_proximity_profile.csv` and `lr_proximity_tier1.csv` hold DIFFERENT enrichments.** The
  profile's `e_*` columns are the analytic estimate over ALL band pairs including cross-FOV
  ones; `enrich_fov` in the tier-1 table is within-FOV only and is what the bootstrap
  resamples and what step 22, the review and this file all quote. Plotting the analytic column
  gave PDCD1LG2-PDCD1 as NA/1.65/1.13/1.03 instead of 2.28/1.54/1.14/1.01 -- losing the 0-10
  band to the expected-count floor, i.e. losing the entire point of the panel. **Use
  `enrich_fov`.**
- **Always give explicit `breaks` on a log axis.** `scale_x_continuous(trans = "log2")` with
  default breaks rendered the odds-ratio forests as `7.6293954.8628123.045...` -- run together
  and unreadable. `OR_BREAKS <- c(0.25, 0.5, 1, 2, 4, 8)` in that script is the fix.
- **`module_contrasts.csv` mixes within-sample and CROSS-SAMPLE contrasts.** Two of its five
  are IB1-vs-RB1, which fact 5 forecloses for the offset-affected genes. Any figure drawing
  from it must mark them; the step-25 panel uses a triangle and says so in the caption.

**G8 is closed, and it qualifies a phrasing.** IB1's composition is much less stable across
the section than RB1's: its two largest fractions have almost completely overlapping IQRs
(myeloid 17.7-34.1% vs B cell 17.5-35.2%), so **"IB1 is myeloid-dominant" is not supportable**
-- though "IB1 23.6% B vs RB1 71.4%" is, RB1's B IQR being 66.1-77.7% and disjoint. The pooled
myeloid estimate also sits 4.4 points above the FOV median, so a few dense fields carry it.
`output/25_figures/composition_fov_summary.csv`.

## Chemokine source cells (steps 27-29)

Who makes CCL17/CCL22: HRS cells or dendritic cells. `scripts/27_script_hotspot_sources.R`
(the 80 hotspots vs RB1's own top-1% bins) and `scripts/28_script_bright_sources.R` (per
cell, sample-wide, by brightness tier) are the analysis; step 29 is the per-cell list, the
figures and the ROI review set. Outputs in `output/27_hotspot_sources/`,
`output/28_bright_sources/` and `output/29_bright_cells/`.

- **The hotspots are a DENSITY detector and miss most bright cells.** Pooled bin cp10k
  dilutes an isolated bright cell (one 100-count cell adds ~52 cp10k to a ~58-cell bin,
  under half the cut), so hotspots are clusters of moderate cells. Only **14% of IB1 cells
  with >= 50 chemokine counts (237 of 1,693) and 18% of those >= 100 (70 of 392) sit in a
  hotspot**; 38 of 80 hotspots contain no >= 100 cell. Step 11's 50/100/200 um sweep kept
  the same pooled rule and so never tested this. For anything about the sources, call
  cells, not bins.
- **Bright cells are Hodgkin-specific; moderate ones are not.** >= 50 counts: IB1 1,693 of
  767k cells vs RB1 **9** of 1.29M; >= 100: IB1 392, RB1 **0**. Cells with 10-49 counts
  exist in both (IB1 ~1.9%, RB1 ~0.1%).
- **Brighter = more HRS-leaning, not more DC.** In IB1, against depth-matched < 10 cells
  (depth EXCLUDES the chemokine counts), BATF3 8.5x -> 25x, IL13 1.9 -> 4.0, TNFRSF8
  1.5 -> 2.9, IRF4 2.3 -> 5.8 from the 10-19 to the 100+ tier (100+ vs 10-19: LTA 3.4,
  BATF3 3.0, IL13 2.1, TNFRSF8 1.9, all q < 1e-10), while DC-lineage genes do not rise
  (LAMP3 1.12 -> 0.74, XCR1 0.88 -> 0.60, CD1C and IDO1 flat). CCR7/CD80/CD83/CD40/CD86 rise
  too -- an activated-APC programme HRS cells also carry, so not discriminating. Area rises
  58 -> 108 um^2 but 10-49 sources are the same size in RB1, so size is partly generic.
- **IL13, BATF3, TNFRSF8 co-occur in the same IB1 sources beyond a matched null**: triple
  O/E **1.36** [1.29, 1.56] vs random detection-matched triples median 1.09 (p = 0.005);
  pairs 1.09-1.18. Lineage controls give 1.47 (CD79A-MS4A1) and 1.57 (CD3D-CD3E), so the
  triple is of lineage-like strength. **CCR7-LAMP3 does NOT co-occur (0.97)** -- no mregDC
  subpopulation among the sources. Two limits: BATF3's matched null pool is only 13 genes
  (its 52% detection among sources is unusual), and O/E shrinks under high detection, so
  do not compare O/E across tiers -- use the k-of-3 shares (>= 2 of 3: 10.6% at 10-19,
  54.6% at 100+; RB1 4.8% at 10-19). RB1's co-occurrence test is underpowered (its lineage
  controls barely clear their null), so it is not a usable contrast.
- **Step 27's pooled source excess is a population average over 23% of hotspot cells**,
  too many to be HRS; read it as "IB1 sources carry an HRS-compatible excess RB1's lack"
  (IL13 4.5x, BATF3 3.5x, CCR7 2.8x IB1/RB1), not as a cell identity.
- **Still NOT established: that any cell is an HRS cell.** Leakage localises a programme to
  a cell OR its contact (fact 6); FUT4 is flat (CD15 is a protein glycan, weak as RNA).
- **A k-of-3 gate is NOT a usable HRS filter; brightness is.** Among >= 50 cells the
  IL13/BATF3/TNFRSF8 triple is 123 observed vs ~108 expected under independence, and
  ">= 2 of 3" is only ~2x more common in IB1 than RB1 at matched brightness, so ~45% of
  passers would pass without HRS biology. Brightness >= 50 is ~300x Hodgkin-specific.
  Rank within the bright cells by the step-29 score; never gate on the genes.
- **Bright cells are NOT excluded from FL neighbourhoods.** Raw rates look enriched near FL
  (>= 50: 2.4x inside the 40 um footprint) but FL tissue is deeper (median depth0 ~500 vs
  ~290 beyond 800 um) and a count cut is easier to reach in deep cells. Depth-standardised,
  no FL-distance band excludes 1 (>= 50 in footprint 1.61 [0.83, 2.40]). On the control
  axis (ANY light-chain neighbourhood) >= 10 sources run 0.65-0.83x at 40-400 um and 1.59x
  beyond 800 um -- a lean away from B-cell clusters in general, not from FL, and >= 50
  cells do not share it. `output/29_bright_cells/fl_exclusion.csv`.

- **Radial CD45 (step 31) -- candidates are CD45-dim at the centre, but NOT rosetted.**
  Segmentation-free: mean CD45 in 5 um rings to 100 um around 500 random >= 50 cells and 500
  area-matched < 10 controls, read from the full-resolution image, each profile over its own
  50-100 um mean (`scripts/31_submit_radial_cd45.sh`). Centre (0-5 um) log2 **-0.44
  [-0.54, -0.34]** vs controls **-0.08 [-0.14, -0.02]** (FOV bootstrap), but no bright ring:
  10-20 um sits at -0.03 vs +0.02. Sorting both groups onto the candidate k-means shapes, the
  dim-centre-bright-ring shape is **equally common** (20.7% vs 19.3%), dim-centre-no-ring
  shapes are 1.8x enriched and the CD45-bright-centre shape 0.55x. And the HRS-leaning score
  does NOT track dimness (median 0.45-0.56 in the dim shapes vs 0.83-0.87 in the ring and
  bright-centre shapes), plausibly because IL13/TNFRSF8 leak in from activated T neighbours.
  **The classes are not discrete**: only k = 5 of k = 2-6 beat the single-Gaussian null, by
  0.001 in silhouette (0.177 vs 0.176), with silhouettes ~0.18 throughout -- a smear sliced
  by k-means. The rule as written had no correction for trying five k; next time require
  the pass at the 99th pct or across adjacent k. `output/31_radial_cd45/`.
- **The step-31 CD45 shapes do NOT separate the RNA (step 33).** Per cluster, the
  controls nearest each centroid reproduce the shape almost exactly -- including cluster 1
  (dim centre, bright ring at 10-20 um; 19% of controls) and cluster 5 (deeply dim centre,
  ring at 25-35 um) -- so the shapes are things ordinary cells do; only their frequency
  differs (clusters 4/5 1.8x, cluster 3 0.55x). Clusters 4/5 are simply LARGE cells (area
  109 um^2 vs 85-94; the DNA dip at 5-10 um is deeper), and cluster 1 is any cell with
  CD45-bright neighbours. In the step-28 gene framework (obs/exp vs depth-matched < 10
  cells, FOV bootstrap), **no HRS-leaning, DC or lineage gene is heterogeneous across the
  five clusters** (label-permutation test, BH q 0.45-0.94): IL13 2.8-4.1x, TNFRSF8 1.7-2.5x,
  LTA 2.3-4.2x everywhere. The one hit is **BATF3** (q = 0.036): 22-25x in the CD45-RICH
  contexts (clusters 1, 3) vs 12.6-15.7x in the dim ones (4, 5) -- the opposite of what an
  HRS-by-morphology reading predicts, and not T leakage (CD3E is < 1 in every cluster and
  CD3D is highest in cluster 4). Power: 63-140 cells per cluster, intervals ~2-fold wide, so
  this excludes large RNA differences only; it was decisive for BATF3 because BATF3 is the
  best-detected gene. `scripts/33_script_cluster_followup.R`, `output/33_cluster_followup/`.
- **Wide CD45 (step 32): candidates sit in a CD45-dim zone >= 400 um across, and it is
  mostly LYMPHOID COMPOSITION, not a CD45-negative niche.** Same cells, 25 um bands to 500 um
  at pyramid level 3, each profile over its own 400-500 um ring
  (`scripts/32_submit_wide_cd45.sh`). Candidate minus control: **-0.33 log2 at 0-25 um,
  -0.24 at 100 um, -0.10 at 250 um, 0 by ~400** (FOV bootstrap excludes 0 to 350 um); DNA
  is slightly HIGHER around candidates (+0.05 to +0.08), so it is not sparse tissue; the
  defect fraction is identical in both groups (~0.46) and adjusting for it changes nothing.
  What explains it: the zone is B-rich (19-23% vs 11%), T-rich (19-26% vs 11%) and
  myeloid-poor (24% vs 31%), decaying over the same ~400 um, and B cells are the dimmest
  CD45 lineage on this stain (step 30: B 1,209 vs T 1,508 vs myeloid 1,789). Adjusting for
  the band's B share cuts the deficit to **-0.14 at 25-50 um and -0.08 at 75-250 um**,
  intervals grazing zero; only the centre band keeps -0.26 (the cell itself, step 31).
  Consistent with fact 4 (chemokines concentrate in B-dense tissue). Two limits: the own-
  ring reference forces zero at 400-500, so ">= 400 um" is a floor not a size; and 100 um
  is plenty for the rosette question (all centre structure is inside 15 um) but any
  environment question needs this wider design. `output/32_wide_cd45/wide_profile.csv`.

- **Nuclear size / shape (step 34): the candidates' nuclei are PALE-FILLED, not larger, and
  the RNA runs the other way.** `scripts/34_script_nuclear_size.R` (vendor NucArea /
  NucAspectRatio / Mean.DAPI, population-wide) + `34_extract_polygons.R` +
  `34_nuclear_in_polygon.py` (DNA measured inside each cell's vendor polygon, the 1,000
  step-31 cells) + `34_polygon_qa_sheet.py` (look at them). `output/34_nuclear_size/`.
  - **Vendor NucArea rises with brightness** (IB1 30.7 -> 50 -> 56 -> 61 -> 70 um^2 across
    tiers) and is 1.2-1.45x larger in >= 50 cells than < 10 cells of the same cell-area
    decile -- but the IMAGE, area-matched, finds no dense-nucleus difference (42.7 vs 40.0
    um^2, ns). The vendor's nuclear segmentation is staining-dependent (nuclear fraction
    0.78 in RB1 vs 0.54 in IB1 at < 10; Spearman with the image area only 0.57), so do not
    read its within-IB1 contrasts as nuclear size.
  - **What IS different, from the picture:** candidate polygons are filled with pale,
    uniform DNA at ~0.7 of dense-nucleus intensity -- pale fraction 0.82 vs 0.72, DNA CV
    0.26 vs 0.33, "mostly pale, not dense" 33% vs 13% of controls [diff 0.15, 0.25] -- with
    fewer dense pixels at the centroid (54% dim vs 36%) and more bright nuclei pressed
    against the outside (0.46 vs 0.39). Lobation (DNA maxima) and nucleoli (dark holes) do
    not differ. **By eye the pale-filled polygons cannot be told from focus blur** in a
    single plane (`qa_polygons_*_pale.png`); vesicular HRS chromatin would look the same.
  - **The HRS-leaning RNA goes with dense, bright nuclei, not pale ones.** Within brightness
    band, brightest vs palest Mean.DAPI quintile: IRF4 1.8-2.1x, BATF3 1.4-1.6x, IL13
    1.3-1.7x (q < 0.01 at 10-49), LYZ 0.36-0.58x; largest vs smallest NucArea quintile:
    LTA 0.5-0.6x, CCR7 0.6-0.7x, IRF4/BATF3/CD40/CD83 0.6-0.8x, nothing up. Among the 500
    candidates, pale-fill vs HRS-leaning score rho = -0.19. Third axis in a row (CD45 shape,
    step 33; nuclear size; DAPI paleness) where morphology-HRS-like and RNA-HRS-like cells
    are different cells. Either the bright population is two populations (pale large
    nuclei, RNA-quiet; dense activated nuclei carrying IRF4/BATF3/IL13) or the pale ones
    are out of focus and lose in-plane transcripts; the data here cannot say which.
  - **Two gate revisions are recorded in the script header.** The first pre-declared gate
    (lymphocyte nuclei smaller than myeloid/fibro; NucArea vs the 5-10 um DNA band) was
    wrong on both premises and was replaced by image agreement BEFORE the full run; and a
    centroid-anchored image measure was abandoned because candidate centroids sit on dim
    pixels (725/1,000 flooded). Polygon-bounded measurement is the usable design.
  - Image reads from scratch are I/O-bound and variable (0.2-2 s per 40 um window).

- **Fragment reconstruction (step 35): chemokine-bright polygons are CONTIGUOUS far beyond
  chance, but the composites are mostly multi-nucleated groups, not single carved-up giant
  cells.** `scripts/35_submit_giant_cells.sh` (select -> per-FOV polygon rasterisation and
  adjacency as a SLURM array -> analyze; ~1 h). Two partitions: connected components of
  touching positive polygons (contact >= 3 px) and a nucleus-anchored assignment (one dense
  nucleus per composite); null = 200 within-FOV permutations of the expression rows.
  - At >= 20 counts **60% of IB1 seed composites have >= 2 members vs 6% under the null,
    20% exceed 300 um^2 vs 0.1%** (p at the 0.005 floor); anchored 20% vs 1.5%; holds at
    >= 50 (22% vs 1.4%) and without seam composites. RB1's brightest cluster only at >= 10
    (19% vs 1.3%), not at >= 20 (1.7%).
  - Nuclei: multi-member composites hold FEWER dense nuclei per member than random
    contiguous sets (0.75 vs 0.84 at >= 20, 0.71 vs 0.85 at >= 50) but 44% contain >= 2
    nuclei and the "<= 1 nucleus" share is not elevated (26% vs 28%). The typical composite
    is several nucleated bright cells plus some nucleus-free fragments.
  - RNA within composites (anchored, 312 composites): fragments carry MORE programme than
    the nucleus member (CCR7 1.84x [1.51, 2.27], BATF3 1.25x, LTA 1.25x, TNFRSF8 1.34x ns)
    and LESS T transcript (CD3E 0.45x, CD3D 0.64x) -- the predicted direction, modest size,
    same under the split/dropped ambiguity variants (8% of fragments ambiguous).
  - **"Nucleus-free" means no DENSE nucleus** (DNA >= the field's Otsu, set by packed small
    lymphocytes). The panels (`panels_IB1_thr20_anch.png`) show many fragments enclosing a
    PALE nucleus, so composites are groups of chemokine-bright pale-nucleus cells around one
    dense nucleus -- step 34's pale phenotype, spatially grouped -- not cytoplasm carved off a
    nucleus. Read "fragment" accordingly. **The panels also draw every CCL17/CCL22 molecule**
    (from step 29's `ib1_chemokine_tx.csv.gz`, IB1 only) with the step-29 density rings and a
    bar per touching pair carrying step 36's class. cid 3505 is the anchor rule failing on
    screen: 787 molecules in three pale-nucleus polygons attached to a 21-count dense-nucleus
    bystander. **Do not judge the mixture from the panel sheet** -- it is the 24 LARGEST
    composites by area, which favours chains, and by eye it read as "mostly chains"; step 36
    scored every pair and found the opposite.

- **One cloud or two (step 36): ~75% of touching member pairs are ONE molecule cloud across
  the polygon boundary, so the typical composite is one source carved into polygons.**
  `scripts/36_script_molecule_dip.R` (~1 min on 8 cores). Statistic: members' molecules
  projected on the centroid axis, 1 um KDE, density at the polygon boundary (from the polygon
  support functions) over the smaller side's peak; cut 0.74 = equal-error point between two
  count-matched references from isolated bright polygons (one cell cut by a random line vs two
  cells placed polygon-tangent), each misclassified ~10% (3-7% once the dimmer member has
  >= 35 counts). Real pairs are bimodal on exactly those two modes. Anchored, >= 20, seed
  composites: 78% of pairs one cloud; 75% of composites one cloud on every edge, 13% separate
  on every edge, 12% mixed; components chain more (61% / 7% / 32%). Flat across thresholds and
  count bins. Nucleus-to-fragment edges are MORE continuous (79%) than fragment-to-fragment
  (65%). **Two earlier statistics failed their own references and must not be revived**: a
  density sampled along the centroid LINE (sigma 2 um merges clouds 10 um apart; sigma 1 um
  misses any cloud off the line) and an unsupervised mixtools EM (spurious optima with a
  near-empty component in a third of the two-cell reference). "One cloud" cannot separate a
  carved cell from a dim neighbour holding only spill-over -- both mean the fragment is not an
  independent source. `output/36_molecule_dip/{cut.csv, summary.csv, pair_class.csv,
  composite_class.csv, pairs_by_nucleus.csv, panel_sheet_pairs.csv}`.
  - **Design gap found, not fixed:** anchors had to be chemokine-positive, so 69% of
    nucleus-free fragments at >= 20 are ORPHANS touching no positive nucleated polygon. A
    giant cell's nucleus polygon has little cytoplasm and may not clear the threshold. The
    follow-up is "any nucleus-bearing polygon can anchor", with step 34's pale-nucleus
    measures on those anchors. `output/35_giant_cells/`; `35_composite_panels.py` draws the
    largest composites on morphology; report: `output/35_summary/hrs_candidates.html`
    (`scripts/35_summary_report.Rmd`, rendered by `scripts/sbatch/render_35_summary.R` after
    `scripts/35_report_figures.py`).

### Step 29: per-cell list, figures, ROI review set

- `source_cells.csv.gz`: every source (>= 10) in both samples with canvas px, tier, the
  HRS-leaning score and k-of-3. **The score** is the mean over IL13/BATF3/TNFRSF8 of
  log2((obs + 0.5) / (exp + 0.5)), exp from the same sample's < 10 cells in the same depth0
  decile -- step 28's reference. It ranks; it does not call.
- `roi_png/<group>/`: 120 full-resolution 150 um PNGs (+ a contact sheet per group) in six
  groups of 20 -- high (100+), medium (50-99), low (10-49) x {random draw, top score} -- in
  the step-15 style (imports `15_roi_pngs.py`: same ranges, Membrane x0.5, G omitted).
  **Rings are CCL17+CCL22 density contours** (sigma 8 um) at 6 / 12 / 24 molecules per
  100 um^2, roughly a 25 / 50 / 100-count cell's own peak, so ring count reads as
  brightness. `roi_png_index.csv` has `own_tx_matches`: the selected cell's assigned
  molecules must equal its chemo count (all 120 do). This is the G2 review set. Each
  composite has a `_gray.png` twin (sorts straight after it): DNA / CD45 / Membrane / PanCK
  as single grey panels over the same fixed ranges, **no** Membrane gain, same rings and
  crosshair, and a "what to look for" caption line for the pathologist.
- `ib1_chemokine_tx.csv.gz`: all 1,272,591 IB1 CCL17/CCL22 molecules (54,945 unassigned) in
  canvas px from the flat tx file via `cosmx_local_to_px`; agrees with the per-FOV transcript
  files (4,921 vs 4,918 in a test box).
- **Two traps from building it.** (a) The smoke run reads the tx file's first lines, which
  cover only FOVs 1-4 and cut the last one part-way; smoke candidates are therefore drawn
  from the fully read FOVs, or their windows carry no molecules and the ring code never runs
  -- which is how the first smoke PNGs came out ringless. (b) The full R run takes ~40 min,
  almost all of it the awk pass over the 13 GB tx file.
- Regenerate: `J=$(sbatch --parsable --mem=48G --time=1:30:00 scripts/sbatch/run_rscript.sh
  scripts/29_script_bright_cells.R)` then `sbatch --dependency=afterok:$J --mem=16G
  scripts/sbatch/run_python.sh scripts/29_bright_cell_rois.py` (~9 min). Both take `smoke`.

### Step 37: the hotspot blind spot on the image

`scripts/37_hotspot_zoom.py` (`sbatch --mem=16G scripts/sbatch/run_python.sh ...`, 3 min)
renders `output/37_hotspot_zoom/hotspots_miss_most_bright_cells.png` (1 mm on hotspot
IB1_55_162, DAPI + molecule density at level 2, clamp 21.1) and `IB1_57_161_zoom.png` (200 um
at full resolution, two panels: A the same overlay, B the step-35 panel style with every
polygon, molecule and ring). Both are embedded at the top of section 4 of the step-35 report.
**`HRS_hotspots_miss_most_bright_cells.png` in the project root is the superseded scratchpad
preview of the first figure** (gitignored); the step-37 copy is the one to use.
`IB1_57_161_box_molecules.csv` is the figure's one number: the box holds 305 molecules, 177
(58%) of them in four >= 20 polygons, 153 in one 90/38/25 three-polygon composite -- so this hotspot DOES contain the
brightest source in its field, and the blind spot is specifically the ISOLATED bright cell
(the 112-count cell at the window's left edge is outside every bin). Density bins find
clusters and grouped sources; they miss most of the isolated ones (fact in step 28 above).

## Conventions

- Region is the unit of every p-value. 154 neighbourhoods are not 400k independent
  observations — step 11 lost a Tfh result to exactly this.
- Circular results are reported twice, with and without the defining cells/genes. Only
  the second is evidence.
- **A percentile bootstrap p needs the `(1 + count)/(1 + n)` correction.** Scored as
  `2*min(mean(d<=0), mean(d>=0))` it returns **exactly 0** when every resample falls one side,
  and BH then reads 0 as infinite significance. In step 21 that alone manufactured **143
  genome-wide "hits"**, every one an artifact. Its floor is `2/(n+1)`, which over ~6,000 genes
  can never reach q < 0.05 — so a genome-wide scan needs a Wald p as well, with the percentile
  interval as the check on it.
- **Composition is a SHARE, never a count.** Total cells per bin falls with chemokine exposure
  inside IB1 FL tissue (rho −0.34 against −0.16 in background), so every cell type's count
  correlates negatively whether or not its share moves. Endothelium is the diagnostic: rho
  −0.378 on counts, **+0.236** on shares. A measure that flips sign with the denominator was
  measuring the denominator.
- **`tgt <- DT[cond]` takes a COPY; a later `DT[, new := ...]` does not reach it.** In step 23
  the target row was extracted before the rank column was added, so `tgt$rank_desc` was NULL,
  a derived value became zero-length, and `data.table()` silently produced a **0-row** table —
  surfacing as "nothing to report" rather than as an error. Add the column first, then
  extract, and `stopifnot()` the row.
- **A permutation null often has a closed form; check before simulating.** Step 22's first
  design was ~24 million permutations over million-element neighbour lists. With A the sparse
  adjacency and P the gene-by-cell positivity matrix, `OBS = P A P^T / sum(A)` and
  `EXP = Pbar (S A S^T) Pbar^T / sum(A)` are exact and take seconds.
- **Do not put an sbatch script or its inputs in `/tmp`.** Compute nodes cannot see the login
  node's node-local `/tmp`; the job fails instantly. `scripts/sbatch/` is the place. This is
  in `~/.claude/CLAUDE.md` and was still walked into during step 23's debugging.
- **An expensive job should be smoke-testable.** `21_script_interface_response.R` takes ~22 min
  and three early runs died on typos in its last few expressions, discarding all the compute.
  It now takes a `smoke` argument (tiny bootstrap, gene subsample, separate output directory)
  that exercises every code path in ~20 s. Run it before every real submit.
- Cluster `i` (over-fixation) is dropped per cell, never via the FOV-level
  `qc_spatial_i_fraction` flag, which is far too aggressive (it cuts IB1 to 196,316).
- `.Rds` objects hold **BPCells pointers** to on-disk matrix directories. They are not
  portable on their own; the `output/01_import_v1/` and `output/03_filter_and_process/bpcells/`
  paths must stay intact.
- Cell ids are `c_<slide>_<fov>_<cell>`; FOV is the second-to-last token.
  **The slide token is counter-intuitive: `c_2_` is IB1 and `c_1_` is RB1.** Never infer
  the sample from the prefix -- join on a sample column (`merge_sample`, `samp`,
  `Run_Tissue_name`) instead. Guessing it the obvious way round silently swaps the two
  samples in every row of the output and nothing errors.
- Sample column is `Run_Tissue_name` in step-01 objects, `merge_sample` in the step-04
  merged object, `sample_id` in the step-02 QC CSVs, `samp` in most derived tables.

## Environment

See `~/.claude/CLAUDE.md` for the cluster rules. In short: this is a SLURM login node,
anything over 1 CPU / 4 GB goes through `sbatch scripts/sbatch/run_rscript.sh`, and
reading any `.Rds` here exceeds that.

## Images: where they are and how to get into their coordinate space

**Assembled OME-TIFFs live on SCRATCH**, not with the project:

```
/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/{IB1,RB1}.ome.tiff   # 87 / 48 GiB
```

That path is `COSMX_TIFF_DIR` in `scripts/lib/cosmx_tiff.R:101`. It is **scratch** — subject
to purge policy, and re-assembling is expensive. Do not assume it is still there; check
before planning image work. `output/15_imagej/README.md` has the geometry (IB1 is
124,852 × 170,561 px, 5 channels) and the pyramid-level table.

`scripts/lib/cosmx_tiff.R` is the only supported way in. Read its "DO NOT SIMPLIFY" block
(lines 44–79) before changing anything in it. It depends on **`TiffPlotR`** (0.1.16) and
**`scPearsonPCA`**, which are local/GitHub packages, not CRAN — a fresh environment will
not have them.

### There are TWO coordinate transforms, and only one of them is analytic

| from | to | how |
|---|---|---|
| FOV-local px (`x_FOV_px`, `y_FOV_px`) | canvas px | `cosmx_local_to_px(sample, fov, x, y)` — analytic, from `input/fov_positions.csv` |
| slide mm (`x_slide_mm`, `y_slide_mm`) | canvas px | `cosmx_fit_mm_affine(sample, seu)` then `cosmx_mm_to_px(mm, aff)` — a **fitted** affine, not a constant |

There is no stored mm→px mapping, so it is regressed from cells that exist in both spaces.
`cosmx_fit_mm_affine()` already checks the three things that make a bad fit look like a
good one, and warns on each: both slopes **negative**, `|slope|` within 1% of
`1/COSMX_MM_PER_PX` (≈ 8,314 px/mm), residual SD under 1 px. If you fit it by hand — step 12
does, because it works from CSVs rather than an object — assert the same three
(`scripts/12_script_fl_vs_hodgkin_overlay.R:78-85`). A residual SD over 1 px means the
affine is wrong; do not trust the figure.

### THE AXIS RULE, in one line

**Canvas px reverses BOTH axes. Slide mm reverses NEITHER.** Nothing else.

Both fitted affine slopes are negative (`ax = ay ≈ -8314` px/mm), so slide mm runs opposite to
canvas px on *both* x and y. Canvas-px panels therefore need `scale_x_reverse()` *and*
`scale_y_reverse()` to render in image orientation; mm panels need neither, because reversing
twice returns you to the default. Getting this wrong mirrors a panel without erroring, and a
mirrored tissue map looks entirely plausible.

| space | x | y | canonical example |
|---|---|---|---|
| canvas px | `scale_x_reverse()` | `scale_y_reverse()` | `.coord_px()` + `.coord_win()`, `scripts/11_script_hodgkin_niche_views.R:96-110` |
| slide mm | plain | plain | the locator, `scripts/12_script_fl_vs_hodgkin_overlay.R:143-151` |

Step 16 got this wrong twice — the region maps and the Tfh view locator both applied
`scale_y_reverse()` in mm space, which mirrored them vertically against the image panels beside
them. Neither errored. If a tissue map looks upside down next to its image, this is why.

### Three further orientation traps, in the order people hit them

1. **Both affine slopes are negative** — the axes are mirrored, so the mm *minimum* maps to
   the px *maximum*. Use `cosmx_px_rect()`, which normalises the bounds, rather than
   passing mapped corners into `TiffRect()` (`cosmx_tiff.R:433`).
2. **`scale_*_reverse()`, never reversed `coord_*` limits.** ggplot2 sorts the limits given
   to `coord_fixed`/`coord_cartesian`, so `coord_fixed(xlim = c(max, min))` silently does
   nothing. Only the scale transform flips the axis. Verified on ggplot2 4.0.3; the one
   canonical definition is `.coord_px()` at `scripts/07_script_lightchain_views.R:576-591`.
3. **A mm panel beside a px panel reads mirrored** unless the rule above is applied to each —
   which is the whole reason the rule exists.

## Colour and channel conventions

These are arbitrary choices that several figures and the ImageJ ROI sets depend on, so
changing one silently breaks the correspondence between a figure and its caption. The
canonical *values* live in code — don't restate them, point at them — but the *decisions*
are these:

- **Channel order** is 1 PanCK (blue), 2 G (green), 3 Membrane (yellow), 4 CD45 (red),
  5 DNA (violet). Defined by `cosmx_channel_names()` / `cosmx_channel_colors()`. The
  OME-TIFF names them only B/G/Y/R/U; the marker names are the project's mapping.
- **Channel G's marker is UNKNOWN and is not CD68.** The vendor names it nowhere (`Mean.G` in
  the metadata). `scripts/15_script_channel_identity.R` ranks genes by per-cell correlation
  with each channel's intensity; raw correlation is a generic brightness axis (all five
  channels share top genes, CD45 ranks PTPRC only ~100-300th), so read the channel-SPECIFIC
  version (r minus mean r of the other four). On that, Membrane -> B2M/MHC I and DAPI ->
  nuclear genes behave, but G has no signal of its own: per-cell r 0.94 (IB1) / 0.83 (RB1)
  with PanCK, specific correlations <= +0.008 in RB1, and CD68 ranks 5080 / 3362 of ~6,170.
  Visually it is sparse ~1 um puncta, not macrophage cytoplasm. Do not caption it as a marker.
  `output/15_imagej/channel_identity/`.
- **Annotated ROI PNGs** (45 regions, full resolution, fixed crop ranges, hotspot/FL outlines,
  **G omitted by default**, `--with-G` restores it):
  `sbatch --mem=8G scripts/sbatch/run_python.sh scripts/15_roi_pngs.py IB1` ->
  `output/15_imagej/roi_png/` (~40 s, 4.2 GB peak). Reads the crop TIFFs, so it survives a
  scratch purge.
- **The two diseases have fixed hues: FL is orange `#eb6834`, Hodgkin is blue `#0b3d7a`.**
  Set in `scripts/12_script_fl_vs_hodgkin_overlay.R:51-52` and carried through step 16.
- **The two samples have fixed hues** in every comparison plot: `SAMP_COLS = c(IB1 = "#2a78d6",
  RB1 = "#eb6834")`.
- **ImageJ ROI colours: cyan = Hodgkin chemokine hotspot, orange = FL neighbourhood.** ROI
  names carry the metric (`c` = CCL17+CCL22 cp10k, `k` = kappa fraction) and flags
  (`_seam`). See `output/15_imagej/README.md`.
- **Diverging scales run blue -> red: red is the HIGH/positive end, blue the LOW/negative
  end.** In a differential comparison red marks the UP group and blue the DOWN group. This
  covers z-scores, log fold changes, and ratios about a midpoint — anything with a
  direction. It does **not** apply to one-hue sequential scales, which carry magnitude with
  no "down" direction and keep whatever hue suits the figure; the CCL17+CCL22 maps in
  `scripts/11_*` and `scripts/12_*` stay on Hodgkin blue, and the ImageJ ROI sets depend on
  that. Where a diverging scale would collide with a fixed categorical hue above, the
  categorical rule wins — flag the clash rather than silently flipping one of them.
- **Sequential scales are clamped to a percentile window, not the data range.** CCL17+CCL22
  per bin runs 0-272 cp10k post-QC with 80% of bins under 20, so a ramp stretched to the max
  renders four-fifths of the map in one pale band. Steps 11 and 12 clamp to the 5th-99th
  percentile with `oob = scales::squish` and say so in the legend.

- House palette and `theme_viz()` / `sf()` / `sc()` are copy-pasted into ~15 scripts rather
  than shared; the canonical copy is `scripts/07_script_lightchain_neighborhoods.R:126-147`.
  Contrast ranges differ between the overview and the crops **on purpose** (30th/99th
  percentile of non-zero pixels per resolution).
