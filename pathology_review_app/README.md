# Pathology review app

Blinded morphology review of chemokine-bright (HRS-candidate) cells and matched controls,
by a pathologist, on their own computer. Answers the G2 gap in `reports/MANUSCRIPT_PLAN.md`
(no morphological confirmation) as far as four IF stains can.

```
app/reviewer.html        the app: ONE self-contained file, no network, no install
app/READ_ME_FIRST.txt    the reader's quick start (copied into every bundle's send/ folder)
build/01_select_cells.R  which cells (candidates, matched controls, range cells)
build/02_render_tiles.py one greyscale JPEG per channel per cell, from the scratch OME-TIFFs
build/03_make_bundle.py  blind ids, stratified random order, repeats, STORED zip, key
build/review_form.json   the checklist / 1-5 score -- edit here, not in the app
build/test_app_headless.py  end-to-end test in headless Firefox over WebDriver BiDi
build/submit_build.sh    runs all of the above as a SLURM chain
ingest/ingest_reviews.R  returned CSVs -> unblinded answers + summaries
keys/<bundle_id>_key.csv item_id -> cell_id, sample, stratum, counts, score. NEVER send.
returns/                 put the reader's exported CSVs here
output/                  (gitignored) tiles, bundles, test screenshots, review summaries
```

## Build and send

```bash
bash pathology_review_app/build/submit_build.sh test    # 50 items; or `full` (~220)
cat pathology_review_app/output/test/apptest/results.txt  # must say ALL PASS
```

Send the three files in `output/<size>/send/` (`reviewer.html`, `P3HRS-*.zip`,
`READ_ME_FIRST.txt`). The zip is too big for most email at full size; use a file share.
When CSVs come back: drop them in `returns/` and run
`Rscript --vanilla pathology_review_app/ingest/ingest_reviews.R`.

## Design decisions (and why)

- **Blind.** The step-29 PNGs print the group, counts and score on every tile; these tiles
  carry nothing. Ids are random 4-character codes; repeats get fresh ids and their own file
  copies. The key stays in `keys/`.
- **No chemokine rings, no hotspot/FL outlines.** Controls would be the only ring-free
  tiles, and RB1 has no outlines, so either would unblind at a glance. The crosshair and
  scale bar are drawn by the app as overlays the reader can hide.
- **Controls** are chemo < 10, defect-free, matched-tier cells of the same sample, each
  matched to a candidate on area (±10%) and depth0 (±15%). Without them a "yes" rate has no
  reference.
- **RB1** contributes candidates (its 9 cells >= 50, plus 10-49 by top score and at random)
  and its own matched controls: bright cells from a section with no Hodgkin component.
- **Per-slide display range** (30th/99th pct of non-zero pixels over 16 random cells per
  slide, fixed per slide, independent of the bundle). A shared range would make IB1 (dim
  membrane defect) recognisable and harder to read. The reader can still adjust brightness
  per channel in the app; that is display-only and not recorded.
- **Order**: stratified interleave, seeded, so every stretch mixes all strata and reading
  order is independent of brightness and sample. Not ranked: a ranked order confounds
  fatigue and anchoring with the score being tested.
- **Repeats** (~10%) re-shown at least max(10, n/3) items later, for intra-reader agreement.
- **Tri-state answers** (Yes / No / Unsure) with an image-adequacy gate first: blur and
  fixation artefact are common here (CLAUDE.md step 34) and a forced binary turns them into
  guesses.

## How saving works (what to tell the reader)

Every Save writes to the browser's localStorage, keyed by bundle id, immediately; it
survives closing the browser. It is lost only by clearing site data, a private window, or
changing browser/computer. Defences: a leave-page warning while anything is unexported;
in Chrome/Edge an optional "Autosave to file" that rewrites a CSV after every save; and
"Import a previous export" restores a wiped browser. Records are append-only: an edit is
revision 2, and the export keeps every revision with an `is_latest` flag. Each record has
`session_id` (one per app opening), `saved_at` (local time with offset) and `saved_at_utc`
for batch identification.

Two browser caveats, not testable on the cluster: Chrome treats every `file://` page as one
origin, so moving `reviewer.html` keeps the answers; Firefox may key storage to the file's
location, so the reader should not move the file between sessions. Safari unzips downloads
automatically, which is why the app also opens an unzipped folder.
