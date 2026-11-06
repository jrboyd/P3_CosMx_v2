"""Pathology review bundle, stage 2 of 3: cut the image windows.

    sbatch --mem=24G --time=1:00:00 scripts/sbatch/run_python.sh \
        pathology_review_app/build/02_render_tiles.py test

Reads pathology_review_app/output/<size>/selection.csv (stage 1) and the assembled
OME-TIFFs on scratch; writes, for each accepted cell, one 8-bit greyscale JPEG per
channel to output/<size>/tiles/<cell_id>/<channel>.jpg, plus rendered.csv and ranges.json.

WHAT IS DELIBERATELY NOT ON THESE TILES, unlike the step-29 ROI PNGs:
  - no caption (it printed counts, score and group -- the answer key);
  - no CCL17/CCL22 density rings (controls would be the only cells without rings, so the
    rings would unblind candidate vs control at a glance);
  - no hotspot / FL outlines (RB1 has none, so they would identify the slide);
  - no crosshair or scale bar -- the app draws both as overlays the reader can hide.
The cell's centroid is the window centre (within half a pixel; exact offset in rendered.csv).

CONTRAST. Per slide, the 30th / 99th percentile of NON-ZERO pixels per channel (the rule of
15_imagej_tiles.py), measured on 16 random matched-tier cells' windows chosen in stage 1 --
so it does not depend on which cells are reviewed or on the bundle size. Per-slide, not one
shared range: IB1 has a dim-membrane staining defect, and a shared range would make the
slide recognisable and its cells harder to read. Channel G is not exported (marker unknown).

ACCEPTANCE, per stratum in priority order until `need` are kept: window inside the canvas,
<= 5% blank pixels in every channel (some IB1 FOVs carry transcripts but no morphology),
and not overlapping a window already kept on that slide (a reader must never see the same
tissue twice except in a planned repeat).
"""
import csv
import json
import os
import sys

import numpy as np
import pyvips
import tifffile
import zarr

SIZE = sys.argv[1] if len(sys.argv) > 1 else "test"
OUT = os.path.join("pathology_review_app/output", SIZE)
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/%s.ome.tiff"
UM_PER_PX = 0.120280945
WIN_UM = 150.0
N = int(round(WIN_UM / UM_PER_PX))                     # 1247 px
MAX_BLANK = 0.05
LO_PCT, HI_PCT = 30.0, 99.0
JPEG_Q = 88
# positional channel order of the OME-TIFF (CLAUDE.md "Colour and channel conventions")
CHANNELS = ["PanCK", "G", "Membrane", "CD45", "DNA"]
EXPORT = ["DNA", "CD45", "Membrane", "PanCK"]


def log(m):
    print(m, flush=True)


def window(z0, cx, cy):
    x0, y0 = int(round(cx)) - N // 2, int(round(cy)) - N // 2
    _, H, W = z0.shape
    if x0 < 0 or y0 < 0 or x0 + N > W or y0 + N > H:
        return None, x0, y0
    return np.asarray(z0[:, y0:y0 + N, x0:x0 + N]), x0, y0


def blank_frac(arr):
    return float((arr == 0).mean(axis=(1, 2)).max())


def slide_ranges(z0, rows):
    samples = [[] for _ in CHANNELS]
    used = 0
    for r in rows:
        arr, _, _ = window(z0, float(r["X_px"]), float(r["Y_px"]))
        if arr is None or blank_frac(arr) > 0.5:
            continue
        used += 1
        for c in range(len(CHANNELS)):
            v = arr[c, ::4, ::4].ravel()
            samples[c].append(v[v > 0])
    if used < 8:
        raise SystemExit("only %d usable range windows" % used)
    rng = {}
    for c, nm in enumerate(CHANNELS):
        nz = np.concatenate(samples[c])
        rng[nm] = (float(np.percentile(nz, LO_PCT)), float(np.percentile(nz, HI_PCT)))
    return rng, used


def save_gray(a, lo, hi, path):
    v = np.clip((a.astype(np.float32) - lo) / (hi - lo), 0.0, 1.0) * 255.0
    a8 = np.ascontiguousarray(v.round().astype(np.uint8))
    pyvips.Image.new_from_memory(a8.tobytes(), N, N, 1, "uchar").jpegsave(
        path, Q=JPEG_Q, strip=True)


def main():
    with open(os.path.join(OUT, "selection.csv"), newline="") as fh:
        sel = list(csv.DictReader(fh))
    for r in sel:
        r["priority"] = int(r["priority"])
        r["need"] = int(r["need"])
    os.makedirs(os.path.join(OUT, "tiles"), exist_ok=True)
    rendered, ranges = [], {}
    for samp in ("IB1", "RB1"):
        rows = [r for r in sel if r["samp"] == samp]
        with tifffile.TiffFile(TIFF % samp) as tf:
            store = tf.series[0].aszarr()
            z0 = zarr.open(store, mode="r")["0"]
            rng, used = slide_ranges(z0, [r for r in rows if r["kind"] == "range"])
            ranges[samp] = rng
            log("%s ranges from %d windows: %s" % (samp, used,
                ", ".join("%s %.0f-%.0f" % (k, *rng[k]) for k in EXPORT)))
            kept_xy = []
            # candidates before controls, so a control never displaces a candidate
            strata = sorted({(r["kind"], r["stratum"]) for r in rows if r["kind"] != "range"})
            for kind, stratum in strata:
                srows = sorted([r for r in rows if r["stratum"] == stratum],
                               key=lambda r: r["priority"])
                need = srows[0]["need"]
                kept, rej = 0, {"canvas": 0, "blank": 0, "overlap": 0}
                for r in srows:
                    if kept == need:
                        break
                    cx, cy = float(r["X_px"]), float(r["Y_px"])
                    if any(abs(cx - x) < N and abs(cy - y) < N for x, y in kept_xy):
                        rej["overlap"] += 1
                        continue
                    arr, x0, y0 = window(z0, cx, cy)
                    if arr is None:
                        rej["canvas"] += 1
                        continue
                    b = blank_frac(arr)
                    if b > MAX_BLANK:
                        rej["blank"] += 1
                        continue
                    d = os.path.join(OUT, "tiles", r["cell_id"])
                    os.makedirs(d, exist_ok=True)
                    for nm in EXPORT:
                        c = CHANNELS.index(nm)
                        save_gray(arr[c], *rng[nm], os.path.join(d, nm + ".jpg"))
                    kept += 1
                    kept_xy.append((cx, cy))
                    rendered.append(dict(r, x0=x0, y0=y0, cx_win=round(cx - x0, 2),
                                         cy_win=round(cy - y0, 2), blank_frac=round(b, 4)))
                log("  %-22s kept %d of need %d (rejected: canvas %d, blank %d, overlap %d)"
                    % (stratum, kept, need, rej["canvas"], rej["blank"], rej["overlap"]))
                if kept < need:
                    log("  WARNING %s short by %d -- add backups in stage 1" % (stratum, need - kept))
            store.close()

    with open(os.path.join(OUT, "rendered.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rendered[0]))
        wr.writeheader()
        wr.writerows(rendered)
    with open(os.path.join(OUT, "ranges.json"), "w") as fh:
        json.dump({"um_per_px": UM_PER_PX, "window_px": N, "lo_pct": LO_PCT, "hi_pct": HI_PCT,
                   "ranges": ranges}, fh, indent=1)
    log("done: %d cells rendered into %s/tiles" % (len(rendered), OUT))


if __name__ == "__main__":
    main()
