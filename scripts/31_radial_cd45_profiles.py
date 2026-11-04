"""Step 31, profile stage: mean CD45 and DNA in 5 um rings around each cell's centroid.

Reads <OUT>/cells.csv (scripts/31_script_radial_cd45.R select) and the IB1 OME-TIFF on
scratch at full resolution, and writes <OUT>/profiles.csv.gz: one row per cell, the mean
intensity in each 5 um band from 0 to 100 um for CD45 and DNA, plus a usable flag.

Pixels come straight from the image -- no segmentation boundary is used, which is the
point (step 30: segmented per-cell CD45 is mostly the neighbours'). A cell is dropped
(ok = False) if its 100 um disc leaves the canvas or more than MAX_BLANK of the disc is
blank in both channels (some IB1 FOVs carry transcripts but no morphology, CLAUDE.md).

  python scripts/31_radial_cd45_profiles.py output/31_radial_cd45
"""
import csv
import gzip
import os
import sys

import numpy as np
import tifffile
import zarr

OUT = sys.argv[1] if len(sys.argv) > 1 else "output/31_radial_cd45"
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff"
UM_PER_PX = 0.120280945
BAND_UM, R_MAX = 5.0, 100.0
NB = int(R_MAX / BAND_UM)
CH = {"cd45": 3, "dna": 4}          # positional, cosmx_channel_names(): PanCK G Membrane CD45 DNA
MAX_BLANK = 0.05


def main():
    with open(os.path.join(OUT, "cells.csv"), newline="") as fh:
        cells = list(csv.DictReader(fh))
    R = int(np.ceil(R_MAX / UM_PER_PX))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    r_um = np.sqrt(xx ** 2 + yy ** 2) * UM_PER_PX
    disc = r_um < R_MAX
    band = np.minimum((r_um / BAND_UM).astype(int), NB - 1)
    rows = []
    with tifffile.TiffFile(TIFF) as tf:
        store = tf.series[0].aszarr()
        z0 = zarr.open(store, mode="r")["0"]
        _, H, W = z0.shape
        for i, c in enumerate(cells):
            cx, cy = int(round(float(c["X_px"]))), int(round(float(c["Y_px"])))
            row = {"cell_id": c["cell_id"], "group": c["group"], "ok": True, "reason": "",
                   "blank_frac": ""}
            if cx - R < 0 or cy - R < 0 or cx + R + 1 > W or cy + R + 1 > H:
                row.update(ok=False, reason="off_canvas")
                rows.append(row)
                continue
            img = {k: np.asarray(z0[ch, cy - R:cy + R + 1, cx - R:cx + R + 1]).astype(np.float64)
                   for k, ch in CH.items()}
            blank = (img["cd45"] == 0) & (img["dna"] == 0)
            bf = float(blank[disc].mean())
            row["blank_frac"] = round(bf, 4)
            if bf > MAX_BLANK:
                row.update(ok=False, reason="blank")
                rows.append(row)
                continue
            use = disc & ~blank
            b = band[use]
            n = np.bincount(b, minlength=NB)
            row.update({"npx_b%02d" % j: int(n[j]) for j in range(NB)})
            for k in CH:
                s = np.bincount(b, weights=img[k][use], minlength=NB)
                row.update({"%s_b%02d" % (k, j): round(s[j] / max(n[j], 1), 3) for j in range(NB)})
            rows.append(row)
            if (i + 1) % 100 == 0:
                print("  %d / %d" % (i + 1, len(cells)), flush=True)
        store.close()
    cols = ["cell_id", "group", "ok", "reason", "blank_frac"] \
        + ["%s_b%02d" % (k, j) for k in CH for j in range(NB)] + ["npx_b%02d" % j for j in range(NB)]
    with gzip.open(os.path.join(OUT, "profiles.csv.gz"), "wt", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, restval="")
        wr.writeheader()
        wr.writerows(rows)
    print("done: %d cells, %d usable -> %s" % (len(rows), sum(r["ok"] for r in rows), OUT), flush=True)


if __name__ == "__main__":
    main()
