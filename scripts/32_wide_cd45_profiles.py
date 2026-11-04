"""Step 32, profile stage: mean CD45 and DNA in 25 um bands to 500 um, pyramid level 3.

Reads <OUT>/cells.csv (scripts/32_script_wide_cd45.R select) and writes
<OUT>/profiles.csv.gz: per cell, per band, the mean of the non-blank pixels, the blank
fraction (the analysis NAs a band over 5% blank), and ok/reason. A cell is dropped if its
400-500 um reference ring is more than half blank or the window leaves the canvas.

  python scripts/32_wide_cd45_profiles.py output/32_wide_cd45
"""
import csv
import gzip
import os
import sys

import numpy as np
import tifffile
import zarr

OUT = sys.argv[1] if len(sys.argv) > 1 else "output/32_wide_cd45"
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff"
LEVEL = 3
UM_PER_PX = 0.120280945 * 2 ** LEVEL
BAND_UM, R_MAX = 25.0, 500.0
NB = int(R_MAX / BAND_UM)
REF = [j for j in range(NB) if (j + 1) * BAND_UM > 400]
CH = {"cd45": 3, "dna": 4}
MAX_REF_BLANK = 0.5


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
        z = zarr.open(store, mode="r")[str(LEVEL)]
        _, H, W = z.shape
        f = 2 ** LEVEL
        for i, c in enumerate(cells):
            cx, cy = int(round(float(c["X_px"]) / f)), int(round(float(c["Y_px"]) / f))
            row = {"cell_id": c["cell_id"], "group": c["group"], "ok": True, "reason": ""}
            if cx - R < 0 or cy - R < 0 or cx + R + 1 > W or cy + R + 1 > H:
                row.update(ok=False, reason="off_canvas")
                rows.append(row)
                continue
            img = {k: np.asarray(z[ch, cy - R:cy + R + 1, cx - R:cx + R + 1]).astype(np.float64)
                   for k, ch in CH.items()}
            blank = (img["cd45"] == 0) & (img["dna"] == 0)
            b_all = band[disc]
            n_all = np.bincount(b_all, minlength=NB)
            n_blank = np.bincount(b_all, weights=blank[disc], minlength=NB)
            bf = n_blank / np.maximum(n_all, 1)
            row.update({"blank_b%02d" % j: round(float(bf[j]), 4) for j in range(NB)})
            if n_blank[REF].sum() / max(n_all[REF].sum(), 1) > MAX_REF_BLANK:
                row.update(ok=False, reason="ref_blank")
                rows.append(row)
                continue
            use = disc & ~blank
            b = band[use]
            n = np.bincount(b, minlength=NB)
            for k in CH:
                s = np.bincount(b, weights=img[k][use], minlength=NB)
                row.update({"%s_b%02d" % (k, j): round(s[j] / max(n[j], 1), 3) for j in range(NB)})
            rows.append(row)
            if (i + 1) % 100 == 0:
                print("  %d / %d" % (i + 1, len(cells)), flush=True)
        store.close()
    cols = ["cell_id", "group", "ok", "reason"] \
        + ["%s_b%02d" % (k, j) for k in CH for j in range(NB)] + ["blank_b%02d" % j for j in range(NB)]
    with gzip.open(os.path.join(OUT, "profiles.csv.gz"), "wt", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, restval="")
        wr.writeheader()
        wr.writerows(rows)
    print("done: %d cells, %d usable -> %s" % (len(rows), sum(r["ok"] for r in rows), OUT), flush=True)


if __name__ == "__main__":
    main()
