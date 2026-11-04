"""Step 34, image stage: nuclear area measured directly from the DNA channel for the
step-31 cells (500 candidates + 500 area-matched controls), to validate the vendor's
NucArea and to compare candidate vs control nuclei from the picture itself.

Per cell: a 40 um full-resolution DNA window at the centroid, Gaussian-smoothed
(sigma 0.6 um); the nucleus = the connected region at or above HALF of the local peak
(max within 2 um of the centroid) that contains the peak. Half-max of the cell's own
peak separates touching nuclei better than a global threshold; a region over MERGE_UM2
is flagged as a probable merge with a neighbour. Measures: area (um^2), moment-based
aspect ratio, mean DNA inside, and the number of DNA local maxima inside the region
(sigma 1.5 um) as a crude lobation count.

numpy + pyvips only (no scipy / scikit-image in this environment).

  python scripts/34_nuclear_image_check.py output/31_radial_cd45/cells.csv output/34_nuclear_size
"""
import csv
import os
import sys

import numpy as np
import pyvips
import tifffile
import zarr

CELLS = sys.argv[1] if len(sys.argv) > 1 else "output/31_radial_cd45/cells.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "output/34_nuclear_size"
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff"
UM = 0.120280945
WIN_UM, SIG_UM, PEAK_UM, LOBE_SIG_UM = 40.0, 0.6, 2.0, 1.5
HALF = 0.5
MERGE_UM2 = 300.0
DNA = 4


def blur(a, sigma_px):
    v = pyvips.Image.new_from_memory(np.ascontiguousarray(a.astype(np.float32)).tobytes(),
                                     a.shape[1], a.shape[0], 1, "float")
    v = v.gaussblur(sigma_px)
    return np.ndarray(buffer=v.write_to_memory(), dtype=np.float32, shape=(v.height, v.width)).copy()


def flood(fg, seed):
    """Connected component of fg (bool) containing seed, via pyvips labelregions."""
    if not fg[seed]:
        return np.zeros_like(fg)
    v = pyvips.Image.new_from_memory(np.ascontiguousarray(fg.astype(np.uint8)).tobytes(),
                                     fg.shape[1], fg.shape[0], 1, "uchar")
    lab = v.labelregions()
    L = np.ndarray(buffer=lab.write_to_memory(), dtype=np.uint32, shape=(lab.height, lab.width))
    return (L == L[seed]) & fg


def local_maxima(a, mask):
    """Strict local maxima of a (8-neighbourhood) inside mask."""
    p = np.pad(a, 1, mode="edge")
    best = np.ones_like(a, dtype=bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy or dx:
                best &= a > p[1 + dy:1 + dy + a.shape[0], 1 + dx:1 + dx + a.shape[1]]
    return int((best & mask).sum())


def main():
    os.makedirs(OUT, exist_ok=True)
    with open(CELLS, newline="") as fh:
        cells = list(csv.DictReader(fh))
    R = int(round(WIN_UM / 2 / UM))
    n = 2 * R + 1
    rr = int(round(PEAK_UM / UM))
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    near = (xx ** 2 + yy ** 2) <= rr ** 2
    rows = []
    with tifffile.TiffFile(TIFF) as tf:
        store = tf.series[0].aszarr()
        z0 = zarr.open(store, mode="r")["0"]
        _, H, W = z0.shape
        for i, c in enumerate(cells):
            cx, cy = int(round(float(c["X_px"]))), int(round(float(c["Y_px"])))
            row = {"cell_id": c["cell_id"], "group": c["group"], "ok": False, "reason": ""}
            if cx - R < 0 or cy - R < 0 or cx + R + 1 > W or cy + R + 1 > H:
                row["reason"] = "off_canvas"; rows.append(row); continue
            a = np.asarray(z0[DNA, cy - R:cy + R + 1, cx - R:cx + R + 1]).astype(np.float32)
            if (a == 0).mean() > 0.5:
                row["reason"] = "blank"; rows.append(row); continue
            s = blur(a, SIG_UM / UM)
            pk = np.unravel_index(np.argmax(np.where(near, s, -np.inf)), s.shape)
            peak = float(s[pk])
            bg = float(np.percentile(s, 5))
            thr = bg + HALF * (peak - bg)
            if peak <= bg:
                row["reason"] = "no_peak"; rows.append(row); continue
            m = flood(s >= thr, pk)
            area_px = int(m.sum())
            ys, xs = np.nonzero(m)
            cov = np.cov(np.vstack([xs, ys])) if area_px > 2 else np.eye(2)
            ev = np.sort(np.linalg.eigvalsh(cov))
            ar = float(np.sqrt(ev[1] / max(ev[0], 1e-9)))
            s2 = blur(a, LOBE_SIG_UM / UM)
            nmax = local_maxima(s2, m)
            touches = bool(m[0, :].any() or m[-1, :].any() or m[:, 0].any() or m[:, -1].any())
            area_um2 = area_px * UM * UM
            row.update(ok=True, img_nuc_area=round(area_um2, 2), img_nuc_ar=round(ar, 3),
                       img_dna_mean=round(float(a[m].mean()), 1), img_dna_peak=round(peak, 1),
                       img_n_maxima=nmax, merge_flag=area_um2 > MERGE_UM2 or touches,
                       peak_offset_um=round(float(np.hypot(pk[1] - R, pk[0] - R) * UM), 2))
            rows.append(row)
            if (i + 1) % 200 == 0:
                print("  %d / %d" % (i + 1, len(cells)), flush=True)
        store.close()
    cols = ["cell_id", "group", "ok", "reason", "img_nuc_area", "img_nuc_ar", "img_dna_mean",
            "img_dna_peak", "img_n_maxima", "merge_flag", "peak_offset_um"]
    with open(os.path.join(OUT, "image_nuclei.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, restval="")
        wr.writeheader(); wr.writerows(rows)
    ok = [r for r in rows if r["ok"]]
    print("done: %d cells, %d measured, %d merge-flagged -> %s/image_nuclei.csv"
          % (len(rows), len(ok), sum(r["merge_flag"] for r in ok), OUT), flush=True)


if __name__ == "__main__":
    main()
