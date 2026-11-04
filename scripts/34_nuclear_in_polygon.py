"""Step 34, image stage (polygon-bounded): nuclear measurements INSIDE each cell's vendor
polygon, from the full-resolution DNA channel.

The centroid-anchored version (34_nuclear_image_check.py) failed because candidate
centroids often sit on DNA-dim pixels between nuclei, so a region grown from the centroid
is either a dim blob or a flood through the background. Bounding by the vendor polygon
asks the right question per cell: how much nucleus does THIS segmented object contain,
and what does it look like?

Per cell (40 um window around the centroid, polygon rasterised by the even-odd rule):
  poly_area        polygon area (um^2)
  nuc_thr          Otsu threshold on the smoothed DNA of the whole window (nuclei vs not)
  nuc_area         DNA >= threshold inside the polygon (um^2); nuc_frac = nuc_area / poly_area
  has_nucleus      nuc_area >= 8 um^2 (a lymphocyte nucleus is ~25-35)
  dna_centroid     smoothed DNA at the centroid, as a fraction of the window's Otsu threshold
  n_maxima         DNA local maxima (sigma 1.5 um) inside the nuclear mask -- lobation
  n_holes          dark components inside the nuclear mask (< 0.6 x the mask's median DNA,
                   >= 0.3 um^2) -- nucleoli / chromatin clearing
  nuc_ar           moment-based aspect ratio of the nuclear mask
  bright_outside   DNA >= threshold within 2 um OUTSIDE the polygon, as a fraction of the
                   polygon's perimeter band -- is a nucleus pressed against this object?
  PALE measures (added after the first QA sheet: half the candidate polygons enclose a large,
  pale, mottled object that the Otsu threshold -- set by dense lymphocyte nuclei -- misses;
  vesicular HRS-like chromatin would look exactly like that, and so would cytoplasm with
  a nucleus elsewhere; these separate the two):
  dna_rel          (median smoothed DNA inside the polygon - window background) /
                   (Otsu - background): 0 = background, 1 = dense-nucleus level
  pale_thr         background + 0.25 x (Otsu - background)
  pale_area        DNA >= pale_thr inside the polygon (um^2); pale_frac = / poly_area
  pale_ar          moment aspect ratio of the pale mask; pale_cv = sd / mean of DNA in it
  pale_n_holes     dark components inside the pale mask (< 0.6 x its median, >= 0.3 um^2)

  python scripts/34_nuclear_in_polygon.py output/34_nuclear_size
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np
import pyvips
import tifffile
import zarr

OUT = sys.argv[1] if len(sys.argv) > 1 else "output/34_nuclear_size"
CELLS = "output/31_radial_cd45/cells.csv"
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff"
UM = 0.120280945
WIN_UM, SIG_UM, LOBE_SIG_UM, RIM_UM = 40.0, 0.6, 1.5, 2.0
MIN_NUC_UM2, HOLE_FRAC, MIN_HOLE_UM2 = 8.0, 0.6, 0.3
DNA = 4


def blur(a, sigma_px):
    v = pyvips.Image.new_from_memory(np.ascontiguousarray(a.astype(np.float32)).tobytes(), a.shape[1], a.shape[0], 1, "float")
    v = v.gaussblur(sigma_px)
    return np.ndarray(buffer=v.write_to_memory(), dtype=np.float32, shape=(v.height, v.width)).copy()


def label(mask):
    v = pyvips.Image.new_from_memory(np.ascontiguousarray(mask.astype(np.uint8)).tobytes(), mask.shape[1], mask.shape[0], 1, "uchar")
    lab = v.labelregions()
    L = np.ndarray(buffer=lab.write_to_memory(), dtype=np.uint32, shape=(lab.height, lab.width)).copy()
    L[~mask] = 0
    return L


def otsu(v):
    h, e = np.histogram(v, bins=256)
    c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h); w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1); m1 = (np.cumsum((h * c)[::-1])[::-1] / np.maximum(w1, 1))
    m1 = np.append(m1[1:], 0)
    between = w0[:-1] * w1[:-1] * (m0[:-1] - m1[:-1]) ** 2
    return float(c[np.argmax(between)])


def rasterise(px, py, shape):
    """Even-odd point-in-polygon on a pixel grid; px, py in window coordinates."""
    h, w = shape
    yy, xx = np.mgrid[0:h, 0:w]
    inside = np.zeros(shape, bool)
    n = len(px)
    for i in range(n):
        x1, y1, x2, y2 = px[i], py[i], px[(i + 1) % n], py[(i + 1) % n]
        if y1 == y2:
            continue
        cond = (yy >= min(y1, y2)) & (yy < max(y1, y2))
        xint = x1 + (yy - y1) * (x2 - x1) / (y2 - y1)
        inside ^= cond & (xx < xint)
    return inside


def dilate(m, r):
    v = pyvips.Image.new_from_memory(np.ascontiguousarray(m.astype(np.uint8) * 255).tobytes(), m.shape[1], m.shape[0], 1, "uchar")
    k = 2 * r + 1
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    mask = pyvips.Image.new_from_array(((xx ** 2 + yy ** 2) <= r * r).astype(np.uint8).tolist())
    mask = (mask * 255).cast("uchar")
    v = v.morph(mask, "dilate")
    return np.ndarray(buffer=v.write_to_memory(), dtype=np.uint8, shape=(v.height, v.width)) > 0


def local_maxima(a, mask):
    p = np.pad(a, 1, mode="edge")
    best = np.ones_like(a, dtype=bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy or dx:
                best &= a > p[1 + dy:1 + dy + a.shape[0], 1 + dx:1 + dx + a.shape[1]]
    return int((best & mask).sum())


def main():
    with open(CELLS, newline="") as fh:
        cells = list(csv.DictReader(fh))
    poly = defaultdict(list)
    with open(os.path.join(OUT, "polygons_canvas.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            poly[r["cell_id"]].append((float(r["X_px"]), float(r["Y_px"])))
    R = int(round(WIN_UM / 2 / UM)); n = 2 * R + 1
    rim = int(round(RIM_UM / UM))
    rows = []
    with tifffile.TiffFile(TIFF) as tf:
        store = tf.series[0].aszarr()
        z0 = zarr.open(store, mode="r")["0"]
        _, H, W = z0.shape
        for i, c in enumerate(cells):
            row = {"cell_id": c["cell_id"], "group": c["group"], "ok": False, "reason": ""}
            if c["cell_id"] not in poly:
                row["reason"] = "no_polygon"; rows.append(row); continue
            P = np.array(poly[c["cell_id"]])
            cx, cy = int(round(float(c["X_px"]))), int(round(float(c["Y_px"])))
            x0, y0 = cx - R, cy - R
            if x0 < 0 or y0 < 0 or x0 + n > W or y0 + n > H:
                row["reason"] = "off_canvas"; rows.append(row); continue
            px, py = P[:, 0] - x0, P[:, 1] - y0
            if px.min() < 0 or py.min() < 0 or px.max() >= n or py.max() >= n:
                row["reason"] = "polygon_outside_window"; rows.append(row); continue
            a = np.asarray(z0[DNA, y0:y0 + n, x0:x0 + n]).astype(np.float32)
            if (a == 0).mean() > 0.5:
                row["reason"] = "blank"; rows.append(row); continue
            s = blur(a, SIG_UM / UM)
            thr = otsu(s[s > 0])
            inside = rasterise(px, py, s.shape)
            poly_area = float(inside.sum()) * UM * UM
            if poly_area < 1:
                row["reason"] = "degenerate_polygon"; rows.append(row); continue
            nuc = inside & (s >= thr)
            nuc_area = float(nuc.sum()) * UM * UM
            ring = dilate(inside, rim) & ~inside
            bright_outside = float((ring & (s >= thr)).sum()) / max(float(ring.sum()), 1.0)
            bg = float(np.percentile(s, 5))
            pale_thr = bg + 0.25 * (thr - bg)
            pale = inside & (s >= pale_thr)
            pale_area = float(pale.sum()) * UM * UM
            dna_rel = (float(np.median(s[inside])) - bg) / max(thr - bg, 1e-6)
            row.update(dna_rel=round(dna_rel, 3), pale_thr=round(pale_thr, 1), pale_area=round(pale_area, 2),
                       pale_frac=round(pale_area / poly_area, 4))
            if pale_area >= MIN_NUC_UM2:
                ys, xs = np.nonzero(pale); cov = np.cov(np.vstack([xs, ys])); ev = np.sort(np.linalg.eigvalsh(cov))
                pm = float(np.median(s[pale]))
                ph = label(pale & (s < HOLE_FRAC * pm)); _, phn = np.unique(ph[ph > 0], return_counts=True)
                row.update(pale_ar=round(float(np.sqrt(ev[1] / max(ev[0], 1e-9))), 3),
                           pale_cv=round(float(s[pale].std() / max(s[pale].mean(), 1e-6)), 4),
                           pale_n_holes=int((phn * UM * UM >= MIN_HOLE_UM2).sum()))
            row.update(ok=True, poly_area=round(poly_area, 2), nuc_thr=round(thr, 1),
                       nuc_area=round(nuc_area, 2), nuc_frac=round(nuc_area / poly_area, 4),
                       has_nucleus=nuc_area >= MIN_NUC_UM2,
                       dna_centroid=round(float(s[R, R]) / thr, 3), bright_outside=round(bright_outside, 3),
                       n_vertex=len(P))
            if nuc_area >= MIN_NUC_UM2:
                s2 = blur(a, LOBE_SIG_UM / UM)
                ys, xs = np.nonzero(nuc)
                cov = np.cov(np.vstack([xs, ys])); ev = np.sort(np.linalg.eigvalsh(cov))
                med = float(np.median(s[nuc]))
                holes = label(nuc & (s < HOLE_FRAC * med))
                hid, hn = np.unique(holes[holes > 0], return_counts=True)
                row.update(n_maxima=local_maxima(s2, nuc), nuc_ar=round(float(np.sqrt(ev[1] / max(ev[0], 1e-9))), 3),
                           n_holes=int((hn * UM * UM >= MIN_HOLE_UM2).sum()),
                           nuc_dna_med=round(med, 1), n_nuc_components=int(len(np.unique(label(nuc))) - 1))
            rows.append(row)
            if (i + 1) % 200 == 0:
                print("  %d / %d" % (i + 1, len(cells)), flush=True)
        store.close()
    cols = ["cell_id", "group", "ok", "reason", "poly_area", "nuc_thr", "nuc_area", "nuc_frac", "has_nucleus",
            "dna_centroid", "bright_outside", "n_vertex", "n_maxima", "nuc_ar", "n_holes", "nuc_dna_med", "n_nuc_components",
            "dna_rel", "pale_thr", "pale_area", "pale_frac", "pale_ar", "pale_cv", "pale_n_holes"]
    with open(os.path.join(OUT, "nuclei_in_polygon.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, restval="")
        wr.writeheader(); wr.writerows(rows)
    ok = [r for r in rows if r["ok"]]
    print("done: %d cells, %d measured, %d with a nucleus (>= %g um^2) -> %s/nuclei_in_polygon.csv"
          % (len(rows), len(ok), sum(r["has_nucleus"] for r in ok), MIN_NUC_UM2, OUT), flush=True)


if __name__ == "__main__":
    main()
