"""Step 35, adjacency stage: one FOV per SLURM array task.

Reads <OUT>/fov_list.csv (row = $SLURM_ARRAY_TASK_ID) and <OUT>/polygons/<samp>_<fov>.csv
(vendor polygon vertices, FOV-local px), rasterises every polygon into a 4256 x 4256
label image (even-odd rule), and writes
  <OUT>/adj/<samp>_<fov>_edges.csv   touching pairs (a < b) with contact_px = number of
                                     label changes across 1-3 px shifts (corners excluded
                                     downstream by a minimum contact)
  <OUT>/adj/<samp>_<fov>_cells.csv   per polygon: area px, dense / pale DNA px (full-
                                     resolution DNA channel, FOV Otsu and background),
                                     mean DNA on a 0 (background) - 1 (Otsu) scale,
                                     distance to the FOV edge (um), fov_blank flag

Local -> canvas is the analytic transform of cosmx_tiff.R::cosmx_local_to_px:
  X = x0 + (4255 - x_local), Y = y0 + y_local  (x is mirrored).

  python scripts/35_fov_adjacency.py <OUT>         # uses $SLURM_ARRAY_TASK_ID
  python scripts/35_fov_adjacency.py <OUT> <task>
"""
import csv
import os
import sys
from collections import defaultdict

import numpy as np
import pyvips
import tifffile
import zarr

OUT = sys.argv[1]
TASK = int(sys.argv[2]) if len(sys.argv) > 2 else int(os.environ["SLURM_ARRAY_TASK_ID"])
TIFF_DIR = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs"
FOV_PX = 4256
UM = 0.120280945
SIG_UM = 0.6
DNA = 4
SHIFTS = [(0, 1), (1, 0), (0, 2), (2, 0), (0, 3), (3, 0)]


def rasterise_into(lab, px, py, value):
    x0, x1 = max(int(np.floor(px.min())), 0), min(int(np.ceil(px.max())) + 1, FOV_PX)
    y0, y1 = max(int(np.floor(py.min())), 0), min(int(np.ceil(py.max())) + 1, FOV_PX)
    if x1 <= x0 or y1 <= y0:
        return 0
    yy, xx = np.mgrid[y0:y1, x0:x1]
    inside = np.zeros(yy.shape, bool)
    n = len(px)
    for i in range(n):
        xa, ya, xb, yb = px[i], py[i], px[(i + 1) % n], py[(i + 1) % n]
        if ya == yb:
            continue
        cond = (yy >= min(ya, yb)) & (yy < max(ya, yb))
        xint = xa + (yy - ya) * (xb - xa) / (yb - ya)
        inside ^= cond & (xx < xint)
    sub = lab[y0:y1, x0:x1]
    free = inside & (sub == 0)            # first polygon wins where vendor polygons overlap
    sub[free] = value
    return int(free.sum())


def otsu(v):
    h, e = np.histogram(v, bins=256)
    c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h); w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1)
    m1 = np.append((np.cumsum((h * c)[::-1])[::-1] / np.maximum(w1, 1))[1:], 0)
    between = w0[:-1] * w1[:-1] * (m0[:-1] - m1[:-1]) ** 2
    return float(c[np.argmax(between)])


def main():
    with open(os.path.join(OUT, "fov_list.csv"), newline="") as fh:
        row = [r for r in csv.DictReader(fh) if int(r["task"]) == TASK][0]
    samp, fov, X0, Y0 = row["samp"], int(row["fov"]), int(float(row["x0"])), int(float(row["y0"]))
    poly = defaultdict(list)
    with open(os.path.join(OUT, "polygons", "%s_%d.csv" % (samp, fov)), newline="") as fh:
        for r in csv.reader(fh):
            poly[int(r[1])].append((FOV_PX - 1 - float(r[2]), float(r[3])))      # (wx, wy)
    lab = np.zeros((FOV_PX, FOV_PX), np.uint32)
    area = {}
    for cid, pts in poly.items():
        P = np.array(pts)
        area[cid] = rasterise_into(lab, P[:, 0], P[:, 1], cid)
    # touching pairs
    pairs = defaultdict(int)
    for dy, dx in SHIFTS:
        a = lab[dy:, dx:]; b = lab[:FOV_PX - dy, :FOV_PX - dx]
        k = (a > 0) & (b > 0) & (a != b)
        if k.any():
            lo = np.minimum(a[k], b[k]).astype(np.int64); hi = np.maximum(a[k], b[k]).astype(np.int64)
            key, cnt = np.unique(lo * (1 << 32) + hi, return_counts=True)
            for kk, cc in zip(key, cnt):
                pairs[(int(kk >> 32), int(kk & 0xFFFFFFFF))] += int(cc)
    # DNA
    with tifffile.TiffFile(os.path.join(TIFF_DIR, "%s.ome.tiff" % samp)) as tf:
        store = tf.series[0].aszarr(); z0 = zarr.open(store, mode="r")["0"]
        _, H, W = z0.shape
        ys, ye, xs, xe = max(Y0, 0), min(Y0 + FOV_PX, H), max(X0, 0), min(X0 + FOV_PX, W)
        dna = np.zeros((FOV_PX, FOV_PX), np.float32)
        dna[ys - Y0:ye - Y0, xs - X0:xe - X0] = np.asarray(z0[DNA, ys:ye, xs:xe]).astype(np.float32)
        store.close()
    blank = bool((dna == 0).mean() > 0.5)
    if not blank:
        v = pyvips.Image.new_from_memory(np.ascontiguousarray(dna).tobytes(), FOV_PX, FOV_PX, 1, "float").gaussblur(SIG_UM / UM)
        s = np.ndarray(buffer=v.write_to_memory(), dtype=np.float32, shape=(FOV_PX, FOV_PX))
        thr = otsu(s[s > 0]); bg = float(np.percentile(s[s > 0], 5)); pale_thr = bg + 0.25 * (thr - bg)
        ids = lab.ravel(); n_lab = int(lab.max()) + 1
        dense_px = np.bincount(ids, weights=(s >= thr).ravel(), minlength=n_lab)
        pale_px = np.bincount(ids, weights=(s >= pale_thr).ravel(), minlength=n_lab)
        sum_s = np.bincount(ids, weights=s.ravel(), minlength=n_lab)
        npx = np.bincount(ids, minlength=n_lab)
    os.makedirs(os.path.join(OUT, "adj"), exist_ok=True)
    with open(os.path.join(OUT, "adj", "%s_%d_cells.csv" % (samp, fov)), "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["fov", "cell_ID", "poly_px", "dense_px", "pale_px", "dna_rel", "edge_um", "fov_blank"])
        for cid, pts in poly.items():
            P = np.array(pts)
            edge = min(P[:, 0].min(), P[:, 1].min(), FOV_PX - 1 - P[:, 0].max(), FOV_PX - 1 - P[:, 1].max()) * UM
            if blank or area[cid] == 0:
                wr.writerow([fov, cid, area[cid], "", "", "", round(float(edge), 2), blank]); continue
            rel = (sum_s[cid] / max(npx[cid], 1) - bg) / max(thr - bg, 1e-6)
            wr.writerow([fov, cid, area[cid], int(dense_px[cid]), int(pale_px[cid]), round(float(rel), 4), round(float(edge), 2), blank])
    with open(os.path.join(OUT, "adj", "%s_%d_edges.csv" % (samp, fov)), "w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["fov", "a", "b", "contact_px"])
        for (a, b), c in sorted(pairs.items()):
            wr.writerow([fov, a, b, c])
    print("%s FOV %d: %d polygons, %d touching pairs, blank=%s" % (samp, fov, len(poly), len(pairs), blank), flush=True)


if __name__ == "__main__":
    main()
