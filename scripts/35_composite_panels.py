"""Step 35 review panels: the largest multi-member seed composites on morphology.

Reads <OUT>/polygon_assignments.csv.gz, seed_composites.csv.gz, fov_list.csv and the
per-FOV polygon files, and renders for each selected composite a full-resolution window
(DNA / CD45 / Membrane composite as in the step-15 ROI PNGs, G omitted, Membrane x0.5)
with every polygon in the window outlined: members of the composite in chemokine-graded
warm colours (the nucleus-bearing member in cyan), other chemokine-positive polygons in
dim orange, the rest grey. A caption carries the composite's members, area, nuclei,
chemokine total and seam flag. One contact sheet per group.

CCL17+CCL22 molecules are drawn as dots when the sample's molecule table exists
(output/29_bright_cells/ib1_chemokine_tx.csv.gz, canvas px, IB1 only): magenta = assigned
to a composite member, violet = assigned to another cell, white = unassigned; plus the
step-29 density contour rings (sigma 8 um; 6 / 12 / 24 molecules per 100 um^2). The
question they answer: do the fragments sit inside the same molecule cloud as the anchor
(one cell carved up) or carry clouds of their own with a gap between (separate cells)?

  python scripts/35_composite_panels.py <OUT> <samp> <thr> <partition: anch|comp> [n]
"""
import csv
import gzip
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pyvips
import tifffile
import zarr

spec = importlib.util.spec_from_file_location("rp", "scripts/15_roi_pngs.py")
rp = importlib.util.module_from_spec(spec); spec.loader.exec_module(rp)
spec2 = importlib.util.spec_from_file_location("np34", "scripts/34_nuclear_in_polygon.py")
np34 = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(np34)
spec3 = importlib.util.spec_from_file_location("r29", "scripts/29_bright_cell_rois.py")
r29 = importlib.util.module_from_spec(spec3); spec3.loader.exec_module(r29)

OUT, SAMP, THR, PART = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
N_TOP = int(sys.argv[5]) if len(sys.argv) > 5 else 24
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/%s.ome.tiff" % SAMP
UM = rp.UM_PER_PX
FOV_PX = 4256
PAD_UM = 25.0
OMIT = ("G",)
WHITE, GREY, AMBER, CYAN = (255, 255, 255), (200, 200, 200), (255, 200, 60), (0, 229, 255)
MAGENTA, VIOLET = (255, 0, 255), (150, 110, 255)
TX_TABLE = {"IB1": "output/29_bright_cells/ib1_chemokine_tx.csv.gz"}
DOT = 3            # half-width of a molecule dot in px (7 px square = 0.84 um)
RING_T = 3


def dot(img, x, y, rgb):
    h, w = img.shape[:2]
    img[max(y - DOT, 0):min(y + DOT + 1, h), max(x - DOT, 0):min(x + DOT + 1, w)] = rgb


def rcsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    fl = {(r["samp"], int(r["fov"])): (int(float(r["x0"])), int(float(r["y0"]))) for r in rcsv(os.path.join(OUT, "fov_list.csv"))}
    asg = [r for r in rcsv(os.path.join(OUT, "polygon_assignments.csv.gz")) if r["samp"] == SAMP and int(r["thr"]) == THR]
    key = "anch" if PART == "anch" else "comp"
    comps = defaultdict(list)
    for r in asg:
        comps[int(r[key])].append(r)
    seedc = [(cid, m) for cid, m in comps.items() if any(r["seed"] == "TRUE" for r in m) and len(m) >= 2]
    seedc.sort(key=lambda t: -sum(float(r["poly_um2"]) for r in t[1]))
    seedc = seedc[:N_TOP]
    print("%s thr %d %s: %d multi-member seed composites, rendering %d" % (SAMP, THR, PART, len([1 for c, m in comps.items() if any(r['seed'] == 'TRUE' for r in m) and len(m) >= 2]), len(seedc)), flush=True)
    ranges = rp.crop_ranges("output/15_imagej/crops/IB1_w00_IB1_nb7269.tif")
    tx = None
    if SAMP in TX_TABLE and os.path.exists(TX_TABLE[SAMP]):
        t = rcsv(TX_TABLE[SAMP])
        tx = {"x": np.array([int(float(r["X_px"])) for r in t]), "y": np.array([int(float(r["Y_px"])) for r in t]),
              "fov": np.array([int(float(r["fov"])) for r in t]), "cell": np.array([int(float(r["cell_ID"])) for r in t])}
        print("molecules:", len(tx["x"]), flush=True)
    tiles = []
    with tifffile.TiffFile(TIFF) as tf:
        store = tf.series[0].aszarr(); z0 = zarr.open(store, mode="r")["0"]; _, H, W = z0.shape
        pcache = {}
        for cid, mem in seedc:
            fov = int(mem[0]["fov"]); x0f, y0f = fl[(SAMP, fov)]
            if fov not in pcache:
                poly = defaultdict(list)
                with open(os.path.join(OUT, "polygons", "%s_%d.csv" % (SAMP, fov)), newline="") as fh:
                    for r in csv.reader(fh):
                        poly[int(r[1])].append((x0f + FOV_PX - 1 - float(r[2]), y0f + float(r[3])))
                pcache[fov] = poly
            poly = pcache[fov]
            ids = [int(r["cell_ID"]) for r in mem]
            pts = np.vstack([np.array(poly[i]) for i in ids if i in poly])
            pad = int(round(PAD_UM / UM))
            x0, y0 = int(pts[:, 0].min()) - pad, int(pts[:, 1].min()) - pad
            x1, y1 = int(pts[:, 0].max()) + pad, int(pts[:, 1].max()) + pad
            side = max(x1 - x0, y1 - y0); cx, cy = (x0 + x1) // 2, (y0 + y1) // 2
            x0, y0 = max(cx - side // 2, 0), max(cy - side // 2, 0); n = side
            if x0 + n > W or y0 + n > H:
                continue
            arr = np.asarray(z0[:, y0:y0 + n, x0:x0 + n])
            img = rp.composite(arr, ranges, OMIT); np.clip(img, 0, 255, out=img)
            chemo = {int(r["cell_ID"]): float(r["chemo"]) for r in asg if int(r["fov"]) == fov}
            memset = {int(r["cell_ID"]): r for r in mem}
            def outline(pid, rgb, t):
                P = np.array(poly[pid]); px, py = P[:, 0] - x0, P[:, 1] - y0
                if px.max() < 0 or py.max() < 0 or px.min() >= n or py.min() >= n: return
                m = np34.rasterise(px, py, (n, n)); er = m.copy()
                for _ in range(t):
                    e2 = er.copy(); e2[1:, :] &= er[:-1, :]; e2[:-1, :] &= er[1:, :]; e2[:, 1:] &= er[:, :-1]; e2[:, :-1] &= er[:, 1:]; er = e2
                img[m & ~er] = rgb
            for pid in poly:
                if pid in memset: continue
                P = np.array(poly[pid])
                if P[:, 0].max() < x0 or P[:, 1].max() < y0 or P[:, 0].min() > x0 + n or P[:, 1].min() > y0 + n: continue
                outline(pid, (200, 120, 40) if chemo.get(pid, 0) >= THR else (110, 110, 110), 2)
            n_in = (0, 0, 0)
            if tx is not None:
                sel = (tx["x"] >= x0 - 2000) & (tx["x"] < x0 + n + 2000) & (tx["y"] >= y0 - 2000) & (tx["y"] < y0 + n + 2000)
                dens = r29.density(tx["x"][sel], tx["y"][sel], x0, y0, n, UM)
                for lev, rgb in r29.LEVELS:
                    r29.ring(img, dens >= lev, rgb, RING_T)
                win = sel & (tx["x"] >= x0) & (tx["x"] < x0 + n) & (tx["y"] >= y0) & (tx["y"] < y0 + n)
                cls = np.where(tx["cell"] == 0, 2, np.where((tx["fov"] == fov) & np.isin(tx["cell"], ids), 0, 1))
                for k in (1, 2, 0):          # members drawn last so they stay on top
                    for xx, yy in zip(tx["x"][win & (cls == k)] - x0, tx["y"][win & (cls == k)] - y0):
                        dot(img, int(xx), int(yy), (MAGENTA, VIOLET, WHITE)[k])
                n_in = tuple(int((win & (cls == k)).sum()) for k in (0, 1, 2))
            for pid, r in memset.items():
                c = float(r["chemo"]); g = min(c / 100.0, 1.0)
                col = CYAN if r["has_nuc"] == "TRUE" else (255, int(120 + 100 * g), int(60 * (1 - g)))
                outline(pid, col, 4)
                P = np.array(poly[pid]); rp.put_text(img, "%d%s" % (c, "*" if r["seed"] == "TRUE" else ""),
                                                     int(P[:, 0].mean() - x0) - 10, int(P[:, 1].mean() - y0) - 10, col, 20)
            area = sum(float(r["poly_um2"]) for r in mem); nnuc = sum(r["has_nuc"] == "TRUE" for r in mem)
            seam = any(float(r["edge_um"]) < 5 for r in mem); amb = sum(1 for r in mem if r["n_parent"] not in ("", "NA") and int(float(r["n_parent"])) >= 2)
            rp.put_text(img, "%s cid %d  %d members  %.0f um2  %d nuclei  chemo %d%s%s" % (
                SAMP, cid, len(mem), area, nnuc, sum(float(r["chemo"]) for r in mem), "  SEAM" if seam else "", "  %d ambiguous" % amb if amb else ""), 8, 8, WHITE, 22)
            if tx is not None:
                rp.put_text(img, "molecules: %d members / %d other cells / %d unassigned" % n_in, 8, 34, GREY, 20)
            bl = int(round(20 / UM)); img[n - 30:n - 24, n - 20 - bl:n - 20] = 255; rp.put_text(img, "20 um", n - 20 - bl, n - 60, WHITE, 20)
            tiles.append(img)
        store.close()
    if not tiles:
        print("nothing rendered"); return
    T = 520; cols = 4; rows = (len(tiles) + cols - 1) // cols
    sheet = np.zeros((rows * (T + 8) + 8, cols * (T + 8) + 8, 3), np.float32)
    for i, t in enumerate(tiles):
        a8 = np.ascontiguousarray(t.clip(0, 255).round().astype(np.uint8))
        v = pyvips.Image.new_from_memory(a8.tobytes(), a8.shape[1], a8.shape[0], 3, "uchar").resize(T / a8.shape[1])
        th = np.ndarray(buffer=v.write_to_memory(), dtype=np.uint8, shape=(v.height, v.width, 3)).astype(np.float32)
        r, k = divmod(i, cols); sheet[8 + r * (T + 8):8 + r * (T + 8) + th.shape[0], 8 + k * (T + 8):8 + k * (T + 8) + th.shape[1]] = th[:T, :T]
    cap = rp.caption_strip(sheet.shape[1], [
        [("%s, threshold >= %d, %s partition: the %d largest multi-member seed composites by area. Composite: DNA / CD45 / Membrane (x0.5), G omitted." % (
            SAMP, THR, "nucleus-anchored" if PART == "anch" else "connected-components", len(tiles)), WHITE)],
        [("Members outlined thick: cyan = nucleus-bearing member, warm = nucleus-free fragment (yellower = more chemokine); number = CCL17+CCL22 counts, * = seed. Thin orange = other chemokine-positive polygons; thin grey = other cells.", GREY)],
        [("CCL17+CCL22 molecules: magenta = assigned to a member, violet = assigned to another cell, white = unassigned. Rings: molecule density (sigma %g um) at" % r29.SIG_UM, GREY)]
        + [("%g" % lev, rgb) for lev, rgb in r29.LEVELS] + [("per 100 um^2.", GREY)],
        [("Not identified HRS cells: composites of touching chemokine-positive vendor polygons. See output/35_giant_cells/composite_tests.csv for whether contiguity exceeds the permutation null.", AMBER)]])
    out = np.ascontiguousarray(np.concatenate([sheet, cap], 0).clip(0, 255).round().astype(np.uint8))
    fn = os.path.join(OUT, "panels_%s_thr%d_%s.png" % (SAMP, THR, PART))
    pyvips.Image.new_from_memory(out.tobytes(), out.shape[1], out.shape[0], 3, "uchar").pngsave(fn)
    print("wrote", fn, flush=True)


if __name__ == "__main__":
    main()
