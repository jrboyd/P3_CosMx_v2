"""Step 35 review panels: the largest multi-member seed composites on morphology.

Reads <OUT>/polygon_assignments.csv.gz, seed_composites.csv.gz, fov_list.csv and the
per-FOV polygon files, and renders for each selected composite a full-resolution window
(DNA / CD45 / Membrane composite as in the step-15 ROI PNGs, G omitted, Membrane x0.5)
with every polygon in the window outlined: members of the composite in chemokine-graded
warm colours (the nucleus-bearing member in cyan), other chemokine-positive polygons in
dim orange, the rest grey. A caption carries the composite's members, area, nuclei,
chemokine total and seam flag. One contact sheet per group.

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

OUT, SAMP, THR, PART = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
N_TOP = int(sys.argv[5]) if len(sys.argv) > 5 else 24
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/%s.ome.tiff" % SAMP
UM = rp.UM_PER_PX
FOV_PX = 4256
PAD_UM = 25.0
OMIT = ("G",)
WHITE, GREY, AMBER, CYAN = (255, 255, 255), (200, 200, 200), (255, 200, 60), (0, 229, 255)


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
            bl = int(round(20 / UM)); img[n - 30:n - 24, n - 20 - bl:n - 20] = 255; rp.put_text(img, "20 um", n - 20 - bl, n - 60, WHITE, 20)
            tiles.append(img)
        store.close()
    if not tiles:
        print("nothing rendered"); return
    T = 420; cols = 4; rows = (len(tiles) + cols - 1) // cols
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
        [("Not identified HRS cells: composites of touching chemokine-positive vendor polygons. See output/35_giant_cells/composite_tests.csv for whether contiguity exceeds the permutation null.", AMBER)]])
    out = np.ascontiguousarray(np.concatenate([sheet, cap], 0).clip(0, 255).round().astype(np.uint8))
    fn = os.path.join(OUT, "panels_%s_thr%d_%s.png" % (SAMP, THR, PART))
    pyvips.Image.new_from_memory(out.tobytes(), out.shape[1], out.shape[0], 3, "uchar").pngsave(fn)
    print("wrote", fn, flush=True)


if __name__ == "__main__":
    main()
