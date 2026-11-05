"""Step 37: the hotspot blind spot, on the image.

Two figures for the HRS-candidate report (scripts/35_summary_report.Rmd, section 4):

  hotspots_miss_most_bright_cells.png
      A 1 mm square centred on hotspot IB1_55_162 at pyramid level 2: DAPI in grey with the
      CCL17+CCL22 molecule density (Gaussian sigma 8 um, molecules per 100 um^2) as an
      inferno overlay, the step-11 hotspot bins as cyan boxes and FL neighbourhoods as orange
      ovals. Dozens of isolated bright sources, three hotspot boxes. This replaces the
      scratchpad preview that was first circulated as HRS_hotspots_miss_most_bright_cells.png
      in the project root; same parameters, same colour clamp (21.1), molecule source now the
      step-29 table so it regenerates.

  IB1_57_161_zoom.png
      A 200 um square on hotspot IB1_57_161 (the seam hotspot in the lower left of the 1 mm
      view) at full resolution, two panels. A: the same DAPI + density overlay as above.
      B: the step-35 panel style -- DNA / CD45 / Membrane composite, every vendor polygon
      outlined (warm = CCL17+CCL22 >= 20 counts with the count printed, thin orange = 10-19,
      grey = other), every molecule as a dot coloured by transcript (CCL17 green, CCL22 cyan)
      and the step-29 density contour rings at
      6 / 12 / 24 molecules per 100 um^2.

Reads the IB1 OME-TIFF on scratch, output/29_bright_cells/ib1_chemokine_tx.csv.gz,
output/35_giant_cells/{fov_list.csv, cells_IB1.csv.gz, polygons/IB1_<fov>.csv} and
output/15_imagej/IB1_regions_canvas_px.csv. ~2 min, < 8 GB.

  sbatch --mem=16G scripts/sbatch/run_python.sh scripts/37_hotspot_zoom.py
"""
import csv
import gzip
import importlib.util
import os
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

OUT = "output/37_hotspot_zoom"
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff"
TX = "output/29_bright_cells/ib1_chemokine_tx.csv.gz"
G35 = "output/35_giant_cells"
REGIONS = "output/15_imagej/IB1_regions_canvas_px.csv"
RANGES_TIF = "output/15_imagej/crops/IB1_w00_IB1_nb7269.tif"
UM = rp.UM_PER_PX
FOV_PX = 4256
CLAMP = 21.1                      # density colour clamp, molecules per 100 um^2 (the 200 um preview's 99.5th pct)
OMIT = ("G",)
THR_BRIGHT, THR_DIM = 20, 10
WHITE, GREY, AMBER, CYAN = (255, 255, 255), (200, 200, 200), (255, 200, 60), (0, 229, 255)
GENE_COL = {"CCL17": (0, 255, 90), "CCL22": (0, 229, 255)}      # green / cyan, by transcript
DOT = 3
STOPS = np.array([[0, 0, 4], [87, 16, 110], [188, 55, 84], [249, 142, 9], [252, 255, 164]], np.float32)   # inferno-ish

# the two views, canvas px
MM_SQUARE = dict(id="IB1_55_162", cx=78709, cy=31204, half=4157, level=2)        # 1 mm, as the circulated preview
ZOOM = dict(id="IB1_57_161", cx=77046, cy=32036, um=200.0, level=0)


def rcsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def ramp(t):
    """inferno-ish colour for t in [0, 1], any array shape -> (..., 3)"""
    idx = np.clip(t, 0, 1) * (len(STOPS) - 1)
    i0 = np.floor(idx).astype(int).clip(0, len(STOPS) - 2); f = (idx - i0)[..., None]
    return STOPS[i0] * (1 - f) + STOPS[i0 + 1] * f


def heat_on_dapi(dna, dens):
    nz = dna[dna > 0]; lo, hi = np.percentile(nz, 30), np.percentile(nz, 99)
    gray = (np.clip((dna - lo) / (hi - lo), 0, 1) * 255 * 0.85)[..., None].repeat(3, axis=2)
    t = np.clip(dens / CLAMP, 0, 1); a = (0.15 + 0.75 * t)[..., None]
    return gray * (1 - a) + ramp(t) * a, lo, hi


def colour_bar(img, x, y, w=30, h=400, size=22):
    tt = 1 - np.arange(h) / h
    img[y:y + h, x:x + w] = ramp(tt)[:, None, :]
    rp.put_text(img, "%.1f" % CLAMP, x - 12, y - 36, WHITE, size); rp.put_text(img, "0", x + 8, y + h + 8, WHITE, size)


def scale_bar(img, um, um_per_px, size=24):
    n = img.shape[0]; bl = int(round(um / um_per_px)); bx1, by = img.shape[1] - 30, n - 40
    img[by - 5:by + 5, bx1 - bl - 14:bx1 + 14] *= 0.35; img[by - 3:by + 3, bx1 - bl:bx1] = 255
    rp.put_text(img, "%g um" % um, bx1 - bl, by - 42, WHITE, size)


def draw_regions(img, X0, Y0, n, F, t=3, size=18):
    shown = []
    for o in [r for r in rp.read_csv(REGIONS) if r["group"] in rp.GROUPS]:
        x1, y1, x2, y2 = ((float(o[k]) - (X0 if k[0] == "x" else Y0)) / F for k in ("x1", "y1", "x2", "y2"))
        if x2 > 0 and x1 < n and y2 > 0 and y1 < n:
            (rp.draw_rect if o["shape"] == "rect" else rp.draw_oval)(img, x1, y1, x2, y2, rp.GROUP_COLOR[o["group"]], t)
            rp.put_text(img, rp.short_label(o), int(max(x1, 0)) + 6, int(max(y1, 0)) + 6, rp.GROUP_COLOR[o["group"]], size)
            shown.append(o["id"])
    return shown


def save(img, path, small_w=None):
    a8 = np.ascontiguousarray(img.clip(0, 255).round().astype(np.uint8))
    pyvips.Image.new_from_memory(a8.tobytes(), a8.shape[1], a8.shape[0], 3, "uchar").pngsave(path)
    print("wrote", path, a8.shape, flush=True)


def dot(img, x, y, rgb):
    h, w = img.shape[:2]
    img[max(y - DOT, 0):min(y + DOT + 1, h), max(x - DOT, 0):min(x + DOT + 1, w)] = rgb


def main():
    os.makedirs(OUT, exist_ok=True)
    t = rcsv(TX)
    tx = {"x": np.array([int(float(r["X_px"])) for r in t]), "y": np.array([int(float(r["Y_px"])) for r in t]),
          "fov": np.array([int(float(r["fov"])) for r in t]), "cell": np.array([int(float(r["cell_ID"])) for r in t]),
          "gene": np.array([r["target"] for r in t])}
    del t
    print("molecules:", len(tx["x"]), flush=True)
    ranges = rp.crop_ranges(RANGES_TIF)
    with tifffile.TiffFile(TIFF) as tf:
        store = tf.series[0].aszarr(); grp = zarr.open(store, mode="r")

        # ---- 1 mm square at level 2 -------------------------------------------------------
        v = MM_SQUARE; F = 2 ** v["level"]; um2 = UM * F
        X0, Y0 = v["cx"] - v["half"], v["cy"] - v["half"]
        x0, y0, n = X0 // F, Y0 // F, (2 * v["half"]) // F
        dna = np.asarray(grp[str(v["level"])][4, y0:y0 + n, x0:x0 + n]).astype(np.float32)
        sel = (tx["x"] >= X0 - 1000) & (tx["x"] < X0 + 2 * v["half"] + 1000) & (tx["y"] >= Y0 - 1000) & (tx["y"] < Y0 + 2 * v["half"] + 1000)
        dens = r29.density((tx["x"][sel] - X0) // F + x0, (tx["y"][sel] - Y0) // F + y0, x0, y0, n, um2)
        nin = int(((tx["x"] >= X0) & (tx["x"] < X0 + 2 * v["half"]) & (tx["y"] >= Y0) & (tx["y"] < Y0 + 2 * v["half"])).sum())
        img, lo, hi = heat_on_dapi(dna, dens)
        print("1 mm: DNA L%d %.0f-%.0f, blank %.1f%%; %d molecules; density median %.2f p99 %.2f max %.2f; %.1f%% above clamp" % (
            v["level"], lo, hi, 100 * (dna == 0).mean(), nin, np.median(dens), np.percentile(dens, 99), dens.max(), 100 * (dens > CLAMP).mean()), flush=True)
        shown = draw_regions(img, X0, Y0, n, F)
        print("regions:", " ".join(shown), flush=True)
        scale_bar(img, 200, um2); colour_bar(img, n - 90, 50)
        lines = [[("IB1, 1 mm square centred on hotspot %s: CCL17+CCL22 molecule density on DAPI" % v["id"], CYAN)],
                 [("Density: Gaussian sigma %g um, molecules per 100 um^2, colour clamped at %.1f. %d molecules in the square. Pyramid level %d, %.3f um/px; DAPI %.0f-%.0f (30th/99th pct of non-zero pixels); canvas px origin (%d, %d)." % (
                     r29.SIG_UM, CLAMP, nin, v["level"], um2, lo, hi, X0, Y0), GREY)],
                 [("cyan box = hotspot 100 um bin (c = CCL17+CCL22 cp10k), orange oval = FL neighbourhood.", GREY),
                  ("Raw molecule density, NOT cp10k: it also rises where cells are packed. Hotspots are chemokine-defined, not identified HRS cells.", AMBER)]]
        save(np.concatenate([img, rp.caption_strip(n, lines)], axis=0), os.path.join(OUT, "hotspots_miss_most_bright_cells.png"))

        # ---- 200 um zoom at full resolution, two panels -----------------------------------
        v = ZOOM; n = int(round(v["um"] / UM)); X0, Y0 = v["cx"] - n // 2, v["cy"] - n // 2
        z0 = grp["0"]; arr = np.asarray(z0[:, Y0:Y0 + n, X0:X0 + n])
        sel = (tx["x"] >= X0 - 2000) & (tx["x"] < X0 + n + 2000) & (tx["y"] >= Y0 - 2000) & (tx["y"] < Y0 + n + 2000)
        dens = r29.density(tx["x"][sel], tx["y"][sel], X0, Y0, n, UM)
        win = sel & (tx["x"] >= X0) & (tx["x"] < X0 + n) & (tx["y"] >= Y0) & (tx["y"] < Y0 + n)
        # A: DAPI + density
        A, lo, hi = heat_on_dapi(arr[4].astype(np.float32), dens)
        draw_regions(A, X0, Y0, n, 1, t=5, size=22)
        scale_bar(A, 50, UM); colour_bar(A, n - 90, 50)
        rp.put_text(A, "A", 8, 8, WHITE, 36)
        # B: composite + polygons + molecules + rings
        B = rp.composite(arr, ranges, OMIT); np.clip(B, 0, 255, out=B)
        fl = {int(r["fov"]): (int(float(r["x0"])), int(float(r["y0"]))) for r in rcsv(os.path.join(G35, "fov_list.csv")) if r["samp"] == "IB1"}
        fovs = [f for f, (fx, fy) in fl.items() if fx < X0 + n and fx + FOV_PX > X0 and fy < Y0 + n and fy + FOV_PX > Y0]
        chemo = {}
        with gzip.open(os.path.join(G35, "cells_IB1.csv.gz"), "rt", newline="") as fh:
            for r in csv.DictReader(fh):
                if int(r["fov"]) in fovs: chemo[(int(r["fov"]), int(r["cell_ID"]))] = float(r["chemo"])
        poly = {}
        for f in fovs:
            fx, fy = fl[f]
            with open(os.path.join(G35, "polygons", "IB1_%d.csv" % f), newline="") as fh:
                for r in csv.reader(fh):
                    poly.setdefault((f, int(r[1])), []).append((fx + FOV_PX - 1 - float(r[2]) - X0, fy + float(r[3]) - Y0))
        def outline(P, rgb, t):
            px, py = P[:, 0], P[:, 1]
            m = np34.rasterise(px, py, (n, n)); er = m.copy()
            for _ in range(t):
                e2 = er.copy(); e2[1:, :] &= er[:-1, :]; e2[:-1, :] &= er[1:, :]; e2[:, 1:] &= er[:, :-1]; e2[:, :-1] &= er[:, 1:]; er = e2
            B[m & ~er] = rgb
        inwin = {k: np.array(P) for k, P in poly.items() if np.array(P)[:, 0].max() >= 0 and np.array(P)[:, 1].max() >= 0 and np.array(P)[:, 0].min() < n and np.array(P)[:, 1].min() < n}
        bright = [k for k in inwin if chemo.get(k, 0) >= THR_BRIGHT]
        for k, P in inwin.items():
            c = chemo.get(k, 0)
            if c >= THR_BRIGHT: continue
            outline(P, (200, 120, 40) if c >= THR_DIM else (110, 110, 110), 2)
        for lev, rgb in r29.LEVELS:
            r29.ring(B, dens >= lev, rgb, 3)
        bright_set = set(bright)
        for gname, rgb in GENE_COL.items():
            m = win & (tx["gene"] == gname)
            for xx, yy in zip(tx["x"][m] - X0, tx["y"][m] - Y0):
                dot(B, int(xx), int(yy), rgb)
        n_in = tuple(int((win & (tx["gene"] == gname)).sum()) for gname in GENE_COL)
        for k in bright:
            P = inwin[k]; c = chemo[k]; g = min(c / 100.0, 1.0); col = (255, int(120 + 100 * g), int(60 * (1 - g)))
            outline(P, col, 4)
            rp.put_text(B, "%d" % c, int(P[:, 0].mean()) - 10, int(P[:, 1].mean()) - 10, col, 20)
        draw_regions(B, X0, Y0, n, 1, t=5, size=22)
        scale_bar(B, 50, UM)
        rp.put_text(B, "B", 8, 8, WHITE, 36)
        rp.put_text(B, "molecules: %d CCL17 / %d CCL22" % n_in, 60, 14, GREY, 20)
        # what the hotspot bin actually holds: molecules inside its rect, and the share in its >= 20 polygons
        reg = [r for r in rp.read_csv(REGIONS) if r["id"] == ZOOM["id"]][0]
        rx1, ry1, rx2, ry2 = (float(reg[k]) for k in ("x1", "y1", "x2", "y2"))
        inbox = (tx["x"] >= rx1) & (tx["x"] <= rx2) & (tx["y"] >= ry1) & (tx["y"] <= ry2)
        n_box = int(inbox.sum())
        per_poly = sorted(((int((inbox & (tx["fov"] == f) & (tx["cell"] == c)).sum()), f, c) for f, c in bright), reverse=True)
        top = [p for p in per_poly if p[0] > 0]
        box_bright = sum(p[0] for p in top)
        print("zoom: FOVs %s, %d polygons in window, %d with >= %d counts; molecules %s; DNA %.0f-%.0f" % (fovs, len(inwin), len(bright), THR_BRIGHT, n_in, lo, hi), flush=True)
        print("hotspot box: %d molecules; %d (%.0f%%) in its %d polygons with >= %d counts: %s" % (
            n_box, box_bright, 100.0 * box_bright / max(n_box, 1), len(top), THR_BRIGHT, ", ".join("c_2_%d_%d=%d" % (f, c, k) for k, f, c in top)), flush=True)
        with open(os.path.join(OUT, "%s_box_molecules.csv" % ZOOM["id"]), "w") as fh:
            fh.write("hotspot,n_molecules_in_box,n_in_bright_polygons,n_bright_polygons,thr,bright_polygons\n")
            fh.write("%s,%d,%d,%d,%d,%s\n" % (ZOOM["id"], n_box, box_bright, len(top), THR_BRIGHT, ";".join("c_2_%d_%d=%d" % (f, c, k) for k, f, c in top)))
        store.close()
    gap = np.zeros((n, 24, 3), np.float32)
    top = np.concatenate([A, gap, B], axis=1)
    lines = [[("Hotspot %s, %g um square at full resolution (%.3f um/px), canvas px origin (%d, %d). FOV seam runs through the window." % (ZOOM["id"], ZOOM["um"], UM, X0, Y0), CYAN)],
             [("A: DAPI + CCL17+CCL22 molecule density (sigma %g um, per 100 um^2, clamped at %.1f) -- the 1 mm view's scale. Cyan box = the hotspot's 100 um bin." % (r29.SIG_UM, CLAMP), GREY)],
             [("B: DNA / CD45 / Membrane (x0.5), G omitted. Polygons: warm thick = CCL17+CCL22 >= %d counts (yellower = more; number = counts), thin orange = %d-%d, grey = other." % (THR_BRIGHT, THR_DIM, THR_BRIGHT - 1), GREY)],
             [("Molecules:", GREY), ("CCL17 green", GENE_COL["CCL17"]), ("CCL22 cyan", GENE_COL["CCL22"]), ("(the cyan box is the hotspot bin, not a molecule). Rings: density at", GREY)]
             + [("%g" % lev, rgb) for lev, rgb in r29.LEVELS] + [("per 100 um^2.", GREY)],
             [("The box holds %d molecules, %d (%.0f%%) of them in its %d polygons with >= %d counts (three touching: one step-35 composite); the bin's cp10k pools them with ~40 other cells. Not identified HRS cells." % (
                 n_box, box_bright, 100.0 * box_bright / max(n_box, 1), len(top), THR_BRIGHT), AMBER)]]
    save(np.concatenate([top, rp.caption_strip(top.shape[1], lines)], axis=0), os.path.join(OUT, "%s_zoom.png" % ZOOM["id"]))


if __name__ == "__main__":
    main()
