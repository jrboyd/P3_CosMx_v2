"""Step 34 QA: DNA windows with the vendor polygon (cyan) and the polygon-bounded nuclear
mask (red) for a random set of candidates and controls, one contact sheet each, so the
measurements in nuclei_in_polygon.csv can be checked against the picture.

  python scripts/34_polygon_qa_sheet.py output/34_nuclear_size [n_per_group] [random|pale]
    pale: the n cells per group with the LARGEST pale-not-dense fraction (pale_frac - nuc_frac),
          i.e. polygons filled with DNA at intermediate intensity -- the phenotype to check
"""
import csv
import os
import random
import sys
from collections import defaultdict

import numpy as np
import pyvips
import tifffile
import zarr

sys.path.insert(0, "scripts")
import importlib.util
spec = importlib.util.spec_from_file_location("np34", "scripts/34_nuclear_in_polygon.py")
np34 = importlib.util.module_from_spec(spec); spec.loader.exec_module(np34)
spec2 = importlib.util.spec_from_file_location("rp", "scripts/15_roi_pngs.py")
rp = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(rp)

OUT = sys.argv[1] if len(sys.argv) > 1 else "output/34_nuclear_size"
N_PER = int(sys.argv[2]) if len(sys.argv) > 2 else 20
MODE = sys.argv[3] if len(sys.argv) > 3 else "random"
random.seed(34)


def main():
    meas = {r["cell_id"]: r for r in csv.DictReader(open(os.path.join(OUT, "nuclei_in_polygon.csv")))}
    cells = {r["cell_id"]: r for r in csv.DictReader(open(np34.CELLS))}
    poly = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUT, "polygons_canvas.csv"))):
        poly[r["cell_id"]].append((float(r["X_px"]), float(r["Y_px"])))
    R = int(round(np34.WIN_UM / 2 / np34.UM)); n = 2 * R + 1
    with tifffile.TiffFile(np34.TIFF) as tf:
        store = tf.series[0].aszarr(); z0 = zarr.open(store, mode="r")["0"]
        for g in ("candidate", "control"):
            ids = [c for c in meas if meas[c]["group"] == g and meas[c]["ok"] == "True"]
            if MODE == "pale":
                ids = sorted(ids, key=lambda c: -(float(meas[c]["pale_frac"]) - float(meas[c]["nuc_frac"])))[:N_PER]
            else:
                random.shuffle(ids); ids = sorted(ids[:N_PER], key=lambda c: -float(meas[c]["nuc_area"]))
            tiles = []
            for c in ids:
                m = meas[c]; cx, cy = int(round(float(cells[c]["X_px"]))), int(round(float(cells[c]["Y_px"])))
                x0, y0 = cx - R, cy - R
                a = np.asarray(z0[np34.DNA, y0:y0 + n, x0:x0 + n]).astype(np.float32)
                s = np34.blur(a, np34.SIG_UM / np34.UM)
                thr = float(m["nuc_thr"])
                P = np.array(poly[c]); inside = np34.rasterise(P[:, 0] - x0, P[:, 1] - y0, s.shape)
                nuc = inside & (s >= thr)
                v = np.clip((a - np.percentile(a, 1)) / (np.percentile(a, 99.5) - np.percentile(a, 1) + 1e-6), 0, 1) * 255
                rgb = np.repeat(v[..., None], 3, 2).astype(np.float32)
                def edge(mk):
                    er = mk.copy(); er[1:, :] &= mk[:-1, :]; er[:-1, :] &= mk[1:, :]; er[:, 1:] &= mk[:, :-1]; er[:, :-1] &= mk[:, 1:]
                    return mk & ~er
                e_in = edge(inside); e_in |= np.roll(e_in, 1, 0) | np.roll(e_in, 1, 1)
                rgb[e_in] = (0, 229, 255)
                rgb[edge(nuc) & ~e_in] = (255, 60, 60)
                rgb[R - 2:R + 3, R - 2:R + 3] = (60, 255, 60)
                rp.put_text(rgb, "%s" % c.replace("c_2_", ""), 6, 6, (255, 255, 255), 18)
                rp.put_text(rgb, "poly %.0f  dense %.0f%%  pale %.0f%%  rel %.2f  cv %s  holes %s" % (
                    float(m["poly_area"]), 100 * float(m["nuc_frac"]), 100 * float(m["pale_frac"]), float(m["dna_rel"]),
                    m["pale_cv"] or "-", m["pale_n_holes"] or "-"), 6, n - 30, (255, 220, 120), 16)
                tiles.append(rgb)
            cols = 5; rows = (len(tiles) + cols - 1) // cols
            sheet = np.zeros((rows * (n + 6), cols * (n + 6), 3), np.float32)
            for i, t in enumerate(tiles):
                r, k = divmod(i, cols); sheet[r * (n + 6):r * (n + 6) + n, k * (n + 6):k * (n + 6) + n] = t
            cap = rp.caption_strip(sheet.shape[1], [
                [("%s cells, 40 um DNA windows, full resolution. Cyan = vendor cell polygon; red = DNA >= window Otsu inside it (the nuclear mask); green = recorded centroid." % g, (255, 255, 255))],
                [("Label: polygon area um2; dense = DNA >= Otsu (%% of polygon); pale = DNA >= background + 25%% of (Otsu minus background); rel = polygon median DNA on a 0 (background) to 1 (Otsu) scale; cv of DNA in the pale mask; dark holes in it. Selection: %s." % MODE, (200, 200, 200))]])
            out = np.concatenate([sheet, cap], 0).clip(0, 255).round().astype(np.uint8)
            out = np.ascontiguousarray(out)
            fn = os.path.join(OUT, "qa_polygons_%s%s.png" % (g, "" if MODE == "random" else "_" + MODE))
            pyvips.Image.new_from_memory(out.tobytes(), out.shape[1], out.shape[0], 3, "uchar").pngsave(fn)
            print("wrote", fn, flush=True)
        store.close()


if __name__ == "__main__":
    main()
