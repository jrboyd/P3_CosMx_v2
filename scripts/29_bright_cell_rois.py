"""Step 29 renderer -- one annotated full-resolution PNG per bright chemokine cell.

Reads output/29_bright_cells/{roi_candidates.csv, ib1_chemokine_tx.csv.gz} (written by
scripts/29_script_bright_cells.R) and the IB1 OME-TIFF on scratch, and writes
output/29_bright_cells/roi_png/<group>/NN_<cell_id>.png, a contact sheet per group, and
roi_png_index.csv.

STYLE is the step-15 ROI PNGs' (scripts/15_roi_pngs.py is imported, not copied): same
composite, same fixed display ranges as the crop TIFFs, Membrane at its GAIN, channel G
omitted (marker unknown, CLAUDE.md), hotspot boxes / FL ovals with the same hues, the
same wrapping caption strip. What is new:

  CONTOURS  CCL17+CCL22 transcript density (Gaussian sigma 8 um, molecules per 100 um^2,
            the measure of the step-15 field previews) drawn as rings at LEVELS. With that
            kernel a single cell's own peak is roughly count / (2 pi sigma^2) x 100, so the
            6 / 12 / 24 rings sit around cells of roughly 25 / 50 / 100 counts: the number
            of rings around a cell reads as its brightness.
  CROSSHAIR marks the selected cell's centroid.

Each group takes candidates in order and keeps the first 20 whose 150 um window sits on
imaged tissue (<= 5% of pixels blank in every channel -- some IB1 FOVs have transcripts
but no morphology, CLAUDE.md) and inside the canvas.

  python scripts/29_bright_cell_rois.py [--smoke]   # --smoke: 2 per group, reads/writes the
                                                  # output/29_bright_cells_smoke/ of the R smoke run
"""
import csv
import gzip
import importlib.util
import os
import sys

import numpy as np
import pyvips
import tifffile
import zarr

spec = importlib.util.spec_from_file_location("rp", "scripts/15_roi_pngs.py")
rp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rp)

SMOKE = "--smoke" in sys.argv
OUT = "output/29_bright_cells_smoke" if SMOKE else "output/29_bright_cells"
PNG_DIR = os.path.join(OUT, "roi_png")
TIFF = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs/IB1.ome.tiff"
CROP_FOR_RANGES = "output/15_imagej/crops/IB1_w00_IB1_nb7269.tif"
N_PER_GROUP = 2 if SMOKE else 20
WIN_UM, MARGIN_UM, SIG_UM = 150.0, 100.0, 8.0
MAX_BLANK = 0.05
OMIT = ("G",)
LEVELS = [(6.0, (188, 55, 84)), (12.0, (249, 142, 9)), (24.0, (252, 255, 164))]  # inferno stops
RING_PX = 4

GROUPS = [
    ("high_random", "High chemokine (100+ counts) - random draw"),
    ("high_hrs", "High chemokine (100+ counts) - top HRS-leaning score"),
    ("medium_random", "Medium chemokine (50-99 counts) - random draw"),
    ("medium_hrs", "Medium chemokine (50-99 counts) - top HRS-leaning score"),
    ("low_random", "Low chemokine (10-49 counts) - random draw"),
    ("low_hrs", "Low chemokine (10-49 counts) - top HRS-leaning score"),
]
WHITE, GREY, AMBER = (255, 255, 255), (200, 200, 200), (255, 200, 60)


def log(m):
    print(m, flush=True)


def shift(m, dy, dx, fill):
    h, w = m.shape
    out = np.full_like(m, fill)
    out[max(dy, 0):h + min(dy, 0), max(dx, 0):w + min(dx, 0)] = \
        m[max(-dy, 0):h + min(-dy, 0), max(-dx, 0):w + min(-dx, 0)]
    return out


def erode(m, r):
    """Binary erosion by a disk; the image border counts as inside, so a region cut by
    the window edge gets no false ring along the edge."""
    out = m.copy()
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            if dy * dy + dx * dx <= r * r and (dy or dx):
                out &= shift(m, dy, dx, True)
    return out


def ring(img, mask, rgb, t):
    edge = mask & ~erode(mask, t)
    halo = ~erode(~edge, 2) & ~edge                  # 2 px dark halo for contrast
    img[halo] *= 0.25
    img[edge] = rgb


def crosshair(img, cx, cy, rgb=WHITE, r1=55, r2=110, t=5):
    h, w = img.shape[:2]
    m = np.zeros((h, w), bool)
    cx, cy = int(round(cx)), int(round(cy))
    for a0, a1 in ((cx + r1, cx + r2), (cx - r2, cx - r1)):
        m[max(cy - t // 2, 0):max(cy + t // 2 + 1, 0), max(a0, 0):max(min(a1, w), 0)] = True
    for a0, a1 in ((cy + r1, cy + r2), (cy - r2, cy - r1)):
        m[max(a0, 0):max(min(a1, h), 0), max(cx - t // 2, 0):max(cx + t // 2 + 1, 0)] = True
    img[~erode(~m, 2) & ~m] *= 0.25
    img[m] = rgb


def density(tx_x, tx_y, x0, y0, n, um_per_px):
    """Gaussian KDE of molecule positions on the window grid, molecules per 100 um^2."""
    M = int(round(MARGIN_UM / um_per_px))
    g = n + 2 * M
    gx = (tx_x - x0).astype(int) + M
    gy = (tx_y - y0).astype(int) + M
    ok = (gx >= 0) & (gx < g) & (gy >= 0) & (gy < g)
    grid = np.zeros((g, g), np.float32)
    np.add.at(grid, (gy[ok], gx[ok]), 1.0)
    sig = SIG_UM / um_per_px
    ky = np.fft.fftfreq(g)[:, None]
    kx = np.fft.rfftfreq(g)[None, :]
    k = np.exp(-2 * (np.pi * sig) ** 2 * (kx ** 2 + ky ** 2))
    dens = np.fft.irfft2(np.fft.rfft2(grid) * k, s=(g, g))
    return dens[M:M + n, M:M + n] / (um_per_px ** 2) * 100.0


def read_csv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh))


def main():
    os.makedirs(PNG_DIR, exist_ok=True)
    um = rp.UM_PER_PX
    n = int(round(WIN_UM / um))
    ranges = rp.crop_ranges(CROP_FOR_RANGES)
    cand = read_csv(os.path.join(OUT, "roi_candidates.csv"))
    tx = read_csv(os.path.join(OUT, "ib1_chemokine_tx.csv.gz"))
    TX = np.array([float(r["X_px"]) for r in tx]); TY = np.array([float(r["Y_px"]) for r in tx])
    # Molecules ASSIGNED to each cell. A selected cell's own count must equal its
    # CCL17+CCL22 from the expression matrix; a mismatch (e.g. the transcript table
    # missing its FOV, as in an early smoke run) would otherwise just draw no rings.
    own = {}
    for r in tx:
        if r["cell_ID"] != "0":
            k = (r["fov"], r["cell_ID"])
            own[k] = own.get(k, 0) + 1
    del tx
    log("molecules loaded: %d" % len(TX))
    regions = [r for r in rp.read_csv("output/15_imagej/IB1_regions_canvas_px.csv")
               if r["group"] in rp.GROUPS]
    for r in regions:
        for k in ("x1", "y1", "x2", "y2"):
            r[k] = float(r[k])

    index = []
    with tifffile.TiffFile(TIFF) as tf:
        store = tf.series[0].aszarr()
        z0 = zarr.open(store, mode="r")["0"]
        _, H, W = z0.shape
        for gkey, gtitle in GROUPS:
            gdir = os.path.join(PNG_DIR, gkey)
            os.makedirs(gdir, exist_ok=True)
            rows = sorted([c for c in cand if c["group"] == gkey], key=lambda c: int(c["rank"]))
            kept, rejected = [], {"blank": 0, "canvas": 0}
            for c in rows:
                if len(kept) == N_PER_GROUP:
                    break
                cx, cy = float(c["X_px"]), float(c["Y_px"])
                x0, y0 = int(round(cx)) - n // 2, int(round(cy)) - n // 2
                if x0 < 0 or y0 < 0 or x0 + n > W or y0 + n > H:
                    rejected["canvas"] += 1
                    continue
                arr = np.asarray(z0[:, y0:y0 + n, x0:x0 + n])
                blank = float((arr.max(axis=0) == 0).mean())
                if blank > MAX_BLANK:
                    rejected["blank"] += 1
                    continue
                img = rp.composite(arr, ranges, OMIT)
                del arr
                np.clip(img, 0, 255, out=img)

                own_n = own.get((c["fov"], c["cell_ID"]), 0)
                if own_n != int(c["chemo"]):
                    log("  WARNING %s: %d assigned molecules in the transcript table vs chemo %s"
                        % (c["cell_id"], own_n, c["chemo"]))
                sel = (TX >= x0 - 900) & (TX < x0 + n + 900) & (TY >= y0 - 900) & (TY < y0 + n + 900)
                dens = density(TX[sel], TY[sel], x0, y0, n, um)
                for lev, rgb in LEVELS:
                    ring(img, dens >= lev, rgb, RING_PX)

                for o in regions:
                    if o["x2"] > x0 and o["x1"] < x0 + n and o["y2"] > y0 and o["y1"] < y0 + n:
                        box = (o["x1"] - x0, o["y1"] - y0, o["x2"] - x0, o["y2"] - y0)
                        (rp.draw_rect if o["shape"] == "rect" else rp.draw_oval)(
                            img, *box, rp.GROUP_COLOR[o["group"]], rp.LINE_PX)
                        if o["shape"] == "rect":
                            lx, ly = int(max(box[0], 0)) + 15, int(max(box[1], 0)) + 15
                        else:
                            lx = int(max(box[0], 8))
                            ly = int(box[1]) - 44 if box[1] >= 50 else int(box[3]) + 12
                        rp.put_text(img, rp.short_label(o), lx, ly, rp.GROUP_COLOR[o["group"]], 24)
                crosshair(img, cx - x0, cy - y0)

                bl = int(round(rp.SCALEBAR_UM / um))
                bx1, by = n - 40, n - 50
                img[by - 6:by + 6, bx1 - bl - 20:bx1 + 20] *= 0.35
                img[by - 4:by + 4, bx1 - bl:bx1] = 255
                rp.put_text(img, "%d µm" % rp.SCALEBAR_UM, bx1 - bl, by - 52, WHITE, 26)

                slot = len(kept) + 1
                lines = [
                    [("%s  #%d/%d" % (gtitle, slot, N_PER_GROUP), WHITE)],
                    [("%s   CCL17+CCL22 %s counts (CCL17 %s, CCL22 %s)   area %.0f µm²   depth %s"
                      % (c["cell_id"], c["chemo"], c["CCL17"], c["CCL22"], float(c["Area.um2"]),
                         c["depth0"]), WHITE)],
                    [("IL13 %s  BATF3 %s  CD30 %s  (%s of 3 detected)   HRS-leaning score %+.2f   IRF4 %s  LTA %s  CCR7 %s  LAMP3 %s"
                      % (c["IL13"], c["BATF3"], c["TNFRSF8"], c["k3"], float(c["hrs_score"]),
                         c["IRF4"], c["LTA"], c["CCR7"], c["LAMP3"]), (150, 190, 255))],
                    [("Channels:", GREY)]
                    + [("%s %.0f-%.0f%s" % (nm, lo, hi, " x%g" % rp.GAIN[nm] if nm in rp.GAIN else ""), rgb)
                       for (nm, rgb), (lo, hi) in zip(rp.CHANNELS, ranges) if nm not in OMIT]
                    + [("(G omitted, marker unknown)", GREY)],
                    [("Rings: CCL17+CCL22 molecule density (sigma %g µm) at" % SIG_UM, GREY)]
                    + [("%g" % lev, rgb) for lev, rgb in LEVELS]
                    + [("per 100 µm² - about a 25 / 50 / 100-count cell's own peak.", GREY),
                       ("Crosshair = selected cell.", WHITE),
                       ("cyan box = hotspot bin, orange oval = FL neighbourhood.", GREY)],
                    [("Full resolution %.4f µm/px, %g µm window, canvas px origin (%d, %d). A candidate ranked on RNA, NOT an identified HRS cell; transcripts leak between touching cells."
                      % (um, WIN_UM, x0, y0), AMBER)],
                ]
                out_img = np.concatenate([img, rp.caption_strip(n, lines)], axis=0)
                png = os.path.join(gdir, "%02d_%s.png" % (slot, c["cell_id"]))
                a8 = np.ascontiguousarray(out_img.clip(0, 255).round().astype(np.uint8))
                pyvips.Image.new_from_memory(a8.tobytes(), a8.shape[1], a8.shape[0], 3,
                                             "uchar").pngsave(png, compression=6)
                kept.append((slot, c, png, img))
                index.append({"group": gkey, "slot": slot, "candidate_rank": c["rank"],
                              "cell_id": c["cell_id"], "png": os.path.relpath(png, PNG_DIR),
                              "chemo": c["chemo"], "hrs_score": c["hrs_score"], "k3": c["k3"],
                              "IL13": c["IL13"], "BATF3": c["BATF3"], "TNFRSF8": c["TNFRSF8"],
                              "in_hotspot": c["hot"], "own_tx": own_n,
                              "own_tx_matches": own_n == int(c["chemo"]),
                              "blank_frac": round(blank, 4),
                              "canvas_x0": x0, "canvas_y0": y0,
                              "peak_density": round(float(dens.max()), 2)})
            log("%-14s kept %d of %d candidates tried (blank %d, off-canvas %d)"
                % (gkey, len(kept), len(kept) + sum(rejected.values()),
                   rejected["blank"], rejected["canvas"]))
            contact_sheet(gkey, gtitle, kept)
            kept.clear()
        store.close()

    with open(os.path.join(PNG_DIR, "roi_png_index.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(index[0]))
        wr.writeheader()
        wr.writerows(index)
    log("done: %d PNGs in %s" % (len(index), PNG_DIR))


def contact_sheet(gkey, gtitle, kept, thumb=360, cols=5):
    """All of a group's ROIs on one page, image area only, with a one-line label each."""
    if not kept:
        return
    rows = (len(kept) + cols - 1) // cols
    lab_h, pad, head = 34, 8, 60
    sheet = np.zeros((head + rows * (thumb + lab_h + pad), cols * (thumb + pad) + pad, 3), np.float32)
    rp.put_text(sheet, gtitle + "   (label: CCL17+CCL22 counts, HRS-leaning score, IL13/BATF3/CD30 detected)",
                pad, 14, WHITE, 24, box=False)
    for i, (slot, c, _, img) in enumerate(kept):
        r, k = divmod(i, cols)
        a8 = np.ascontiguousarray(img.clip(0, 255).round().astype(np.uint8))
        v = pyvips.Image.new_from_memory(a8.tobytes(), a8.shape[1], a8.shape[0], 3, "uchar")
        v = v.resize(thumb / v.width)
        t = np.ndarray(buffer=v.write_to_memory(), dtype=np.uint8,
                       shape=(v.height, v.width, 3)).astype(np.float32)
        y, x = head + r * (thumb + lab_h + pad), pad + k * (thumb + pad)
        sheet[y:y + t.shape[0], x:x + t.shape[1]] = t[:thumb, :thumb]
        rp.put_text(sheet, "#%d  %s ct  score %+.2f  %s/3" % (slot, c["chemo"], float(c["hrs_score"]), c["k3"]),
                    x, y + thumb + 4, GREY, 18, box=False)
    a8 = np.ascontiguousarray(sheet.clip(0, 255).round().astype(np.uint8))
    pyvips.Image.new_from_memory(a8.tobytes(), a8.shape[1], a8.shape[0], 3, "uchar").pngsave(
        os.path.join(PNG_DIR, "%s_contact_sheet.png" % gkey))


if __name__ == "__main__":
    main()
