"""One full-resolution annotated PNG per ImageJ-export region.

The alternative to opening the crop TIFFs in Fiji: a flat RGB composite per
region (hodgkin_hotspot + fl_nb, the 45 regions the crops were built for), with
every region that falls in the window outlined and labelled, and a caption strip
underneath carrying what a reader needs without the repo.

Inputs, all from scripts/15_submit_imagej_export.sh:
  output/15_imagej/IB1_regions_canvas_px.csv     region bounds, full-res canvas px
  output/15_imagej/IB1_regions_L5.csv            per-region blank_frac
  output/15_imagej/crops/IB1_crop_manifest.csv   crop window origins
  output/15_imagej/crops/IB1_w##_<seed>.tif      level-0 pixels (5 ch, uint16)

Writes output/15_imagej/roi_png/<id>.png and roi_png_index.csv.

PIXELS. Read from the crop TIFFs, not the OME on scratch: they are level-0
extracts, so nothing is resampled, and this runs even after scratch is purged.
Crop pixels are canvas px in image orientation and the regions are canvas px too,
so local = canvas - crop origin, with no axis flip anywhere (the flips in
CLAUDE.md's axis rule apply to slide mm, which nothing here touches).

CONTRAST. Fixed: the display ranges embedded in the crop TIFFs (one shared set,
30th/99th pct of non-zero level-0 pixels, see 15_imagej_tiles.py). Read from the
first crop's IJ metadata and asserted identical in every crop used, so every PNG
is on the same scale and matches what Fiji shows, except that GAIN (below) dims
Membrane to half. Channels are blended
additively with clipping, as Fiji's composite mode does.

CHANNEL G IS OMITTED by default: its marker is unknown and it carries no signal
of its own (CLAUDE.md, scripts/15_script_channel_identity.R). --with-G puts it back.

  python scripts/15_roi_pngs.py [SAMPLE] [--smoke] [--with-G]   # --smoke: 2 regions
"""
import csv
import os
import sys

import numpy as np
import pyvips
import tifffile

OUT_ROOT = "output/15_imagej"
UM_PER_PX = 0.120280945
PAD_UM = 50.0                       # context around each region, every side
LINE_PX = 5                         # outline thickness (~0.6 um)
SCALEBAR_UM = 50.0
GROUPS = ("hodgkin_hotspot", "fl_nb")

# Positional, matching cosmx_channel_names() / cosmx_channel_colors()
# (scripts/lib/cosmx_tiff.R) and CHANNELS in 15_imagej_tiles.py.
CHANNELS = [("PanCK", (0, 0, 255)),
            ("G", (0, 255, 0)),
            ("Membrane", (255, 255, 0)),
            ("CD45", (255, 0, 0)),
            ("DNA", (238, 130, 238))]

# Display gain applied after the range stretch (1.0 = as in Fiji). Membrane is
# halved because at full weight its stromal staining swamps the other channels.
GAIN = {"Membrane": 0.5}

# Same outline hues as the ImageJ ROI sets (15_imagej_tiles.py GROUP_COLOR).
GROUP_COLOR = {"hodgkin_hotspot": (0, 229, 255), "fl_nb": (235, 104, 52)}
GROUP_NAME = {"hodgkin_hotspot": "Hodgkin chemokine hotspot",
              "fl_nb": "FL kappa-restricted neighbourhood"}


def log(msg):
    print(msg, flush=True)


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def crop_ranges(path):
    """The (lo, hi) per channel that 15_imagej_tiles.py embedded."""
    with tifffile.TiffFile(path) as t:
        flat = t.imagej_metadata["Ranges"]
    return [(float(flat[2 * i]), float(flat[2 * i + 1])) for i in range(len(CHANNELS))]


def short_label(r):
    if r["group"] == "hodgkin_hotspot":
        s = "%s  c%.1f" % (r["id"], float(r["metric"]))
    else:
        s = "%s  k%.2f" % (r["id"], float(r["metric"]))
    return s + ("  SEAM" if r["flag"] == "seam" else "")


def composite(arr, ranges, omit=()):
    """(C, Y, X) uint16 -> (Y, X, 3) float32 in 0..255, additive, unclipped."""
    out = np.zeros(arr.shape[1:] + (3,), np.float32)
    for c, (name, rgb) in enumerate(CHANNELS):
        if name in omit:
            continue
        lo, hi = ranges[c]
        v = np.clip((arr[c].astype(np.float32) - lo) / (hi - lo), 0.0, 1.0) * GAIN.get(name, 1.0)
        for k in range(3):
            if rgb[k]:
                out[..., k] += v * rgb[k]
    return out


def paint(img, mask, rgb, alpha=1.0):
    """Blend a solid colour into img where mask (bool or 0..1 float) is set."""
    a = mask.astype(np.float32)[..., None] * alpha
    img[:] = img * (1.0 - a) + np.asarray(rgb, np.float32) * a


def draw_rect(img, x1, y1, x2, y2, rgb, t):
    h, w = img.shape[:2]
    m = np.zeros((h, w), bool)
    xa, xb = int(round(x1)), int(round(x2))
    ya, yb = int(round(y1)), int(round(y2))
    cx = slice(max(xa, 0), min(xb, w))
    cy = slice(max(ya, 0), min(yb, h))
    for yy in (ya, yb - t):
        m[max(yy, 0):max(min(yy + t, h), 0), cx] = True
    for xx in (xa, xb - t):
        m[cy, max(xx, 0):max(min(xx + t, w), 0)] = True
    paint(img, m, rgb)


def draw_oval(img, x1, y1, x2, y2, rgb, t):
    h, w = img.shape[:2]
    cx, cy, rx, ry = (x1 + x2) / 2, (y1 + y2) / 2, (x2 - x1) / 2, (y2 - y1) / 2
    xa, xb = max(int(x1) - t, 0), min(int(x2) + t + 1, w)
    ya, yb = max(int(y1) - t, 0), min(int(y2) + t + 1, h)
    if xa >= xb or ya >= yb:
        return
    yy, xx = np.mgrid[ya:yb, xa:xb].astype(np.float32)
    r = np.sqrt(((xx - cx) / rx) ** 2 + ((yy - cy) / ry) ** 2)
    ring = np.abs(r - 1.0) * min(rx, ry) < t / 2.0
    paint(img[ya:yb, xa:xb], ring, rgb)


def text_mask(s, size):
    t = pyvips.Image.text(s, font="sans bold %d" % size, dpi=72)
    return np.ndarray(buffer=t.write_to_memory(), dtype=np.uint8,
                      shape=(t.height, t.width)).astype(np.float32) / 255.0


def put_text(img, s, x, y, rgb, size, box=True):
    """Draw text with its top-left at (x, y); clipped to the image. Returns width."""
    m = text_mask(s, size)
    h, w = img.shape[:2]
    pad = 6
    if box:
        bx0, by0 = max(x - pad, 0), max(y - pad, 0)
        bx1, by1 = min(x + m.shape[1] + pad, w), min(y + m.shape[0] + pad, h)
        if bx0 < bx1 and by0 < by1:
            img[by0:by1, bx0:bx1] *= 0.35
    x0, y0 = max(x, 0), max(y, 0)
    x1, y1 = min(x + m.shape[1], w), min(y + m.shape[0], h)
    if x0 < x1 and y0 < y1:
        paint(img[y0:y1, x0:x1], m[y0 - y:y1 - y, x0 - x:x1 - x], rgb)
    return m.shape[1]


def caption_strip(width, lines):
    """lines: list of [(text, rgb), ...] runs. Black strip under the image.
    Runs wrap word by word, since a hotspot PNG is only ~1660 px wide."""
    size, lh, margin, gap = 26, 38, 16, 14
    placed, row = [], 0
    for runs in lines:
        x = margin
        for s, rgb in runs:
            for word in s.split(" "):
                if not word:
                    continue
                wpx = text_mask(word, size).shape[1]
                if x > margin and x + wpx > width - margin:
                    row, x = row + 1, margin + 2 * gap
                placed.append((word, rgb, x, row))
                x += wpx + 9
            x += gap
        row += 1
    strip = np.zeros((margin * 2 + lh * row, width, 3), np.float32)
    for word, rgb, x, r in placed:
        put_text(strip, word, x, margin + r * lh, rgb, size, box=False)
    return strip


def main(argv):
    sample = next((a for a in argv if not a.startswith("--")), "IB1")
    smoke = "--smoke" in argv
    omit = () if "--with-G" in argv else ("G",)
    crops_dir = os.path.join(OUT_ROOT, "crops")
    out_dir = os.path.join(OUT_ROOT, "roi_png_smoke" if smoke else "roi_png")
    os.makedirs(out_dir, exist_ok=True)

    regions = [r for r in read_csv(os.path.join(OUT_ROOT, "%s_regions_canvas_px.csv" % sample))
               if r["group"] in GROUPS]
    for r in regions:
        for k in ("x1", "y1", "x2", "y2"):
            r[k] = float(r[k])
    blank = {r["id"]: float(r["blank_frac"])
             for r in read_csv(os.path.join(OUT_ROOT, "%s_regions_L5.csv" % sample))}
    crops = read_csv(os.path.join(crops_dir, "%s_crop_manifest.csv" % sample))
    for c in crops:
        for k in ("x0", "y0", "width", "height"):
            c[k] = int(c[k])
    assert len(regions) == 45 or sample != "IB1", "expected 45 IB1 regions, got %d" % len(regions)

    ranges = crop_ranges(os.path.join(crops_dir, crops[0]["crop_file"]))
    log("shared display ranges: " + ", ".join(
        "%s %.0f-%.0f" % (n, lo, hi) for (n, _), (lo, hi) in zip(CHANNELS, ranges)))
    log("omitted channels: %s" % (", ".join(omit) or "none"))

    pad = int(round(PAD_UM / UM_PER_PX))
    todo = regions[:2] if smoke else regions
    index = []
    for n, r in enumerate(todo, 1):
        # Padded window in canvas px, then the crop that holds it with the most room.
        wx0, wy0 = int(np.floor(r["x1"])) - pad, int(np.floor(r["y1"])) - pad
        wx1, wy1 = int(np.ceil(r["x2"])) + pad, int(np.ceil(r["y2"])) + pad

        def margin(c):
            return min(wx0 - c["x0"], wy0 - c["y0"],
                       c["x0"] + c["width"] - wx1, c["y0"] + c["height"] - wy1)
        c = max(crops, key=margin)
        clipped = margin(c) < 0
        wx0, wy0 = max(wx0, c["x0"]), max(wy0, c["y0"])
        wx1, wy1 = min(wx1, c["x0"] + c["width"]), min(wy1, c["y0"] + c["height"])

        tif = os.path.join(crops_dir, c["crop_file"])
        assert crop_ranges(tif) == ranges, "%s has different display ranges" % tif
        mm = tifffile.memmap(tif, mode="r")
        assert mm.shape == (len(CHANNELS), c["height"], c["width"]), mm.shape
        arr = np.asarray(mm[:, wy0 - c["y0"]:wy1 - c["y0"], wx0 - c["x0"]:wx1 - c["x0"]])
        del mm
        img = composite(arr, ranges, omit)
        del arr
        np.clip(img, 0, 255, out=img)

        # Every exported region overlapping the window, the subject drawn last.
        inside = [o for o in regions
                  if o["x2"] > wx0 and o["x1"] < wx1 and o["y2"] > wy0 and o["y1"] < wy1]
        inside.sort(key=lambda o: o["id"] == r["id"])
        for o in inside:
            box = (o["x1"] - wx0, o["y1"] - wy0, o["x2"] - wx0, o["y2"] - wy0)
            draw = draw_rect if o["shape"] == "rect" else draw_oval
            draw(img, *box, GROUP_COLOR[o["group"]], LINE_PX)
            # Hotspot bins tile and share edges, so a label outside its box sits on
            # the neighbour's edge and reads as the neighbour's; put it inside the
            # visible top-left corner instead. Ovals keep theirs outside the curve.
            if o["shape"] == "rect":
                lx, ly = int(max(box[0], 0)) + LINE_PX + 10, int(max(box[1], 0)) + LINE_PX + 10
            else:
                lx = int(max(box[0], 8))
                ly = int(box[1]) - 44 if box[1] >= 50 else int(box[3]) + 12
            put_text(img, short_label(o), lx, ly,
                     GROUP_COLOR[o["group"]], 30 if o["id"] == r["id"] else 24)

        # Scale bar, bottom right.
        h, w = img.shape[:2]
        bl = int(round(SCALEBAR_UM / UM_PER_PX))
        bx1, by = w - 40, h - 50
        img[by - 6:by + 6, bx1 - bl - 20:bx1 + 20] *= 0.35
        img[by - 4:by + 4, bx1 - bl:bx1] = 255
        put_text(img, "%d µm" % SCALEBAR_UM, bx1 - bl, by - 52, (255, 255, 255), 26)

        is_hot = r["group"] == "hodgkin_hotspot"
        head = "%s   %s" % (r["id"], GROUP_NAME[r["group"]])
        head += ("   CCL17+CCL22 %.1f cp10k, 100 µm bin" % float(r["metric"]) if is_hot
                 else "   kappa fraction %.2f, n = %s cells" % (float(r["metric"]), r["n_cells"]))
        warn = []
        if r["flag"] == "seam":
            warn.append("SEAM: within 25 µm of an FOV edge")
        if blank.get(r["id"], 0) > 0.01:
            warn.append("%.0f%% of the region has no morphology signal (FOV blank in the image, "
                        "transcripts still real)" % (100 * blank[r["id"]]))
        if clipped:
            warn.append("context clipped at the crop edge")
        lines = [[(head, GROUP_COLOR[r["group"]])]]
        if warn:
            lines.append([("  |  ".join(warn), (255, 200, 60))])
        lines.append([("Channels:", (200, 200, 200))]
                     + [("%s %.0f-%.0f%s" % (nm, lo, hi,
                                             " x%g" % GAIN[nm] if nm in GAIN else ""), rgb)
                        for (nm, rgb), (lo, hi) in zip(CHANNELS, ranges) if nm not in omit]
                     + [("(fixed ranges and gains, same in every PNG)", (200, 200, 200))]
                     + [("%s omitted (marker unknown)" % "/".join(omit), (200, 200, 200))
                        for _ in [0] if omit])
        lines.append([("Outlines:", (200, 200, 200)),
                       ("cyan box = hotspot 100 µm bin (c = CCL17+CCL22 cp10k)",
                        GROUP_COLOR["hodgkin_hotspot"]),
                       ("orange oval = FL neighbourhood, radius = span (k = kappa fraction)",
                        GROUP_COLOR["fl_nb"])])
        lines.append([("Full resolution %.4f µm/px; canvas px origin (%d, %d); from crops/%s. "
                       "Hotspots are defined by chemokine expression, not by identified HRS cells."
                       % (UM_PER_PX, wx0, wy0, c["crop_file"]), (200, 200, 200))])
        cap = caption_strip(w, lines)
        img = np.concatenate([img, cap], axis=0)

        out_f = os.path.join(out_dir, "%s_%s.png" % (sample, r["id"].replace("%s_" % sample, "", 1))
                             if not r["id"].startswith(sample) else "%s.png" % r["id"])
        rgb8 = np.ascontiguousarray(img.round().astype(np.uint8))
        pyvips.Image.new_from_memory(rgb8.tobytes(), rgb8.shape[1], rgb8.shape[0], 3,
                                     "uchar").pngsave(out_f, compression=6)
        index.append({"id": r["id"], "group": r["group"], "flag": r["flag"],
                      "metric": r["metric"], "blank_frac": blank.get(r["id"], ""),
                      "png": os.path.basename(out_f), "source_crop": c["crop_file"],
                      "canvas_x0": wx0, "canvas_y0": wy0, "width_px": w, "height_px": h,
                      "context_clipped": clipped,
                      "other_regions": " ".join(o["id"] for o in inside if o["id"] != r["id"])})
        log("[%2d/%d] %s  %dx%d px  from %s%s"
            % (n, len(todo), os.path.basename(out_f), w, h, c["crop_file"],
               "  CLIPPED" if clipped else ""))

    with open(os.path.join(out_dir, "roi_png_index.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(index[0]))
        wr.writeheader()
        wr.writerows(index)
    log("done: %d PNGs in %s" % (len(index), out_dir))


if __name__ == "__main__":
    main(sys.argv[1:])
