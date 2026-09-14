"""STAGE 2 of the ImageJ/Fiji export -- pixels and ROIs.

Reads output/15_imagej/<SAMPLE>_regions_canvas_px.csv (stage 1,
scripts/15_script_imagej_export.R) and writes, under output/15_imagej/:

  <SAMPLE>_overview_L<N>.tif        whole slide at pyramid level N, 5 channels
  <SAMPLE>_overview_L<N>_RoiSet.zip every region, in level-N coordinates
  <SAMPLE>_overview_L<N>.ijm        the same ROIs as a macro
  <SAMPLE>_regions_L<N>.csv         the scaled coordinates, for inspection
  crops/<SAMPLE>_w##_<seed>.tif     full-resolution crop windows
  crops/<SAMPLE>_w##_<seed>_RoiSet.zip / .ijm   crop-local ROIs
  crops/<SAMPLE>_crop_manifest.csv  which regions landed in which crop

WHY NO RE-PYRAMIDING. IB1.ome.tiff is already a 10-level power-of-2 pyramid
(512x512 tiles, Deflate+predictor BigTIFF). A downsampled overview is a direct
read of an existing level, and a full-resolution crop is a lazy read of ~400
tiles via the zarr view -- neither ever materialises the 213 GB level 0.

WHY ImageJ-FORMAT TIFF RATHER THAN OME. `imagej=True` writes a native Fiji
hyperstack: it opens instantly with channel sliders and a correct scale bar, with
no Bio-Formats round-trip and no chance of the importer handing back the wrong
pyramid level. The cost is no pyramid inside the output, which is fine because
every output here is already sized to open whole.

CONTRAST. Ranges are the 30th/99th percentile of NON-ZERO pixels, per channel --
the same intent as .cosmx_gquants_to_precalc(min_q=0.30, max_q=0.99) in
scripts/lib/cosmx_tiff.R:163, whose comment records that the 30th percentile is 0
in every group and so resolves to "smallest non-zero quantile". That function
reads a cached quantile grid; this computes the two percentiles on pixels
directly, which is equivalent in intent and needs no cache.
Ranges are computed PER RESOLUTION (overview from its own level, crops from a
sample of level-0 windows) because downsampling averages away the peaks -- the
reason the project's precalc table is keyed by resolution level at all.

Usage:
  python scripts/15_imagej_tiles.py [SAMPLE] [LEVEL] [CROP_PX] [GROUPS]
    SAMPLE   IB1 | RB1                              default IB1
    LEVEL    overview pyramid level                  default 5
    CROP_PX  full-res crop window side, in px        default 4256 (one FOV)
    GROUPS   comma-separated groups to export        default hodgkin_hotspot,fl_nb
  --no-crops        overview only
  --no-overview     crops only
"""

import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))

import numpy as np
import tifffile
import zarr

import imagej_roi as ij

# ---------------------------------------------------------------- config ----
TIFF_DIR = "/gpfs2/scratch/pi/avolaric/P3_CosMx_assemble_FOVs"
OUT_DIR = "output/15_imagej"
MM_PER_PX = 120.280945e-6
UM_PER_PX_L0 = MM_PER_PX * 1000.0

# Positional, matching cosmx_channel_names() / cosmx_channel_colors()
# (scripts/lib/cosmx_tiff.R:136-144). Reordering silently recolours everything.
CHANNELS = [("PanCK", (0, 0, 255)),
            ("G", (0, 255, 0)),
            ("Membrane", (255, 255, 0)),
            ("CD45", (255, 0, 0)),
            ("DNA", (238, 130, 238))]

# Hodgkin blue / FL orange are the project's disease hues
# (12_script_fl_vs_hodgkin_overlay.R). #0b3d7a is invisible on a dark
# fluorescence image, so the ROI outline uses a bright cyan in its place.
GROUP_COLOR = {"hodgkin_hotspot": "#00e5ff",
               "fl_nb": "#eb6834",
               "fl_focus": "#ffcc00"}

LO_PCT, HI_PCT = 30.0, 99.0
BLANK_LEVEL = 6            # coarse level used only for the blank-morphology check
RANGE_SAMPLE_WINDOWS = 12   # level-0 windows sampled to set the crop contrast
RANGE_SAMPLE_STRIDE = 4


def log(msg):
    print(msg, flush=True)


def parse_args(argv):
    pos, flags = [], set()
    for a in argv:
        (flags.add(a) if a.startswith("--") else pos.append(a))
    return (pos[0] if len(pos) > 0 else "IB1",
            int(pos[1]) if len(pos) > 1 else 5,
            int(pos[2]) if len(pos) > 2 else 4256,
            (pos[3] if len(pos) > 3 else "hodgkin_hotspot,fl_nb").split(","),
            flags)


def luts():
    """Per-channel 256-entry ramps, as ImageJ wants them: (3, 256) uint8."""
    out = []
    ramp = np.arange(256, dtype=np.float32) / 255.0
    for _, (r, g, b) in CHANNELS:
        out.append(np.stack([(ramp * c).astype(np.uint8) for c in (r, g, b)]))
    return out


def percentile_ranges(sample_per_channel):
    """(lo, hi) per channel from non-zero pixels; falls back when a channel is
    entirely zero or flat, since ImageJ rejects a zero-width display range."""
    ranges = []
    for ci, vals in enumerate(sample_per_channel):
        nz = vals[vals > 0]
        if nz.size == 0:
            lo, hi = 0.0, 1.0
        else:
            lo = float(np.percentile(nz, LO_PCT))
            hi = float(np.percentile(nz, HI_PCT))
            if hi <= lo:
                hi = lo + 1.0
        ranges.append((lo, hi))
        log("    %-9s range %8.1f - %8.1f  (non-zero %.1f%% of %d px)"
            % (CHANNELS[ci][0], lo, hi,
               100.0 * (vals > 0).mean() if vals.size else 0.0, vals.size))
    return ranges


def write_hyperstack(path, arr, um_per_px, ranges):
    """arr is (C, Y, X) uint16 -> a Fiji composite hyperstack."""
    flat = tuple(float(v) for lohi in ranges for v in lohi)
    tifffile.imwrite(
        path, arr,
        imagej=True,
        resolution=(1.0 / um_per_px, 1.0 / um_per_px),
        metadata={"axes": "CYX", "unit": "um", "mode": "composite",
                  "Ranges": flat, "LUTs": luts(),
                  "Labels": [n for n, _ in CHANNELS]},
    )
    return os.path.getsize(path)


def blank_fractions(grp, regions, level, n_levels):
    """Fraction of each region's box where EVERY channel is zero.

    Some FOVs in the assembled canvas carry no morphology signal at all -- they
    show up as clean black squares. A bin can still be a real expression call
    there (the transcript data is independent of the morphology channels), so
    these are not errors, but a crop that is half black is worth knowing about
    before you open 34 of them. A whole missing FOV is 4256 px, which is 66 px at
    level 6, so a coarse level resolves this fine and costs ~50 MB.
    """
    lv = min(level, n_levels - 1)
    z = grp[str(lv)]
    scale = 1.0 / (2 ** lv)
    arr = np.asarray(z[:])
    blank = (arr == 0).all(axis=0)
    del arr
    h, w = blank.shape
    out = {}
    for r in regions:
        x0 = max(0, int(r["x1"] * scale))
        y0 = max(0, int(r["y1"] * scale))
        x1 = min(w, int(np.ceil(r["x2"] * scale)))
        y1 = min(h, int(np.ceil(r["y2"] * scale)))
        sub = blank[y0:y1, x0:x1]
        out[r["id"]] = float(sub.mean()) if sub.size else 1.0
    return out


def read_regions(sample, groups):
    f = os.path.join(OUT_DIR, "%s_regions_canvas_px.csv" % sample)
    if not os.path.exists(f):
        sys.exit("missing %s -- run scripts/15_script_imagej_export.R first" % f)
    keep, skipped = [], 0
    with open(f) as fh:
        for row in csv.DictReader(fh):
            if row["group"] not in groups:
                skipped += 1
                continue
            for k in ("x1", "y1", "x2", "y2", "cx", "cy", "w", "h", "metric"):
                row[k] = float(row[k])
            keep.append(row)
    log("read %d regions from %s (%d rows in other groups, skipped)"
        % (len(keep), f, skipped))
    if not keep:
        sys.exit("no regions in groups %s" % groups)
    return keep


def roi_dicts(regions, scale=1.0, ox=0.0, oy=0.0, label_metric=True):
    """Regions -> ROI dicts in the target image's coordinate space.

    scale is 1/2**level for the overview, 1.0 for a full-resolution crop; ox/oy
    is the crop window origin in canvas px.
    """
    out = []
    for r in regions:
        x = (r["x1"] - ox) * scale
        y = (r["y1"] - oy) * scale
        w = max((r["x2"] - r["x1"]) * scale, 1.0)
        h = max((r["y2"] - r["y1"]) * scale, 1.0)
        name = r["id"]
        if label_metric:
            name = "%s_%s%g" % (r["id"], "k" if r["group"].startswith("fl") else "c",
                                round(r["metric"], 2))
        if r["flag"]:
            name += "_" + r["flag"]
        out.append(dict(x=x, y=y, width=w, height=h,
                        shape=r["shape"], name=name,
                        color=GROUP_COLOR.get(r["group"]), stroke_width=2))
    return out


def plan_windows(regions, crop_px, canvas_w, canvas_h):
    """Greedy: seed a window on the top-ranked unassigned region, absorb every
    other unassigned region that fits entirely inside it, repeat.

    Guarantees each region sits wholly within exactly one window, so a crop never
    shows half a hotspot. A region wider than crop_px (an FL neighbourhood with a
    large span_mm) grows its own window rather than being clipped.
    """
    order = sorted(range(len(regions)),
                   key=lambda i: (regions[i]["group"], -regions[i]["metric"]))
    unassigned = list(order)
    windows = []
    while unassigned:
        seed = unassigned[0]
        s = regions[seed]
        side = max(crop_px, int(np.ceil(s["w"])) + 16, int(np.ceil(s["h"])) + 16)
        x0 = int(round(s["cx"] - side / 2.0))
        y0 = int(round(s["cy"] - side / 2.0))
        # clip to canvas, keeping the requested size where there is room
        x0 = max(0, min(x0, canvas_w - side)) if side <= canvas_w else 0
        y0 = max(0, min(y0, canvas_h - side)) if side <= canvas_h else 0
        x1 = min(canvas_w, x0 + side)
        y1 = min(canvas_h, y0 + side)

        members, rest = [], []
        for i in unassigned:
            r = regions[i]
            if r["x1"] >= x0 and r["x2"] <= x1 and r["y1"] >= y0 and r["y2"] <= y1:
                members.append(i)
            else:
                rest.append(i)
        if seed not in members:      # can only happen if the seed exceeds the canvas
            members.append(seed)
            rest = [i for i in rest if i != seed]
        windows.append(dict(x0=x0, y0=y0, x1=x1, y1=y1,
                            seed=regions[seed]["id"],
                            group=regions[seed]["group"],
                            members=members))
        unassigned = rest
    return windows


def main():
    sample, level, crop_px, groups, flags = parse_args(sys.argv[1:])
    tiff_path = os.path.join(TIFF_DIR, "%s.ome.tiff" % sample)
    if not os.path.exists(tiff_path):
        sys.exit("no TIFF at %s (scratch may have been cleaned)" % tiff_path)
    crops_dir = os.path.join(OUT_DIR, "crops")
    os.makedirs(crops_dir, exist_ok=True)

    log("=" * 60)
    log("step 15 stage 2 -- %s, overview L%d, crop %d px" % (sample, level, crop_px))
    log("groups: %s" % ", ".join(groups))
    log("=" * 60)

    regions = read_regions(sample, groups)
    for g in groups:
        n = sum(1 for r in regions if r["group"] == g)
        log("  %-16s %3d regions" % (g, n))

    with tifffile.TiffFile(tiff_path) as tf:
        series = tf.series[0]
        store = series.aszarr()
        try:
            grp = zarr.open(store, mode="r")
            n_levels = len(series.levels)
            if not 0 <= level < n_levels:
                sys.exit("level %d out of range (0-%d)" % (level, n_levels - 1))
            z0 = grp["0"]
            _, canvas_h, canvas_w = z0.shape
            log("\ncanvas %d x %d px, %d channels, %d pyramid levels"
                % (canvas_w, canvas_h, z0.shape[0], n_levels))

            log("\n-- blank-morphology check (level %d)" % BLANK_LEVEL)
            blank = blank_fractions(grp, regions, BLANK_LEVEL, n_levels)
            flagged = sorted((v, k) for k, v in blank.items() if v > 0.01)
            if flagged:
                log("   %d of %d regions overlap FOVs with no morphology signal."
                    % (len(flagged), len(regions)))
                log("   These are still real expression calls -- the transcript data "
                    "does not depend on\n   the morphology channels -- but their crops "
                    "will be partly black:")
                for v, k in reversed(flagged):
                    log("     %-18s %5.1f%% blank" % (k, v * 100))
            else:
                log("   all %d regions sit on imaged tissue" % len(regions))

            # -------------------------------------------------- overview ----
            if "--no-overview" not in flags:
                zl = grp[str(level)]
                scale = 1.0 / (2 ** level)
                um = UM_PER_PX_L0 * (2 ** level)
                log("\n-- overview: L%d %s, %.1f MB, %.3f um/px"
                    % (level, zl.shape, zl.nbytes / 1e6, um))
                ov = np.asarray(zl[:])
                log("  computing display ranges")
                ranges = percentile_ranges([ov[c].ravel() for c in range(ov.shape[0])])
                ov_f = os.path.join(OUT_DIR, "%s_overview_L%d.tif" % (sample, level))
                sz = write_hyperstack(ov_f, ov, um, ranges)
                log("  wrote %s (%.1f MB)" % (ov_f, sz / 1e6))
                del ov

                rois = roi_dicts(regions, scale=scale)
                too_big = [r for r in rois
                           if max(r["x"] + r["width"], r["y"] + r["height"]) > ij.INT16_MAX]
                if too_big:
                    log("  WARNING: %d ROI(s) exceed the int16 limit at L%d -- "
                        "macro only. Use a coarser level for binary ROIs."
                        % (len(too_big), level))
                else:
                    zf = os.path.join(OUT_DIR,
                                      "%s_overview_L%d_RoiSet.zip" % (sample, level))
                    log("  wrote %s (%d ROIs)" % (zf, ij.write_roi_zip(zf, rois)))
                mf = os.path.join(OUT_DIR, "%s_overview_L%d.ijm" % (sample, level))
                ij.write_macro(mf, rois,
                               title="%s overview, pyramid level %d (%.3f um/px)"
                                     % (sample, level, um),
                               note="Coordinates are level-%d pixels.\n"
                                    "Open %s first." % (level, os.path.basename(ov_f)))
                log("  wrote %s" % mf)

                cf = os.path.join(OUT_DIR, "%s_regions_L%d.csv" % (sample, level))
                with open(cf, "w", newline="") as fh:
                    wr = csv.writer(fh)
                    wr.writerow(["group", "id", "shape", "x", "y", "width", "height",
                                 "metric", "metric_name", "n_cells", "flag",
                                 "blank_frac", "note"])
                    for r, d in zip(regions, rois):
                        wr.writerow([r["group"], r["id"], r["shape"],
                                     round(d["x"], 2), round(d["y"], 2),
                                     round(d["width"], 2), round(d["height"], 2),
                                     r["metric"], r["metric_name"], r["n_cells"],
                                     r["flag"], round(blank[r["id"]], 4), r["note"]])
                log("  wrote %s" % cf)

            # ----------------------------------------------------- crops ----
            if "--no-crops" in flags:
                log("\n-- crops skipped (--no-crops)")
            else:
                windows = plan_windows(regions, crop_px, canvas_w, canvas_h)
                total = sum((w["x1"] - w["x0"]) * (w["y1"] - w["y0"])
                            * z0.shape[0] * 2 for w in windows)
                log("\n-- crops: %d regions -> %d windows, ~%.1f GB total"
                    % (len(regions), len(windows), total / 1e9))
                sides = sorted({(w["x1"] - w["x0"]) for w in windows})
                log("   window widths: %s px" % ", ".join(str(s) for s in sides))

                # One contrast for every crop, so they are comparable. Sampled
                # from level-0 windows spread across the list.
                log("  sampling level-0 windows for a shared display range")
                pick = np.linspace(0, len(windows) - 1,
                                   min(RANGE_SAMPLE_WINDOWS, len(windows))
                                   ).round().astype(int)
                acc = [[] for _ in range(z0.shape[0])]
                for i in sorted(set(pick.tolist())):
                    w = windows[i]
                    sub = np.asarray(z0[:, w["y0"]:w["y1"]:RANGE_SAMPLE_STRIDE,
                                        w["x0"]:w["x1"]:RANGE_SAMPLE_STRIDE])
                    for c in range(sub.shape[0]):
                        acc[c].append(sub[c].ravel())
                    del sub
                ranges = percentile_ranges([np.concatenate(a) for a in acc])
                del acc

                man_f = os.path.join(crops_dir, "%s_crop_manifest.csv" % sample)
                with open(man_f, "w", newline="") as fh:
                    wr = csv.writer(fh)
                    wr.writerow(["crop_file", "window", "seed_group", "seed_id",
                                 "x0", "y0", "width", "height", "n_regions",
                                 "max_blank_frac", "region_ids"])
                    for wi, w in enumerate(windows):
                        tag = "%s_w%02d_%s" % (sample, wi, w["seed"])
                        cw, ch = w["x1"] - w["x0"], w["y1"] - w["y0"]
                        arr = np.asarray(z0[:, w["y0"]:w["y1"], w["x0"]:w["x1"]])
                        tif_f = os.path.join(crops_dir, tag + ".tif")
                        sz = write_hyperstack(tif_f, arr, UM_PER_PX_L0, ranges)
                        del arr

                        mem = [regions[i] for i in w["members"]]
                        rois = roi_dicts(mem, scale=1.0, ox=w["x0"], oy=w["y0"])
                        zf = os.path.join(crops_dir, tag + "_RoiSet.zip")
                        ij.write_roi_zip(zf, rois)
                        ij.write_macro(
                            os.path.join(crops_dir, tag + ".ijm"), rois,
                            title="%s crop %d, full resolution (%.4f um/px)"
                                  % (sample, wi, UM_PER_PX_L0),
                            note="Canvas-px window origin (%d, %d), %dx%d.\n"
                                 "Coordinates are crop-local. Open %s first."
                                 % (w["x0"], w["y0"], cw, ch, tag + ".tif"))
                        mb = max(blank[m["id"]] for m in mem)
                        wr.writerow([tag + ".tif", wi, w["group"], w["seed"],
                                     w["x0"], w["y0"], cw, ch, len(mem),
                                     round(mb, 4),
                                     " ".join(m["id"] for m in mem)])
                        log("  [%2d/%d] %s  %dx%d  %d region(s)  %.0f MB%s"
                            % (wi + 1, len(windows), tag, cw, ch, len(mem), sz / 1e6,
                               "  %.0f%% BLANK" % (mb * 100) if mb > 0.01 else ""))
                log("  wrote %s" % man_f)
        finally:
            store.close()

    with open(os.path.join(OUT_DIR, "%s_tiles.done" % sample), "w") as fh:
        fh.write("ok\n")
    log("\ndone.")


if __name__ == "__main__":
    main()
