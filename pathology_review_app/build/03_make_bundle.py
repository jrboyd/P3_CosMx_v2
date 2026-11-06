"""Pathology review bundle, stage 3 of 3: blind, order, repeat, zip.

    python pathology_review_app/build/03_make_bundle.py test     # stdlib only, light

Reads output/<size>/{rendered.csv, ranges.json, tiles/} and build/review_form.json; writes

  output/<size>/send/<bundle_id>.zip      the data blob the reader opens in the app
  output/<size>/send/reviewer.html        copy of app/reviewer.html, to send alongside
  output/<size>/send/READ_ME_FIRST.txt
  keys/<bundle_id>_key.csv                item_id -> cell, kept HERE, never in the zip

BLINDING. Items get random 4-character ids; nothing in the zip names a cell, sample,
stratum, count or score. Repeats get their own fresh id and their own copy of the files,
so a repeat cannot be linked to its original from the file listing.

ORDER. Stratified interleave, seeded: within each stratum the cells are shuffled and
spread evenly across the sequence ((rank + U) / n), so every stretch of the review mixes
all strata and reading order is independent of brightness, score, sample and kind.

REPEATS (intra-reader agreement). A random draw of originals from the first 60% of the
sequence, each re-inserted at least max(10, n/3) positions after its original.

The zip is STORED (no compression; JPEGs do not compress) so the app can read each file
by offset without a zip library and without loading the whole blob into memory.
"""
import csv
import datetime
import json
import os
import random
import shutil
import sys
import zipfile

SIZE = sys.argv[1] if len(sys.argv) > 1 else "test"
ROOT = "pathology_review_app"
OUT = os.path.join(ROOT, "output", SIZE)
N_REPEAT = {"test": 5, "full": 20}[SIZE]
SEED = 38
CH_DISPLAY = [  # positional colours of the project composite (CLAUDE.md); G omitted
    {"name": "DNA",      "rgb": [238, 130, 238], "gain": 1.0},
    {"name": "CD45",     "rgb": [255, 0, 0],     "gain": 1.0},
    {"name": "Membrane", "rgb": [255, 255, 0],   "gain": 0.5},
    {"name": "PanCK",    "rgb": [0, 0, 255],     "gain": 1.0},
]
ID_ALPHABET = "ACDEFGHJKLMNPQRTUVWXY34679"     # no 0/O, 1/I/L-lookalikes, 5/S, 8/B, 2/Z


def main():
    rnd = random.Random(SEED)
    with open(os.path.join(OUT, "rendered.csv"), newline="") as fh:
        cells = list(csv.DictReader(fh))
    with open(os.path.join(OUT, "ranges.json")) as fh:
        ranges = json.load(fh)
    with open(os.path.join(ROOT, "build", "review_form.json")) as fh:
        form = json.load(fh)

    # ---- order: stratified interleave ---------------------------------------------
    by = {}
    for c in cells:
        by.setdefault(c["stratum"], []).append(c)
    keyed = []
    for stratum in sorted(by):
        grp = by[stratum][:]
        rnd.shuffle(grp)
        n = len(grp)
        for i, c in enumerate(grp):
            keyed.append(((i + rnd.random()) / n, c))
    keyed.sort(key=lambda t: t[0])
    seq = [dict(c, repeat_of="") for _, c in keyed]

    # ---- ids -----------------------------------------------------------------------
    used = set()

    def new_id():
        while True:
            s = "".join(rnd.choice(ID_ALPHABET) for _ in range(4))
            if s not in used:
                used.add(s)
                return s

    for c in seq:
        c["item_id"] = new_id()

    # ---- repeats -------------------------------------------------------------------
    n0 = len(seq)
    gap = max(10, n0 // 3)
    eligible = [c for c in seq[: int(n0 * 0.6)]]
    originals = rnd.sample(eligible, min(N_REPEAT, len(eligible)))
    for o in originals:
        pos0 = seq.index(o)
        lo = min(pos0 + gap, len(seq))
        pos = rnd.randint(lo, len(seq))
        seq.insert(pos, dict(o, item_id=new_id(), repeat_of=o["item_id"]))
    for o in originals:  # insertions can shift things; verify the spacing actually held
        rep = next(c for c in seq if c["repeat_of"] == o["item_id"])
        assert seq.index(rep) - seq.index(o) >= gap - 1, "repeat too close to its original"

    # ---- zip + manifest ------------------------------------------------------------
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M")
    bundle_id = "P3HRS-%s%d-%s" % (SIZE, len(seq), stamp)
    send = os.path.join(OUT, "send")
    if os.path.isdir(send):
        shutil.rmtree(send)
    os.makedirs(send)
    items = []
    for i, c in enumerate(seq, 1):
        items.append({
            "id": c["item_id"], "seq": i,
            "cx": float(c["cx_win"]), "cy": float(c["cy_win"]),
            "files": {ch["name"]: "images/%s_%s.jpg" % (c["item_id"], ch["name"])
                      for ch in CH_DISPLAY},
        })
    manifest = {
        "bundle_id": bundle_id,
        "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "project": "P3_CosMx_v2 HRS morphology review",
        "n_items": len(items),
        "um_per_px": ranges["um_per_px"],
        "window_px": ranges["window_px"],
        "channels": CH_DISPLAY,
        "form": form,
        "items": items,
    }
    zpath = os.path.join(send, bundle_id + ".zip")
    with zipfile.ZipFile(zpath, "w", compression=zipfile.ZIP_STORED, allowZip64=False) as z:
        z.writestr("manifest.json", json.dumps(manifest, indent=1))
        for c in seq:
            for ch in CH_DISPLAY:
                z.write(os.path.join(OUT, "tiles", c["cell_id"], ch["name"] + ".jpg"),
                        "images/%s_%s.jpg" % (c["item_id"], ch["name"]))
    shutil.copy(os.path.join(ROOT, "app", "reviewer.html"), os.path.join(send, "reviewer.html"))
    shutil.copy(os.path.join(ROOT, "app", "READ_ME_FIRST.txt"), os.path.join(send, "READ_ME_FIRST.txt"))

    # ---- key (stays in the repo) ---------------------------------------------------
    os.makedirs(os.path.join(ROOT, "keys"), exist_ok=True)
    kpath = os.path.join(ROOT, "keys", bundle_id + "_key.csv")
    kcols = ["item_id", "seq", "repeat_of", "cell_id", "samp", "kind", "stratum", "src_group",
             "chemo", "hrs_score", "k3", "Area.um2", "depth0", "matched_to", "fov_id",
             "X_px", "Y_px", "x0", "y0", "cx_win", "cy_win", "blank_frac"]
    with open(kpath, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=kcols, extrasaction="ignore")
        wr.writeheader()
        for i, c in enumerate(seq, 1):
            wr.writerow(dict(c, seq=i))

    mb = os.path.getsize(zpath) / 1e6
    print("bundle %s: %d items (%d cells + %d repeats), %.1f MB -> %s"
          % (bundle_id, len(seq), n0, len(originals), mb, zpath))
    print("key -> %s" % kpath)
    strata = {}
    for c in seq:
        strata[c["stratum"]] = strata.get(c["stratum"], 0) + 1
    for k in sorted(strata):
        print("  %-22s %d" % (k, strata[k]))


if __name__ == "__main__":
    main()
