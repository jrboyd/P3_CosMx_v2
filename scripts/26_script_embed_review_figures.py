#!/usr/bin/env python3
"""
step 26 -- embed the figure-list PNGs into reports/manuscript_plan.html as data: URIs.

WHY A SCRIPT AND NOT A ONE-OFF EDIT. The review is a living document: panels get rebuilt,
statuses change, and the whole thing has been regenerated twice already. Hand-pasting 30
base64 blobs would be unmaintainable and would rot the moment a panel is redrawn. This is
re-runnable -- it strips any previously embedded strip before inserting the current one, so
running it after rebuilding panels just refreshes them.

WHY data: URIs. The house rule is that a report is a SELF-CONTAINED local HTML file with no
external references, so it renders correctly opened straight off GPFS with no network. An
<img src="../output/..."> would break the moment the file is copied off the cluster, which is
exactly how these get shared.

SIZE DISCIPLINE. Sources run to 5 MB apiece and there are ~30 of them; embedded raw that is a
150 MB page. Each is downscaled to a review-appropriate width and re-encoded:
  - plots (line art)      -> PNG, which keeps text crisp
  - morphology snapshots  -> JPEG, because a photograph in PNG is enormous and JPEG artefacts
                             are invisible at thumbnail scale
The full-resolution path stays printed under every thumbnail, so the original is always
findable; the embedded copy is a preview, not a replacement.

PATH RESOLUTION. The review cites paths in four forms -- full (`output/25_figures/x.png`),
project-relative (`19_celltypes/figures/x.png`), bare (`fig_11_niche_map_IB1.png`), and with
shell-style expansion (`fig_22_snap_CXCL13CXCR5_{cherry,repr}.png`, `tfh_*.png`). All four are
resolved; anything that cannot be found is reported and skipped rather than silently dropped.

Usage: python3 scripts/26_script_embed_review_figures.py
"""
import base64, glob, os, re, subprocess, sys, tempfile

ROOT = "/gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2"
HTML = os.path.join(ROOT, "reports/manuscript_plan.html")
W_MAIN, W_SUPP = 1000, 660
MAX_PER_REF, MAX_PER_BLOCK = 2, 8
## Directories holding NON-PRIMARY variants of a figure. A bare filename like
## `fig_11_niche_map_IB1.png` exists in both 11_hodgkin_niche/ and 11_hodgkin_niche_preqc/,
## and `fig_16_regions_*.png` exists at the matched, strict and preqc tiers. Embedding a preqc
## panel unlabelled beside the primary one would put a map CLAUDE.md says is for POSITION ONLY
## next to one used for magnitude -- so they are excluded, and the matched tier wins over
## strict wherever both exist.
EXCLUDE = ("preqc", "_prev", "prev_", "SUPERSEDED", "noplasmacut")
JPEG_Q = 86
MARK_A, MARK_B = "<!--EMBED:START-->", "<!--EMBED:END-->"

def expand(ref):
    """brace + glob expansion, then resolve against the project."""
    refs = [ref]
    m = re.search(r"\{([^}]*)\}", ref)
    if m:
        refs = [ref[:m.start()] + o + ref[m.end():] for o in m.group(1).split(",")]
    out = []
    for r in refs:
        cands = [os.path.join(ROOT, r), os.path.join(ROOT, "output", r)]
        hits = []
        for c in cands:
            hits += sorted(glob.glob(c))
        if not hits:                       # bare filename: search under output/
            hits = sorted(glob.glob(os.path.join(ROOT, "output", "**", os.path.basename(r)),
                                    recursive=True))
        out += hits
    seen, uniq = set(), []
    for h in out:
        if h not in seen and os.path.isfile(h):
            seen.add(h); uniq.append(h)
    keep = [u for u in uniq if not any(x in u for x in EXCLUDE)]
    if keep: uniq = keep
    if any("/matched/" in u for u in uniq):
        uniq = [u for u in uniq if "/strict/" not in u]
    return uniq[:MAX_PER_REF]

def thumb(path, width):
    """downscale + re-encode; returns (mime, bytes)."""
    photo = "snap" in os.path.basename(path)
    ext = "jpg" if photo else "png"
    with tempfile.NamedTemporaryFile(suffix="." + ext, delete=False) as tf:
        tmp = tf.name
    cmd = ["convert", path, "-resize", f"{width}>", "-strip"]
    if photo:
        cmd += ["-quality", str(JPEG_Q)]
    cmd.append(tmp)
    subprocess.run(cmd, check=True, capture_output=True)
    data = open(tmp, "rb").read(); os.unlink(tmp)
    return ("image/jpeg" if photo else "image/png"), data

def strip_block(html, tag):
    return re.sub(re.escape(MARK_A) + tag + r".*?" + re.escape(MARK_B), "", html, flags=re.S)

def img_html(path, width, tag):
    mime, data = thumb(path, width)
    b64 = base64.b64encode(data).decode()
    rel = os.path.relpath(path, ROOT)
    return (f'<figure class="pv"><img alt="{os.path.basename(path)}" '
            f'src="data:{mime};base64,{b64}"><figcaption>{rel}</figcaption></figure>'), len(data)

def main():
    html = open(HTML, encoding="utf-8").read()
    html = strip_block(html, "")                       # refresh: drop any previous embed
    if "pv-grid" not in html:                          # inject the gallery CSS once
        html = html.replace("  hr{border:0;", """  .pv-grid{display:flex;flex-wrap:wrap;gap:.7rem;margin:.9rem 0 1.1rem}
  .pv{margin:0;flex:1 1 30rem;max-width:100%;background:#fff;border:1px solid var(--rule);
      border-radius:5px;padding:.45rem}
  .pv img{width:100%;height:auto;display:block;border-radius:3px}
  .pv figcaption{font:0.68rem/1.35 ui-monospace,SFMono-Regular,Menlo,monospace;
      color:var(--ink3);margin-top:.35rem;word-break:break-all}
  .pv-note{font-size:.78rem;color:var(--ink3);margin:.2rem 0 .8rem}
  hr{border:0;""")

    total, missing, n = 0, [], 0
    # ---- main figure blocks: insert a preview strip between the panel list and the legend
    def do_fig(m):
        nonlocal total, n
        block = m.group(0)
        refs = re.findall(r"<code>([^<]*\.png[^<]*)</code>", block)
        files = []
        for r in dict.fromkeys(refs):
            hits = expand(r)
            if not hits:
                missing.append(r)
            files += hits
        files = list(dict.fromkeys(files))[:MAX_PER_BLOCK]
        if not files:
            return block
        parts = []
        for f in files:
            h, sz = img_html(f, W_MAIN, ""); parts.append(h); total += sz; n += 1
        strip = (MARK_A + '<div class="pv-grid">' + "".join(parts) + "</div>"
                 '<p class="pv-note">Previews are downscaled; the paths above each caption are '
                 'the full-resolution originals.</p>' + MARK_B)
        return block.replace('<p class="legend">', strip + '<p class="legend">', 1)

    html = re.sub(r'<div class="fig">.*?</div>', do_fig, html, flags=re.S)

    # ---- supplementary: one gallery after the table, at a smaller width
    srefs = []
    mt = re.search(r"<h3>Supplementary figures</h3>.*?</table>\s*</div>", html, re.S)
    if mt:
        for r in dict.fromkeys(re.findall(r"<code>([^<]*\.png[^<]*)</code>", mt.group(0))):
            hits = expand(r)
            if not hits: missing.append(r)
            srefs += hits[:1]
        srefs = list(dict.fromkeys(srefs))
        if srefs:
            parts = []
            for f in srefs:
                h, sz = img_html(f, W_SUPP, ""); parts.append(h); total += sz; n += 1
            gal = (MARK_A + '<h4>Supplementary panel previews</h4><div class="pv-grid">'
                   + "".join(parts) + "</div>" + MARK_B)
            html = html[:mt.end()] + gal + html[mt.end():]

    open(HTML, "w", encoding="utf-8").write(html)
    print(f"embedded {n} images, {total/1e6:.1f} MB of image data")
    print(f"page is now {os.path.getsize(HTML)/1e6:.1f} MB")
    if missing:
        print("COULD NOT RESOLVE (skipped):")
        for r in dict.fromkeys(missing): print("   ", r)

if __name__ == "__main__":
    main()
