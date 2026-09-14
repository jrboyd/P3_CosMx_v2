"""Write ImageJ/Fiji .roi files and RoiSet zips -- rectangles and ovals only.

WHY THIS EXISTS. Nothing installed here writes ImageJ ROIs: RImageJROI and ijtiff
are absent, and TiffPlotR only READS them (TiffPlotR::readImageJRois). The format
is small enough to write directly, and the reader in TiffPlotR gives us a
round-trip test -- see scripts/15_verify_imagej_export.R.

THE FORMAT (ij/io/RoiEncoder.java, ij/io/RoiDecoder.java). All big-endian.
A 64-byte header, then a 64-byte "header2", then the name as UTF-16BE.

    0   4  magic "Iout"
    4   2  version              -- 227
    6   1  type                 -- 1 = rect, 2 = oval
    7   1  unused
    8   2  top     (int16)
   10   2  left    (int16)
   12   2  bottom  (int16)
   14   2  right   (int16)
   16   2  nCoordinates         -- 0 for rect/oval
   18  16  x1,y1,x2,y2 float32  -- reused as xd,yd,widthd,heightd when the
                                   SUB_PIXEL_RESOLUTION option is set
   34   2  strokeWidth
   36   4  shapeRoiSize
   40   4  strokeColor
   44   4  fillColor
   48   2  subtype
   50   2  options
   52   1  arrowStyle / aspectRatio
   53   1  arrowHeadSize
   54   2  roundedRectArcSize
   56   4  position
   60   4  header2 offset       -- 64
   64  ..  header2, offsets RELATIVE to the header2 offset:
        +4  4  C position
        +8  4  Z position
       +12  4  T position
       +16  4  name offset      -- 128
       +20  4  name length      -- in CHARACTERS, not bytes
       +24  4  overlay label colour
       +28  2  overlay font size
       +30  1  group
       +31  1  image opacity
       +32  4  image size
       +36  4  float stroke width
       +40  4  roi props offset
       +44  4  roi props length
       +48  4  counters offset
  128  ..  name, UTF-16BE

THE int16 TRAP -- the reason stage 2 works in DOWNSAMPLED or CROP-LOCAL
coordinates and never in full-resolution canvas coordinates. The bounding box at
offsets 8-15 is signed 16-bit, so it cannot express a coordinate above 32767.
IB1's canvas is 124852 x 170561 px, so a full-resolution canvas ROI is simply not
representable. We therefore:
  * write the float32 sub-pixel fields as well, with SUB_PIXEL_RESOLUTION set, so
    a reader that honours them gets the exact value; and
  * still refuse, loudly, to write anything above INT16_MAX, because the int16
    box is what ImageJ falls back on and a silent wrap puts the ROI in the wrong
    place. Downsample further, or use the .ijm macro, which has no such limit.
"""

import struct
import zipfile

HEADER_SIZE = 64
HEADER2_SIZE = 64
NAME_OFFSET = HEADER_SIZE + HEADER2_SIZE  # 128

TYPE_RECT = 1
TYPE_OVAL = 2
_TYPES = {"rect": TYPE_RECT, "oval": TYPE_OVAL}

VERSION = 227
SUB_PIXEL_RESOLUTION = 128
INT16_MAX = 32767
INT16_MIN = -32768


def _argb(color):
    """'#rrggbb' / 'rrggbb' / (r,g,b) / (a,r,g,b) / int -> packed ARGB int.

    0 means "unset", which is how ImageJ spells "use the default ROI colour", so
    an opaque black stroke has to be nudged to alpha 255 rather than left at 0.
    """
    if color is None:
        return 0
    if isinstance(color, int):
        return color
    if isinstance(color, str):
        s = color.lstrip("#")
        if len(s) == 6:
            a, r, g, b = 255, int(s[0:2], 16), int(s[2:4], 16), int(s[4:6], 16)
        elif len(s) == 8:
            a, r, g, b = (int(s[i:i + 2], 16) for i in (0, 2, 4, 6))
        else:
            raise ValueError("colour %r must be #rrggbb or #aarrggbb" % color)
    else:
        vals = tuple(int(v) for v in color)
        if len(vals) == 3:
            a, (r, g, b) = 255, vals
        elif len(vals) == 4:
            a, r, g, b = vals
        else:
            raise ValueError("colour tuple must be (r,g,b) or (a,r,g,b)")
    return ((a << 24) | (r << 16) | (g << 8) | b) & 0xFFFFFFFF


def encode_roi(x, y, width, height, shape="rect", name="roi",
               position=0, stroke_width=0, stroke_color=None):
    """One ROI as bytes. x/y is the top-left corner, in pixels."""
    try:
        roi_type = _TYPES[shape]
    except KeyError:
        raise ValueError(
            "shape must be one of %s, got %r" % (sorted(_TYPES), shape))

    if width <= 0 or height <= 0:
        raise ValueError("ROI %r has non-positive size %rx%r"
                         % (name, width, height))

    left, top = int(round(x)), int(round(y))
    right, bottom = int(round(x + width)), int(round(y + height))

    for label, v in (("left", left), ("top", top),
                     ("right", right), ("bottom", bottom)):
        if not INT16_MIN <= v <= INT16_MAX:
            raise ValueError(
                "ROI %r: %s = %d exceeds the int16 range of the ImageJ ROI "
                "format (max %d). ROIs must be written in downsampled or "
                "crop-local coordinates, not full-resolution canvas "
                "coordinates. Use the .ijm macro instead if you need these "
                "coordinates verbatim." % (name, label, v, INT16_MAX))

    hdr = bytearray(NAME_OFFSET)
    hdr[0:4] = b"Iout"
    struct.pack_into(">h", hdr, 4, VERSION)
    hdr[6] = roi_type
    hdr[7] = 0
    struct.pack_into(">hhhh", hdr, 8, top, left, bottom, right)
    struct.pack_into(">h", hdr, 16, 0)                       # nCoordinates
    # xd, yd, widthd, heightd -- honoured when SUB_PIXEL_RESOLUTION is set
    struct.pack_into(">ffff", hdr, 18, float(x), float(y),
                     float(width), float(height))
    struct.pack_into(">h", hdr, 34, int(stroke_width))
    # ">I", not ">i": an opaque colour has alpha 0xFF set, so the packed ARGB
    # exceeds signed int32. ImageJ reads the field with a signed getInt(), but
    # only the bit pattern matters, and unsigned packing produces the same bytes.
    struct.pack_into(">I", hdr, 40, _argb(stroke_color) & 0xFFFFFFFF)
    struct.pack_into(">h", hdr, 50, SUB_PIXEL_RESOLUTION)    # options
    struct.pack_into(">i", hdr, 56, int(position))
    struct.pack_into(">i", hdr, 60, HEADER_SIZE)             # header2 offset

    h2 = HEADER_SIZE
    struct.pack_into(">i", hdr, h2 + 16, NAME_OFFSET)
    struct.pack_into(">i", hdr, h2 + 20, len(name))          # CHARACTERS
    if stroke_width:
        struct.pack_into(">f", hdr, h2 + 36, float(stroke_width))

    return bytes(hdr) + name.encode("utf-16-be")


def write_roi_zip(path, rois):
    """Write a RoiSet zip that ImageJ's ROI Manager > More >> Open... accepts.

    `rois` is an iterable of dicts with keys x, y, width, height and optionally
    shape, name, position, stroke_width, color. Duplicate names get a _2, _3 ...
    suffix, because a zip with two identical entry names silently loses one in
    the ROI Manager.
    """
    seen = {}
    n = 0
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for r in rois:
            name = str(r.get("name", "roi"))
            seen[name] = seen.get(name, 0) + 1
            if seen[name] > 1:
                name = "%s_%d" % (name, seen[name])
            z.writestr(name + ".roi", encode_roi(
                r["x"], r["y"], r["width"], r["height"],
                shape=r.get("shape", "rect"), name=name,
                position=r.get("position", 0),
                stroke_width=r.get("stroke_width", 0),
                stroke_color=r.get("color")))
            n += 1
    return n


def write_macro(path, rois, title=None, note=None):
    """An .ijm macro that rebuilds the same ROIs via makeRectangle/makeOval.

    The macro route has no int16 limit and is human-readable, so it is the
    fallback whenever the binary ROIs will not do. Run it with the target image
    already open (Plugins > Macros > Run...).
    """
    lines = ["// ImageJ macro -- rebuilds the ROI set in the ROI Manager.",
             "// Open the matching image FIRST, then Plugins > Macros > Run...",
             ""]
    if title:
        lines.insert(0, "// %s" % title)
    if note:
        lines += ["// %s" % l for l in note.splitlines()] + [""]
    lines += ['if (nImages == 0) exit("Open the image first.");',
              "roiManager(\"reset\");", ""]
    for r in rois:
        maker = "makeOval" if r.get("shape") == "oval" else "makeRectangle"
        lines.append("%s(%.2f, %.2f, %.2f, %.2f);"
                     % (maker, r["x"], r["y"], r["width"], r["height"]))
        lines.append('Roi.setName("%s");' % str(r.get("name", "roi")))
        if r.get("color"):
            # the macro parser takes bare hex ("eb6834") or a colour name,
            # but not a leading "#"
            lines.append('Roi.setStrokeColor("%s");'
                         % str(r["color"]).lstrip("#"))
        lines.append('roiManager("add");')
    lines += ["", "roiManager(\"show all with labels\");",
              'print("' + str(len(rois)) + ' ROIs loaded.");', ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return len(rois)
