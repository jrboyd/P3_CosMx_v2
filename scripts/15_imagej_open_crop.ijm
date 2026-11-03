// Normalise an IB1 crop (output/15_imagej/crops/*.tif) in Fiji: 5-channel
// composite, house LUTs, and the ONE display range all 34 crops share.
//
// A plain File > Open of these TIFFs already gives all of this -- the header
// declares channels=5 / mode=composite and embeds the LUTs and ranges below.
// This macro is for when an opener ignores that (e.g. Bio-Formats), which shows
// the five channels as Z slices with auto contrast.
//
// Use: open a crop, then Plugins > Macros > Run... this file (or Install it to
// get a menu entry). Channel order is positional and matches
// cosmx_channel_names() in scripts/lib/cosmx_tiff.R.
//
// Ranges are the shared crop ranges embedded by scripts/15_imagej_tiles.py
// (30th/99th pct of non-zero level-0 pixels). They are NOT valid for
// IB1_overview_L5.tif, which has its own lower ranges -- re-open that one plainly.

names = newArray("PanCK", "G", "Membrane", "CD45", "DNA");
luts  = newArray("Blue", "Green", "Yellow", "Red", "Magenta");  // Magenta stands in for violet
lo    = newArray(220, 120,  964,  752,  924);
hi    = newArray(1428, 352, 8464, 9924, 7304);

getDimensions(w, h, nc, nz, nt);
if (w != 4256 || h != 4256)
    exit("Not a 4256 px crop (" + w + " x " + h + "); these ranges are for the crops only.");
if (nc * nz * nt != 5)
    exit("Expected 5 planes, found " + (nc * nz * nt) + " -- this is not an intact crop TIFF.");
if (nc != 5)
    run("Stack to Hyperstack...", "order=xyczt(default) channels=5 slices=1 frames=1 display=Composite");

Stack.setDisplayMode("composite");
for (i = 0; i < 5; i++) {
    Stack.setChannel(i + 1);
    run(luts[i]);
    setMinAndMax(lo[i], hi[i]);
    setMetadata("Label", names[i]);
}
Stack.setChannel(1);
