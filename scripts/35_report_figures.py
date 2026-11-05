"""Downscaled copies of the figures scripts/35_summary_report.Rmd embeds, so the
self-contained HTML stays small. Each source PNG is resized to at most MAX_W px wide and
written to output/35_summary/fig_small/<name>. Missing sources are skipped and listed.

  python scripts/35_report_figures.py
"""
import os
import sys

import pyvips

OUT = "output/35_summary/fig_small"
MAX_W = 1500
SOURCES = {
    "fig_11_hotspot_myeloid_vs_dc.png": "output/11_hodgkin_niche/figures/fig_11_hotspot_myeloid_vs_dc.png",
    "fig29_chemokine_rank.png": "output/29_bright_cells/figures/fig29_chemokine_rank.png",
    "fig29_brightness_gene_excess.png": "output/29_bright_cells/figures/fig29_brightness_gene_excess.png",
    "fig29_cooccurrence_vs_null.png": "output/29_bright_cells/figures/fig29_cooccurrence_vs_null.png",
    "roi_example.png": "output/29_bright_cells/roi_png/high_hrs/01_c_2_39_2721.png",
    "roi_example_gray.png": "output/29_bright_cells/roi_png/high_hrs/01_c_2_39_2721_gray.png",
    "fig31_average_profile.png": "output/31_radial_cd45/figures/fig31_average_profile.png",
    "fig32_wide_profile.png": "output/32_wide_cd45/figures/fig32_wide_profile.png",
    "fig33_gene_excess_by_cluster.png": "output/33_cluster_followup/figures/fig33_gene_excess_by_cluster.png",
    "qa_polygons_candidate.png": "output/34_nuclear_size/qa_polygons_candidate.png",
    "fig35_contiguity.png": "output/35_giant_cells/figures/fig35_contiguity.png",
    "fig35_nuclei_per_member.png": "output/35_giant_cells/figures/fig35_nuclei_per_member.png",
    "fig35_rna_fragment_vs_nucleus.png": "output/35_giant_cells/figures/fig35_rna_fragment_vs_nucleus.png",
    "panels_IB1_thr20_anch.png": "output/35_giant_cells/panels_IB1_thr20_anch.png",
    "fig36_statistic_distributions.png": "output/36_molecule_dip/figures/fig36_statistic_distributions.png",
    "fig36_composite_class.png": "output/36_molecule_dip/figures/fig36_composite_class.png",
    "hotspots_miss_most_bright_cells.png": "output/37_hotspot_zoom/hotspots_miss_most_bright_cells.png",
    "IB1_57_161_zoom.png": "output/37_hotspot_zoom/IB1_57_161_zoom.png",
}


def main():
    os.makedirs(OUT, exist_ok=True)
    missing = []
    for name, src in SOURCES.items():
        if not os.path.exists(src):
            missing.append(src); continue
        im = pyvips.Image.new_from_file(src)
        if im.width > MAX_W:
            im = im.resize(MAX_W / im.width)
        im.pngsave(os.path.join(OUT, name), compression=8)
        print("wrote", name, im.width, "x", im.height, flush=True)
    if missing:
        print("missing (skipped):\n  " + "\n  ".join(missing), flush=True)


if __name__ == "__main__":
    main()
