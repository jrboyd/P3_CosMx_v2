## Render scripts/35_summary_report.Rmd -> output/35_summary/hrs_candidates.html
## Run scripts/35_report_figures.py first (downscaled figure copies the Rmd embeds).
## Pulls in scripts/lib/methods_notes.Rmd as a knitr child, so run from the project root.
dir.create("output/35_summary", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render("scripts/35_summary_report.Rmd",
                  output_file = "hrs_candidates.html",
                  output_dir  = "output/35_summary",
                  knit_root_dir = normalizePath("."),
                  quiet = FALSE)
cat("rendered output/35_summary/hrs_candidates.html\n")
