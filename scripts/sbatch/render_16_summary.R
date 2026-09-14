## Render scripts/16_summary_report.Rmd -> output/16_summary/step16_FL_comparison.html
## root.dir is the project root so every rd()/img() path in the Rmd is project-relative.
dir.create("output/16_summary", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render("scripts/16_summary_report.Rmd",
                  output_file = "step16_FL_comparison.html",
                  output_dir  = "output/16_summary",
                  knit_root_dir = normalizePath("."),
                  quiet = FALSE)
cat("rendered output/16_summary/step16_FL_comparison.html\n")
