## Render scripts/19_summary_report.Rmd -> output/19_summary/step19_cell_identity.html
## knit_root_dir is the project root so every rd()/img() path in the Rmd is project-relative.
dir.create("output/19_summary", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render("scripts/19_summary_report.Rmd",
                  output_file = "step19_cell_identity.html",
                  output_dir  = "output/19_summary",
                  knit_root_dir = normalizePath("."),
                  quiet = FALSE)
cat("rendered output/19_summary/step19_cell_identity.html\n")
