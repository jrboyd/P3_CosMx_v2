## Render scripts/17_summary_report.Rmd -> output/17_summary/step17_PD1_PDL1.html
## root.dir is the project root so every rd()/img() path in the Rmd is project-relative.
## NOTE: the Rmd pulls in scripts/lib/methods_notes.Rmd as a knitr child, so that path is
## resolved relative to knit_root_dir too -- do not run this from another directory.
dir.create("output/17_summary", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render("scripts/17_summary_report.Rmd",
                  output_file = "step17_PD1_PDL1.html",
                  output_dir  = "output/17_summary",
                  knit_root_dir = normalizePath("."),
                  quiet = FALSE)
cat("rendered output/17_summary/step17_PD1_PDL1.html\n")
