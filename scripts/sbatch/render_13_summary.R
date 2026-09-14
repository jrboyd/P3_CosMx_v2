## Render scripts/13_summary_report.Rmd -> output/13_summary/IB1_summary.html
## Filename is stable on purpose: CLAUDE.md and both later reports point at it. The revision
## is signalled by the dated banner and subtitle inside the document, not by a new filename.
## Pulls in scripts/lib/methods_notes.Rmd as a knitr child, so run from the project root.
dir.create("output/13_summary", recursive = TRUE, showWarnings = FALSE)
rmarkdown::render("scripts/13_summary_report.Rmd",
                  output_file = "IB1_summary.html",
                  output_dir  = "output/13_summary",
                  knit_root_dir = normalizePath("."),
                  quiet = FALSE)
cat("rendered output/13_summary/IB1_summary.html\n")
