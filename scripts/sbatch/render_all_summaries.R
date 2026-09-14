## Re-render steps 13, 16, 17 after the methods/caption pass. Each Rmd pulls in
## scripts/lib/methods_notes.Rmd as a knitr child, so that path must resolve from the root.
for (x in list(c("13_summary_report.Rmd","output/13_summary","IB1_summary.html"),
               c("16_summary_report.Rmd","output/16_summary","step16_FL_comparison.html"),
               c("17_summary_report.Rmd","output/17_summary","step17_PD1_PDL1.html"))) {
  dir.create(x[2], recursive = TRUE, showWarnings = FALSE)
  cat("\n#### rendering", x[1], "\n")
  rmarkdown::render(file.path("scripts", x[1]), output_file = x[3], output_dir = x[2],
                    knit_root_dir = normalizePath("."), quiet = FALSE)
}
cat("\nall three rendered\n")
