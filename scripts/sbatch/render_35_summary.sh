#!/bin/bash
#SBATCH --job-name=35_render
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --output=logs/35_render_%j.out
#SBATCH --error=logs/35_render_%j.err
# Rebuild the downscaled figure copies, then render scripts/35_summary_report.Rmd.
set -euo pipefail
cd /gpfs1/pi/avolaric/files_jrboyd/P3_CosMx_v2
export R_LIBS_USER=/gpfs1/home/j/r/jrboyd/R/x86_64-pc-linux-gnu-library/4.5
/gpfs1/home/j/r/jrboyd/miniforge3/envs/cosmx-img/bin/python scripts/35_report_figures.py
/gpfs1/sw/rh9/pkgs/stacks/gcc/13.3.0/R/4.5.1/bin/Rscript --vanilla scripts/sbatch/render_35_summary.R
