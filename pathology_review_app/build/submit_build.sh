#!/usr/bin/env bash
# Build a pathology review bundle end to end, as a chain of SLURM jobs.
#
#   bash pathology_review_app/build/submit_build.sh test    # 50 items, ~30 min
#   bash pathology_review_app/build/submit_build.sh full    # ~217 items, ~70 min (image reads dominate)
#
# select (R, 16-24G) -> render (python, image reads from scratch) -> bundle (stdlib) ->
# headless-Firefox end-to-end test of the bundle. Run from the project root.
# Result: pathology_review_app/output/<size>/send/ holds exactly what to send the reader
# (reviewer.html, the zip, READ_ME_FIRST.txt); the key goes to pathology_review_app/keys/.
set -euo pipefail
SIZE="${1:-test}"
cd "$(git rev-parse --show-toplevel)"
B=pathology_review_app/build
mkdir -p logs

J1=$(sbatch --parsable --job-name=pr_select --mem=24G --time=0:30:00 \
  --output=logs/pr_select_%j.out --error=logs/pr_select_%j.err \
  scripts/sbatch/run_rscript.sh $B/01_select_cells.R "$SIZE")
J2=$(sbatch --parsable --dependency=afterok:$J1 --job-name=pr_render --mem=24G --time=2:00:00 \
  --output=logs/pr_render_%j.out --error=logs/pr_render_%j.err \
  scripts/sbatch/run_python.sh $B/02_render_tiles.py "$SIZE")
J3=$(sbatch --parsable --dependency=afterok:$J2 --job-name=pr_bundle --mem=4G --time=0:20:00 \
  --output=logs/pr_bundle_%j.out --error=logs/pr_bundle_%j.err \
  scripts/sbatch/run_python.sh $B/03_make_bundle.py "$SIZE")
J4=$(sbatch --parsable --dependency=afterok:$J3 --job-name=pr_apptest --mem=8G --cpus-per-task=2 --time=0:20:00 \
  --output=logs/pr_apptest_%j.out --error=logs/pr_apptest_%j.err \
  scripts/sbatch/run_python.sh $B/test_app_headless.py "$SIZE")
echo "select $J1 -> render $J2 -> bundle $J3 -> app test $J4"
echo "check: sacct -j $J1,$J2,$J3,$J4 --format=JobID,JobName,State,Elapsed,MaxRSS"
