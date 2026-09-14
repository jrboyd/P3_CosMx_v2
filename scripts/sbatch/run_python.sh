#!/usr/bin/env bash
#SBATCH --job-name=py-adhoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=logs/pyadhoc_%j.out
#SBATCH --error=logs/pyadhoc_%j.err

# Generic sbatch wrapper for the imaging-side Python work (tifffile / zarr /
# imagecodecs / pyvips). The counterpart of scripts/sbatch/run_rscript.sh.
#
# WHY A WRAPPER. None of this stack is on the default PATH and there is no
# python/vips/libtiff module on this cluster -- the only place it exists is the
# cosmx-img conda env below. Invoking that interpreter by absolute path is enough
# for tifffile/numpy/zarr/imagecodecs; PATH is prepended as well so the env's
# vips/tiffinfo/tiffdump CLIs are reachable from the same job.
#
# Usage, from the project root:
#   sbatch scripts/sbatch/run_python.sh path/to/script.py [args...]
#   sbatch scripts/sbatch/run_python.sh -c 'import tifffile; print(tifffile.__version__)'
#
# Override the defaults above on the command line as needed, e.g.
#   sbatch --mem=64G --time=4:00:00 scripts/sbatch/run_python.sh foo.py
#
# NOTE: compute nodes cannot see the login node's node-local /tmp, so the script
# passed here (and anything it reads or writes) must live on shared GPFS.

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch $0 <script.py> [args...] | sbatch $0 -c '<python code>'" >&2
  exit 1
fi

CONDA_ENV="/gpfs1/home/j/r/jrboyd/miniforge3/envs/cosmx-img"
PYTHON="$CONDA_ENV/bin/python"
[[ -x "$PYTHON" ]] || { echo "no interpreter at $PYTHON" >&2; exit 1; }

export PATH="$CONDA_ENV/bin:$PATH"
# The env's ImageMagick loader module is broken (libMagickCore .so missing), which
# makes vips print a harmless warning on startup. Nothing here uses the magick
# path -- native tiffload is what reads the OME-TIFFs.
export VIPS_WARNING=0

mkdir -p logs

exec "$PYTHON" "$@"
