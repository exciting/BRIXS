#!/usr/bin/env bash
# Runs the BRIXS regression suite. Invoked by CTest (see test/CMakeLists.txt);
# can also be run directly from within test/.
#
# If a `pybrixs_test` conda environment exists (the local/HPC convention),
# it's activated first. Otherwise this just runs `python3` as found on PATH,
# which is what CI uses (h5py/numpy installed directly via pip there).
set -euo pipefail

if command -v conda >/dev/null 2>&1 && conda env list | grep -q '^pybrixs_test '; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate pybrixs_test
fi

python3 "$(dirname "$0")/test.py"
