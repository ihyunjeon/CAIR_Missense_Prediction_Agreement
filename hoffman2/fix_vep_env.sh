#!/bin/bash
# `conda create --clone` leaves a pip whose internals are inconsistent
# (ImportError: get_runnable_pip). Repair pip and add protobuf via conda,
# bypassing pip entirely.
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
module load anaconda3 2>/dev/null || module load anaconda3/2020.11
conda install -p "$HOME/venvs/vep" -y protobuf pip 2>&1 | tail -6
echo "[verify]"
"$HOME/venvs/vep/bin/python" - <<'PY'
import torch, transformers, google.protobuf as pb
print("  torch       :", torch.__version__, "| cuda", torch.version.cuda)
print("  transformers:", transformers.__version__)
print("  protobuf    :", pb.__version__)
PY
"$HOME/venvs/vep/bin/python" -m pip --version 2>&1 | tail -1
