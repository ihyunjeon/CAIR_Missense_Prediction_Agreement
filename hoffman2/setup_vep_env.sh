#!/bin/bash
# One-time: build the `vep` conda env on Hoffman2 by cloning the known-good
# condenseq env (proven torch 2.6.0+cu118), then adding protobuf.
# Cloning rather than installing torch fresh avoids re-resolving the CUDA build.
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2   # login node is thread-capped

module load anaconda3 2>/dev/null || module load anaconda3/2020.11

if [ -d "$HOME/venvs/vep" ]; then
    echo "[skip] $HOME/venvs/vep already exists"
else
    echo "[clone] condenseq -> vep (hardlinked, same filesystem)"
    conda create --clone "$HOME/venvs/condenseq" -p "$HOME/venvs/vep" -y 2>&1 | tail -4
fi

echo "[pip] protobuf -- EsmTokenizer.from_pretrained needs it and the error message"
echo "      does not say so. NOT sentencepiece (no wheel builds on this cluster)."
"$HOME/venvs/vep/bin/pip" install -q protobuf 2>&1 | tail -3

echo "[verify]"
"$HOME/venvs/vep/bin/python" - <<'PY'
import torch, transformers
import google.protobuf as pb
print("  torch       :", torch.__version__, "| cuda", torch.version.cuda)
print("  transformers:", transformers.__version__)
print("  protobuf    :", pb.__version__)
PY
