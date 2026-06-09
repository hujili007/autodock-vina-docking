#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-docking}"
CHANNEL="${CONDA_CHANNEL:-conda-forge}"

if ! command -v conda >/dev/null 2>&1; then
  echo "Error: conda was not found. Install Miniconda, Anaconda, or Mambaforge first." >&2
  exit 1
fi

PACKAGES=(
  "python=3.11"
  "meeko"
  "openbabel"
  "vina"
  "rdkit"
  "prody"
  "pymol-open-source"
  "pillow"
)

if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
  echo "Conda environment '${ENV_NAME}' already exists; installing/updating docking packages."
  conda install -n "${ENV_NAME}" -c "${CHANNEL}" "${PACKAGES[@]}" -y
else
  echo "Creating conda environment '${ENV_NAME}'."
  conda create -n "${ENV_NAME}" -c "${CHANNEL}" "${PACKAGES[@]}" -y
fi

echo "Verifying docking environment imports."
conda run -n "${ENV_NAME}" python - <<'PY'
import prody
import rdkit
import vina
import meeko
import pymol
from PIL import Image
from openbabel import openbabel

print("Environment verification OK")
print("PyMOL:", pymol.__file__)
PY

echo "Done. Activate with: conda activate ${ENV_NAME}"
