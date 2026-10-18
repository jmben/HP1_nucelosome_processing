#!/bin/bash
# =============================================================================
# 00_setup_env.sh
# Creates the conda environment for chromatin polymer simulations.
#
# Usage:
#   bash 00_setup_env.sh
#
# After running, activate with:
#   conda activate chromatin_sim
# =============================================================================

set -euo pipefail

ENV_NAME="chromatin_sim"

echo "============================================="
echo " Setting up conda environment: $ENV_NAME"
echo "============================================="

# ── Check conda is available ──────────────────────────────────────────────────
if ! command -v conda &>/dev/null; then
    echo ""
    echo "ERROR: conda not found. Install Miniconda first:"
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh"
    echo "  bash miniconda.sh -b -p \$HOME/miniconda3"
    echo "  eval \"\$(\$HOME/miniconda3/bin/conda shell.bash hook)\""
    exit 1
fi

# ── Create environment ────────────────────────────────────────────────────────
echo ""
echo "[1/3] Creating conda environment (python=3.10)..."
conda create -y -n "$ENV_NAME" python=3.10

# ── Install packages ──────────────────────────────────────────────────────────
echo ""
echo "[2/3] Installing packages (openmm, numpy, scipy, matplotlib, tqdm)..."
conda install -y -n "$ENV_NAME" \
    -c conda-forge \
    openmm \
    numpy \
    scipy \
    matplotlib \
    tqdm

echo ""
echo "[3/3] Installing pyBigWig via pip..."
conda run -n "$ENV_NAME" pip install pyBigWig

# ── Verify installation ───────────────────────────────────────────────────────
echo ""
echo "============================================="
echo " Verifying installation..."
echo "============================================="

conda run -n "$ENV_NAME" python3 - <<'PYCHECK'
import sys

checks = [
    ("openmm",    "import openmm as mm; print(f'  OpenMM {mm.__version__}')"),
    ("pyBigWig",  "import pyBigWig; print('  pyBigWig OK')"),
    ("numpy",     "import numpy as np; print(f'  numpy {np.__version__}')"),
    ("scipy",     "import scipy; print(f'  scipy {scipy.__version__}')"),
    ("matplotlib","import matplotlib; print(f'  matplotlib {matplotlib.__version__}')"),
    ("tqdm",      "import tqdm; print(f'  tqdm {tqdm.__version__}')"),
]

failed = []
for name, code in checks:
    try:
        exec(code)
    except Exception as e:
        print(f"  MISSING: {name} — {e}")
        failed.append(name)

if not failed:
    print("\nAll packages installed successfully.")
else:
    print(f"\n{len(failed)} package(s) failed: {', '.join(failed)}")
    sys.exit(1)
PYCHECK

# ── Check GPU support ─────────────────────────────────────────────────────────
echo ""
echo "Checking OpenMM platforms (GPU availability)..."
conda run -n "$ENV_NAME" python3 -c "
import openmm as mm
platforms = [mm.Platform.getPlatform(i).getName()
             for i in range(mm.Platform.getNumPlatforms())]
print('  Available platforms:', ', '.join(platforms))
if 'CUDA' in platforms:
    print('  ✓ CUDA platform available — GPU simulations supported')
elif 'OpenCL' in platforms:
    print('  ✓ OpenCL platform available — GPU simulations supported')
else:
    print('  ℹ  No GPU platform found — simulations will run on CPU')
    print('     For CUDA support, install the matching cudatoolkit version.')
"

echo ""
echo "============================================="
echo " Setup complete!"
echo " Activate with:  conda activate $ENV_NAME"
echo "============================================="
