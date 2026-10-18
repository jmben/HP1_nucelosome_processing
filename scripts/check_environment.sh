#!/bin/bash
# =============================================================================
# check_environment.sh
# Verifies all required tools and Python packages are available
# =============================================================================

set -euo pipefail

PASS=0
FAIL=0

check_tool() {
    local name=$1
    local cmd=$2
    if command -v "$cmd" &>/dev/null; then
        echo "  [OK]  $name ($(command -v $cmd))"
        ((PASS++)) || true
    else
        echo "  [MISSING] $name — install with: $3"
        ((FAIL++)) || true
    fi
}

check_python_pkg() {
    local pkg=$1
    if python3 -c "import $pkg" 2>/dev/null; then
        echo "  [OK]  python3: $pkg"
        ((PASS++)) || true
    else
        echo "  [MISSING] python3: $pkg — install with: pip install $pkg"
        ((FAIL++)) || true
    fi
}

echo "=============================="
echo " Environment Check"
echo "=============================="
echo ""
echo "--- Command-line tools ---"
check_tool "samtools"    samtools   "conda install -c bioconda samtools"
check_tool "bamCoverage" bamCoverage "conda install -c bioconda deeptools"
check_tool "bedtools"    bedtools   "conda install -c bioconda bedtools"
check_tool "python3"     python3    "conda install python=3.10"

echo ""
echo "--- Python packages ---"
check_python_pkg numpy
check_python_pkg scipy
check_python_pkg pandas
check_python_pkg pyBigWig

echo ""
echo "=============================="
if [[ $FAIL -eq 0 ]]; then
    echo " All checks passed ($PASS/$((PASS+FAIL)))"
else
    echo " $FAIL missing, $PASS passed"
    echo ""
    echo " Quick install (conda recommended):"
    echo "   conda create -n condensability python=3.10"
    echo "   conda activate condensability"
    echo "   conda install -c bioconda -c conda-forge samtools deeptools bedtools"
    echo "   pip install numpy scipy pandas pyBigWig"
fi
echo "=============================="
