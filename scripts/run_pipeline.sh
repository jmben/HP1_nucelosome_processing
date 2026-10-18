#!/bin/bash
# =============================================================================
# Nucleosome Condensability Pipeline
# =============================================================================
# Samples:
#   20U input:       Dennis14_S846_L007_R1_paired_001.bam
#   20U supernatant rep1: Dennis8_S844_L007_R1_paired_001.bam
#   20U supernatant rep2: Dennis7_S843_L007_R1_paired_001.bam
#
# Requirements: samtools, deeptools (bamCoverage), bedtools, python3
#               (numpy, scipy, pandas, pyBigWig)
# =============================================================================

set -euo pipefail

# ── USER CONFIG ───────────────────────────────────────────────────────────────
BAM_DIR="${BAM_DIR:-./bams}"           # directory containing your BAM files
CHROM_SIZES="${CHROM_SIZES:-./chm13v2.0.chrom.sizes}"  # T2T chrom sizes file
OUTDIR="${OUTDIR:-./condensability_output}"
THREADS="${THREADS:-8}"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "$OUTDIR"/{coverage,nucleosomes,scores,bigwig,logs}

echo "================================================"
echo " Nucleosome Condensability Pipeline"
echo " Output: $OUTDIR"
echo " Threads: $THREADS"
echo "================================================"

# =============================================================================
# STEP 1: Generate coverage BEDGraphs from BAM files
# =============================================================================
echo ""
echo "[STEP 1] Generating coverage profiles..."

declare -A SAMPLES=(
    ["20U_input"]="Dennis14_S846_L007_R1_paired_001"
    ["20U_supernatant_rep1"]="Dennis8_S844_L007_R1_paired_001"
    ["20U_supernatant_rep2"]="Dennis7_S843_L007_R1_paired_001"
)

for LABEL in "${!SAMPLES[@]}"; do
    BAM="${BAM_DIR}/${SAMPLES[$LABEL]}.bam"
    BEDGRAPH="$OUTDIR/coverage/${LABEL}.bedgraph"
    BIGWIG="$OUTDIR/coverage/${LABEL}.bw"

    if [[ ! -f "$BAM" ]]; then
        echo "  WARNING: BAM not found: $BAM — skipping $LABEL"
        continue
    fi

    echo "  Processing $LABEL..."
    bamCoverage \
        --bam "$BAM" \
        --outFileName "$BIGWIG" \
        --outFileFormat bigwig \
        --binSize 1 \
        --normalizeUsing CPM \
        --extendReads \
        --ignoreDuplicates \
        --numberOfProcessors "$THREADS" \
        2>"$OUTDIR/logs/${LABEL}_bamCoverage.log"

    # Also produce raw (non-normalised) bedgraph for condensability calculation
    bamCoverage \
        --bam "$BAM" \
        --outFileName "$BEDGRAPH" \
        --outFileFormat bedgraph \
        --binSize 1 \
        --normalizeUsing None \
        --extendReads \
        --ignoreDuplicates \
        --numberOfProcessors "$THREADS" \
        2>>"$OUTDIR/logs/${LABEL}_bamCoverage.log"

    echo "    -> $BIGWIG"
    echo "    -> $BEDGRAPH"
done

# =============================================================================
# STEP 2 & 3: Call nucleosome peaks, compute condensability scores
# =============================================================================
echo ""
echo "[STEP 2+3] Calling nucleosome peaks and computing condensability scores..."

INPUT_BG="$OUTDIR/coverage/20U_input.bedgraph"

for REP in "rep1" "rep2"; do
    FRACTION_BG="$OUTDIR/coverage/20U_supernatant_${REP}.bedgraph"

    if [[ ! -f "$INPUT_BG" || ! -f "$FRACTION_BG" ]]; then
        echo "  WARNING: Missing bedgraph for 20U supernatant ${REP} — skipping"
        continue
    fi

    echo "  Computing condensability: 20U supernatant ${REP} vs input..."

    python3 condensability.py \
        --input-bedgraph    "$INPUT_BG" \
        --fraction-bedgraph "$FRACTION_BG" \
        --chrom-sizes       "$CHROM_SIZES" \
        --titration         "20U" \
        --fraction          "supernatant_${REP}" \
        --outdir            "$OUTDIR" \
        --threads           "$THREADS" \
        2>"$OUTDIR/logs/20U_supernatant_${REP}_condensability.log"
done

echo ""
echo "================================================"
echo " Pipeline complete!"
echo " Outputs in: $OUTDIR"
echo "   coverage/     — per-sample CPM bigWigs and raw BEDGraphs"
echo "   nucleosomes/  — peak calls per titration"
echo "   bigwig/       — condensability BEDGraph + bigWig tracks"
echo "================================================"
