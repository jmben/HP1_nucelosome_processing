#!/usr/bin/env python3
"""
condensability.py
=================
Computes genome-wide nucleosome condensability scores from:
  - Input (pre-condensation) coverage BEDGraph
  - Supernatant coverage BEDGraph

Steps:
  1. Load coverage into per-chromosome arrays
  2. Call nucleosome peaks on input (local maxima, max 40 bp overlap between 147 bp windows)
  3. Per-nucleosome condensability: -ln((supernatant_area + 1) / (input_area + 1)) in 171 bp window
  4. Sliding-window condensability: 171 bp window, 25 bp step
  5. Binned condensability: 1 kb and 10 kb bins
  6. Write BEDGraph + bigWig for each resolution
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import pyBigWig

# ── Constants (from the paper) ────────────────────────────────────────────────
NUC_SIZE        = 147   # nucleosome footprint (bp)
WINDOW_SIZE     = 171   # coverage window for area calculation (bp)
MAX_OVERLAP     = 40    # maximum allowed overlap between 147 bp windows
SLIDE_STEP      = 25    # sliding window step size (bp)
BIN_1KB         = 1_000
BIN_10KB        = 10_000
PSEUDOCOUNT     = 1     # added to both numerator and denominator before log
# ─────────────────────────────────────────────────────────────────────────────


def parse_args():
    p = argparse.ArgumentParser(description="Compute nucleosome condensability scores")
    p.add_argument("--input-bedgraph",    required=True, help="Input control BEDGraph (raw counts)")
    p.add_argument("--fraction-bedgraph", required=True, help="Supernatant BEDGraph (raw counts)")
    p.add_argument("--chrom-sizes",       required=True, help="Chromosome sizes file (2-column TSV)")
    p.add_argument("--titration",         required=True, help="Titration label, e.g. 20U")
    p.add_argument("--fraction",          required=True, help="Fraction label: supernatant")
    p.add_argument("--outdir",            required=True, help="Output directory")
    p.add_argument("--threads",           type=int, default=4)
    return p.parse_args()


# =============================================================================
# I/O helpers
# =============================================================================

def load_chrom_sizes(path):
    """Returns OrderedDict: chrom -> size"""
    sizes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def bedgraph_to_arrays(path, chrom_sizes):
    """
    Parses a BEDGraph into per-chromosome numpy float32 arrays.
    Only chromosomes present in chrom_sizes are kept.
    """
    print(f"  Loading BEDGraph: {path}", flush=True)
    arrays = {chrom: np.zeros(size, dtype=np.float32)
              for chrom, size in chrom_sizes.items()}

    with open(path) as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#")):
                continue
            parts = line.split()
            chrom, start, end, val = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            if chrom in arrays:
                arrays[chrom][start:end] = val

    return arrays


# =============================================================================
# Nucleosome peak calling
# =============================================================================

def call_nucleosome_peaks(coverage_array, chrom, chrom_size):
    """
    Calls nucleosome peak positions on a single chromosome.

    Algorithm (faithful to methods):
      1. Find all local maxima in the coverage profile.
      2. Starting from a random peak, greedily select peaks in both
         directions, rejecting any peak whose 147 bp window overlaps
         an already-selected peak by more than 40 bp.

    Returns sorted array of peak centre positions.
    """
    # find_peaks returns indices of local maxima
    peaks, _ = find_peaks(coverage_array, distance=NUC_SIZE - MAX_OVERLAP)

    if len(peaks) == 0:
        return np.array([], dtype=np.int64)

    # Greedy selection with overlap filter
    # Sort by position (find_peaks already returns sorted indices)
    selected = []
    last_end = -1  # right edge of last accepted 147 bp window

    for p in peaks:
        win_start = p - NUC_SIZE // 2
        win_end   = p + NUC_SIZE // 2

        if win_start < 0 or win_end >= chrom_size:
            continue

        overlap = max(0, last_end - win_start)
        if overlap <= MAX_OVERLAP:
            selected.append(p)
            last_end = win_end

    return np.array(selected, dtype=np.int64)


# =============================================================================
# Coverage area calculation
# =============================================================================

def coverage_area(cov_array, centres, half_win, chrom_size):
    """
    For each centre position, sums coverage in [centre-half_win, centre+half_win].
    Returns float32 array of areas.
    """
    areas = np.empty(len(centres), dtype=np.float64)
    for i, c in enumerate(centres):
        s = max(0, c - half_win)
        e = min(chrom_size, c + half_win + 1)
        areas[i] = cov_array[s:e].sum()
    return areas


# =============================================================================
# Condensability score: -ln( (fraction + pseudo) / (input + pseudo) )
# =============================================================================

def condensability_score(input_areas, fraction_areas):
    """
    score = +ln( (fraction_area + 1) / (input_area + 1) )
    Positive → nucleosome enriched in fraction (high condensability)
    Negative → nucleosome depleted in fraction
    Matches published convention (e.g. GM NCP condensability tracks).
    """
    return np.log((fraction_areas + PSEUDOCOUNT) / (input_areas + PSEUDOCOUNT))


# =============================================================================
# Write BEDGraph
# =============================================================================

def write_bedgraph(path, records):
    """records: list of (chrom, start, end, score)"""
    with open(path, "w") as fh:
        fh.write('track type=bedGraph\n')
        for chrom, start, end, score in records:
            fh.write(f"{chrom}\t{start}\t{end}\t{score:.6f}\n")
    print(f"  Written: {path}", flush=True)


def bedgraph_to_bigwig(bedgraph_path, bw_path, chrom_sizes_path):
    """
    Convert BEDGraph to bigWig by:
      1. Sorting the BEDGraph by chrom then position (required by bedGraphToBigWig)
      2. Calling UCSC bedGraphToBigWig
    Falls back to pyBigWig if bedGraphToBigWig is not in PATH.
    """
    import subprocess, shutil, tempfile

    sorted_bg = bedgraph_path + ".sorted.tmp"
    try:
        # Sort: skip track header line, sort by chrom then start
        with open(bedgraph_path) as fh:
            lines = [l for l in fh if not l.startswith("track")]
        lines.sort(key=lambda l: (l.split()[0], int(l.split()[1])))
        with open(sorted_bg, "w") as fh:
            fh.writelines(lines)

        if shutil.which("bedGraphToBigWig"):
            subprocess.run(
                ["bedGraphToBigWig", sorted_bg, chrom_sizes_path, bw_path],
                check=True, capture_output=True
            )
        else:
            # pyBigWig fallback — feed entries chromosome by chromosome
            import pyBigWig
            # Read chrom sizes from file
            chrom_sizes_dict = {}
            with open(chrom_sizes_path) as fh:
                for line in fh:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        chrom_sizes_dict[parts[0]] = int(parts[1])

            bw = pyBigWig.open(bw_path, "w")
            bw.addHeader(list(chrom_sizes_dict.items()))

            current_chrom = None
            chroms, starts, ends, vals = [], [], [], []

            def flush_chrom(bw, chroms, starts, ends, vals):
                if chroms:
                    bw.addEntries(chroms, starts, ends=ends, values=vals)

            with open(sorted_bg) as fh:
                for line in fh:
                    parts = line.split()
                    c, s, e, v = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
                    if c != current_chrom:
                        flush_chrom(bw, chroms, starts, ends, vals)
                        chroms, starts, ends, vals = [], [], [], []
                        current_chrom = c
                    chroms.append(c); starts.append(s); ends.append(e); vals.append(v)
            flush_chrom(bw, chroms, starts, ends, vals)
            bw.close()

        print(f"  Written: {bw_path}", flush=True)
    finally:
        if os.path.exists(sorted_bg):
            os.remove(sorted_bg)


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    label = f"{args.titration}_{args.fraction}"
    nuc_dir   = os.path.join(args.outdir, "nucleosomes")
    score_dir = os.path.join(args.outdir, "scores")
    bw_dir    = os.path.join(args.outdir, "bigwig")
    for d in [nuc_dir, score_dir, bw_dir]:
        os.makedirs(d, exist_ok=True)

    # ── Load data ──────────────────────────────────────────────────────────
    chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    input_cov    = bedgraph_to_arrays(args.input_bedgraph,    chrom_sizes)
    fraction_cov = bedgraph_to_arrays(args.fraction_bedgraph, chrom_sizes)

    half_win = WINDOW_SIZE // 2  # 85 bp

    # We'll accumulate records for three resolutions
    records_per_nuc   = []   # per-nucleosome peak
    records_sliding   = []   # 171 bp / 25 bp step
    records_1kb       = []
    records_10kb      = []

    # Also save per-nucleosome peak table
    peak_rows = []

    print(f"\n  Processing {len(chrom_sizes)} chromosomes...", flush=True)

    for chrom, chrom_size in chrom_sizes.items():
        inp = input_cov[chrom]
        frac = fraction_cov[chrom]

        # ── Per-nucleosome ────────────────────────────────────────────────
        peaks = call_nucleosome_peaks(inp, chrom, chrom_size)
        if len(peaks) > 0:
            inp_areas  = coverage_area(inp,  peaks, half_win, chrom_size)
            frac_areas = coverage_area(frac, peaks, half_win, chrom_size)
            scores     = condensability_score(inp_areas, frac_areas)

            for p, ia, fa, sc in zip(peaks, inp_areas, frac_areas, scores):
                s = max(0, p - half_win)
                e = min(chrom_size, p + half_win + 1)
                records_per_nuc.append((chrom, s, e, sc))
                peak_rows.append({
                    "chrom": chrom, "start": s, "end": e,
                    "peak_pos": p,
                    "input_area": ia, "fraction_area": fa,
                    "condensability": sc
                })

        # ── Sliding window 171 bp / 25 bp step ───────────────────────────
        centres_slide = np.arange(half_win, chrom_size - half_win, SLIDE_STEP)
        if len(centres_slide) > 0:
            inp_areas  = coverage_area(inp,  centres_slide, half_win, chrom_size)
            frac_areas = coverage_area(frac, centres_slide, half_win, chrom_size)
            scores     = condensability_score(inp_areas, frac_areas)
            for c, sc in zip(centres_slide, scores):
                s = c - half_win
                e = c + half_win + 1
                records_sliding.append((chrom, s, e, sc))

        # ── 1 kb bins ────────────────────────────────────────────────────
        for bin_start in range(0, chrom_size, BIN_1KB):
            bin_end = min(chrom_size, bin_start + BIN_1KB)
            ia  = inp[bin_start:bin_end].sum()
            fa  = frac[bin_start:bin_end].sum()
            sc  = condensability_score(ia, fa)
            records_1kb.append((chrom, bin_start, bin_end, float(sc)))

        # ── 10 kb bins ───────────────────────────────────────────────────
        for bin_start in range(0, chrom_size, BIN_10KB):
            bin_end = min(chrom_size, bin_start + BIN_10KB)
            ia  = inp[bin_start:bin_end].sum()
            fa  = frac[bin_start:bin_end].sum()
            sc  = condensability_score(ia, fa)
            records_10kb.append((chrom, bin_start, bin_end, float(sc)))

        print(f"    {chrom}: {len(peaks)} nucleosomes called", flush=True)

    # ── Write outputs ──────────────────────────────────────────────────────

    # Per-nucleosome peak TSV
    peak_tsv = os.path.join(nuc_dir, f"{label}_nucleosomes.tsv")
    pd.DataFrame(peak_rows).to_csv(peak_tsv, sep="\t", index=False)
    print(f"\n  Peak table: {peak_tsv}", flush=True)

    resolutions = [
        ("per_nuc",  records_per_nuc),
        ("slide171_25", records_sliding),
        ("1kb",      records_1kb),
        ("10kb",     records_10kb),
    ]

    for res_label, records in resolutions:
        bg_path = os.path.join(bw_dir, f"{label}_{res_label}.bedgraph")
        bw_path = os.path.join(bw_dir, f"{label}_{res_label}.bw")
        write_bedgraph(bg_path, records)
        try:
            bedgraph_to_bigwig(bg_path, bw_path, args.chrom_sizes)
        except Exception as e:
            print(f"  WARNING: bigWig conversion failed for {res_label}: {e}", file=sys.stderr)

    print(f"\n  Done: {label}", flush=True)


if __name__ == "__main__":
    main()
