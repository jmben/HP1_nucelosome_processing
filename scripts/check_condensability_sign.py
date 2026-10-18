#!/usr/bin/env python3
"""
check_condensability_sign.py
============================
Verifies that a condensability bigWig has the correct sign convention
for use as the φ vector in chromatin polymer simulations.

Convention:
    positive values  →  B compartment (condensed / heterochromatin)
    negative values  →  A compartment (accessible / euchromatin)

The script computes a condensability-like score S from a numerator and
denominator bigWig:

    S = ln( (numerator + 1) / (denominator + 1) )

and reports the Pearson r between S and a compartment reference track.
A negative r means the track is A-enriched (needs inversion before use);
a positive r means B-enriched (correct orientation for the simulation).

Usage:
    python check_condensability_sign.py \\
        --numerator   20U_supernatant_rep1_1kb_inverted.bw \\
        --denominator 20U_input_1kb.bw \\
        --compartment MCF10A_subcompartments_20000_T2T.bw \\
        --binsize     1000000

    # Or check an already-computed track directly:
    python check_condensability_sign.py \\
        --track       20U_supernatant_rep1_1kb_inverted.bw \\
        --compartment MCF10A_subcompartments_20000_T2T.bw \\
        --binsize     1000000
"""

import argparse
import sys
import numpy as np
import pyBigWig
from scipy.stats import pearsonr, spearmanr

# Autosomes + chrX of T2T CHM13v2.0
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def fetch_binned(bw_path: str, chrom: str, chrom_len: int, binsize: int) -> np.ndarray:
    """Fetch mean signal per bin for one chromosome."""
    n_bins = chrom_len // binsize
    if n_bins == 0:
        return np.array([])
    bw = pyBigWig.open(bw_path)
    vals = np.array(
        bw.stats(chrom, 0, n_bins * binsize, type="mean", nBins=n_bins, exact=False),
        dtype=float,
    )
    bw.close()
    return vals


def get_chrom_sizes(bw_path: str) -> dict:
    bw = pyBigWig.open(bw_path)
    sizes = dict(bw.chroms())
    bw.close()
    return sizes


def compute_score(numerator_vals: np.ndarray, denominator_vals: np.ndarray) -> np.ndarray:
    """S = ln( (numerator + 1) / (denominator + 1) )"""
    return np.log((numerator_vals + 1.0) / (denominator_vals + 1.0))


def parse_args():
    p = argparse.ArgumentParser(
        description="Check condensability bigWig sign convention vs compartment track"
    )
    # Option A: pass numerator + denominator bigWigs to compute S on the fly
    p.add_argument("--numerator",   help="Numerator bigWig (e.g. 20U_supernatant_rep1_1kb_inverted.bw)")
    p.add_argument("--denominator", help="Denominator bigWig (e.g. 20U_input_1kb.bw)")
    # Option B: pass an already-computed condensability track
    p.add_argument("--track",       help="Pre-computed condensability bigWig to check directly")
    # Required
    p.add_argument("--compartment", required=True,
                   help="Compartment reference bigWig (positive=A convention, e.g. MCF10A_subcompartments_20000_T2T.bw)")
    p.add_argument("--binsize",     type=int, default=1_000_000,
                   help="Bin size in bp for comparison (default: 1 Mb)")
    return p.parse_args()


def main():
    args = parse_args()

    # Validate inputs
    if args.track is None and (args.numerator is None or args.denominator is None):
        print("ERROR: provide either --track or both --numerator and --denominator")
        sys.exit(1)

    mode = "track" if args.track else "ratio"

    # Get chromosome sizes from compartment reference
    comp_sizes = get_chrom_sizes(args.compartment)

    all_track_vals = []
    all_comp_vals  = []

    print(f"\nBin size   : {args.binsize:,} bp ({args.binsize // 1_000_000} Mb)")
    if mode == "ratio":
        print(f"Numerator  : {args.numerator}")
        print(f"Denominator: {args.denominator}")
        print(f"Score      : S = ln((numerator+1)/(denominator+1))")
    else:
        print(f"Track      : {args.track}")
    print(f"Compartment: {args.compartment}")
    print(f"            (positive=A, negative=B convention assumed)")
    print()

    for chrom in CHROMS:
        if chrom not in comp_sizes:
            continue
        chrom_len = comp_sizes[chrom]

        comp_vals = fetch_binned(args.compartment, chrom, chrom_len, args.binsize)
        if len(comp_vals) == 0:
            continue

        if mode == "ratio":
            num_vals  = fetch_binned(args.numerator,   chrom, chrom_len, args.binsize)
            den_vals  = fetch_binned(args.denominator, chrom, chrom_len, args.binsize)
            # Align lengths
            n = min(len(comp_vals), len(num_vals), len(den_vals))
            comp_vals = comp_vals[:n]
            num_vals  = num_vals[:n]
            den_vals  = den_vals[:n]
            track_vals = compute_score(num_vals, den_vals)
        else:
            track_vals = fetch_binned(args.track, chrom, chrom_len, args.binsize)
            n = min(len(comp_vals), len(track_vals))
            comp_vals  = comp_vals[:n]
            track_vals = track_vals[:n]

        # Keep only bins with finite values in both tracks
        mask = np.isfinite(track_vals) & np.isfinite(comp_vals)
        if mask.sum() < 3:
            continue

        all_track_vals.append(track_vals[mask])
        all_comp_vals.append(comp_vals[mask])

    if not all_track_vals:
        print("ERROR: no valid bins found — check chromosome names match between files.")
        sys.exit(1)

    track_all = np.concatenate(all_track_vals)
    comp_all  = np.concatenate(all_comp_vals)

    pr, pp = pearsonr(track_all, comp_all)
    sr, sp = spearmanr(track_all, comp_all)

    print(f"{'Genome-wide comparison':}")
    print(f"  n bins     : {len(track_all):,}")
    print(f"  Pearson  r : {pr:+.4f}  (p={pp:.2e})")
    print(f"  Spearman r : {sr:+.4f}  (p={sp:.2e})")
    print()

    # Interpret result
    # Compartment reference: positive=A (euchromatin), negative=B (heterochromatin)
    # For the simulation, we want: positive φ = B-enriched (condensed)
    # → track should correlate NEGATIVELY with the A-positive compartment reference
    if pr < -0.1:
        print("✓  CORRECT orientation:")
        print("   High track values co-occur with B compartment (negative compartment score).")
        print("   This track can be used directly as φ (positive = condensed = B).")
        print("   No inversion needed.")
    elif pr > 0.1:
        print("✗  INVERTED orientation:")
        print("   High track values co-occur with A compartment (positive compartment score).")
        print("   This track needs sign inversion before use as φ.")
        print("   Run: python invert_bigwig.py input.bw output_inverted.bw")
    else:
        print("?  WEAK correlation — track may be noisy or at wrong resolution.")
        print("   Check that both files use the same genome assembly (T2T CHM13v2.0).")

    print()


if __name__ == "__main__":
    main()
