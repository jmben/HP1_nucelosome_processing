#!/usr/bin/env python3
"""
01_bin_bigwig.py  —  Bin condensability bigwig to 10 kb beads
=============================================================================
Reads A_Inp_rep1_sensitivity.bw over chr19:11,003,182-13,275,070,
bins values to 10 kb resolution (mean of 10 bp values per bin),
normalises them into the φ vector used by the simulation energy function,
and saves a .npz file for downstream use.

Usage:
    python 01_bin_bigwig.py \
        --bw A_Inp_rep1_sensitivity.bw \
        --chrom chr19 \
        --start 11003182 \
        --end 13275070 \
        --binsize 10000 \
        --out data/phi_vector.npz
=============================================================================
"""

import argparse
import numpy as np
import pyBigWig
from pathlib import Path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Bin condensability bigwig to beads.")
    p.add_argument("--bw",      required=True,  help="Path to bigwig file")
    p.add_argument("--chrom",   default="chr19", help="Chromosome name")
    p.add_argument("--start",   type=int, default=11_003_182, help="Region start (0-based)")
    p.add_argument("--end",     type=int, default=13_275_070, help="Region end")
    p.add_argument("--binsize", type=int, default=10_000,     help="Bin size in bp (default 10 kb)")
    p.add_argument("--phi_bar",  type=float, default=0.6,      help="φ̄ centering constant (paper default 0.6)")
    p.add_argument("--clip_pct", type=float, default=1.0,      help="Percentile clip for outliers before normalisation (default 1.0)")
    p.add_argument("--flip",     action="store_true",           help="Negate raw scores before normalising (use when signal is anti-correlated with condensability, e.g. MNase sensitivity)")
    p.add_argument("--out",     default="data/phi_vector.npz", help="Output .npz file")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------
def bin_bigwig(bw_path: str, chrom: str, start: int, end: int, binsize: int) -> np.ndarray:
    """
    Fetch mean condensability score per bin.

    For each bin [b_start, b_end) we query the bigwig at full resolution
    and take the nanmean.  Bins at the right edge that are shorter than
    binsize are included (they will be slightly less precise but avoid
    dropping the last partial bin).

    Returns
    -------
    raw_scores : np.ndarray, shape (n_bins,)
        Mean condensability per bin.  NaN where no data.
    bin_starts : np.ndarray, shape (n_bins,)
        Genomic start coordinate of each bin.
    """
    bw = pyBigWig.open(bw_path)

    # Validate chromosome exists
    chroms = bw.chroms()
    if chrom not in chroms:
        raise ValueError(
            f"Chromosome '{chrom}' not found in bigwig.\n"
            f"Available: {list(chroms.keys())}"
        )

    region_len = end - start
    n_bins = int(np.ceil(region_len / binsize))
    raw_scores = np.full(n_bins, np.nan)
    bin_starts = np.zeros(n_bins, dtype=int)

    print(f"  Region  : {chrom}:{start:,}-{end:,}  ({region_len/1e3:.1f} kb)")
    print(f"  Bin size: {binsize/1e3:.0f} kb  →  {n_bins} beads")

    for i in range(n_bins):
        b_start = start + i * binsize
        b_end   = min(b_start + binsize, end)
        bin_starts[i] = b_start

        try:
            vals = bw.values(chrom, b_start, b_end, numpy=True)
            if vals is not None and len(vals) > 0:
                finite = vals[np.isfinite(vals)]
                raw_scores[i] = np.mean(finite) if len(finite) > 0 else np.nan
        except RuntimeError:
            pass  # bin outside chromosome bounds → stays NaN

    bw.close()

    n_nan = np.sum(np.isnan(raw_scores))
    if n_nan > 0:
        print(f"  Warning : {n_nan}/{n_bins} bins have no data (NaN) — filling with region mean.")
        region_mean = np.nanmean(raw_scores)
        raw_scores = np.where(np.isnan(raw_scores), region_mean, raw_scores)

    return raw_scores, bin_starts


def normalise_to_phi(raw_scores: np.ndarray, phi_bar: float = 0.6,
                     clip_pct: float = 1.0) -> np.ndarray:
    """
    Reproduce the paper's normalisation (Supplementary Note 2):

        q  = (p - min p) / (max p - min p)
        φ  = q - <q> + φ̄

    where p = raw condensability scores and φ̄ = 0.6 (empirical centering
    constant from the paper).

    clip_pct : float
        Percentile clip applied BEFORE min-max normalisation to remove
        outlier bins (e.g. mappability artifacts).  Default 1.0 clips the
        bottom 1% and top 1%.  Set to 0.0 to disable.
    """
    p = raw_scores.copy().astype(float)

    if clip_pct > 0:
        lo = np.percentile(p, clip_pct)
        hi = np.percentile(p, 100 - clip_pct)
        n_clipped = int(np.sum(p < lo) + np.sum(p > hi))
        if n_clipped > 0:
            print(f"  Clipping {n_clipped} outlier bins "
                  f"({clip_pct}–{100-clip_pct}th pct: [{lo:.4f}, {hi:.4f}])")
        p = np.clip(p, lo, hi)

    q = (p - p.min()) / (p.max() - p.min())
    phi = q - q.mean() + phi_bar
    return phi


def summarise(scores: np.ndarray, label: str):
    print(f"  {label}: min={scores.min():.4f}  max={scores.max():.4f}  "
          f"mean={scores.mean():.4f}  std={scores.std():.4f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("\n[Step 1] Reading and binning bigwig ...")
    raw_scores, bin_starts = bin_bigwig(
        args.bw, args.chrom, args.start, args.end, args.binsize
    )
    summarise(raw_scores, "raw scores (p)")

    print("\n[Step 2] Normalising to φ vector ...")
    if args.flip:
        print("  Flipping sign of raw scores (MNase sensitivity → condensability)")
        raw_scores = -raw_scores
    phi = normalise_to_phi(raw_scores, phi_bar=args.phi_bar, clip_pct=args.clip_pct)
    summarise(phi, "phi (φ)")

    print(f"\n[Step 3] Saving to {out_path} ...")
    np.savez(
        out_path,
        phi=phi,
        raw_scores=raw_scores,
        bin_starts=bin_starts,
        chrom=np.array([args.chrom]),
        region_start=np.array([args.start]),
        region_end=np.array([args.end]),
        binsize=np.array([args.binsize]),
        phi_bar=np.array([args.phi_bar]),
    )

    print(f"\n✓  Saved {len(phi)} beads  →  {out_path}\n")


if __name__ == "__main__":
    main()
