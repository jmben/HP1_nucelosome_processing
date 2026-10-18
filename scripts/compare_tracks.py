#!/usr/bin/env python3
"""
compare_tracks.py  —  Direct comparison of bigwig tracks vs Hi-C compartments
=============================================================================
Computes per-bin correlation between each bigwig and the MCF10A compartment
track across a defined genomic region, then produces a summary plot.

Usage:
    python compare_tracks.py \
        --compartment MCF10A_subcompartments_20000_T2T.bw \
        --chrom chr2 --start 10000000 --end 14000000 \
        --binsize 20000
=============================================================================
"""

import argparse
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

try:
    import pyBigWig
except ImportError:
    raise ImportError("pip install pyBigWig")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
except ImportError:
    raise ImportError("pip install matplotlib")


# ---------------------------------------------------------------------------
# Tracks to compare — edit flip=True if track is anti-correlated with B
# ---------------------------------------------------------------------------
TRACKS = [
    # (label,                    filename,                                          flip)
    # All files are correctly oriented: positive=B-enriched (condensed)
    # No flipping needed
    ("GM_condense-seq\n(paper)", "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T.bw", False),
    ("20U supernatant rep1",     "20U_supernatant_rep1_1kb_inverted.bw",            False),
    ("Genomic rep1\nsensitivity","MCF10A_genomic_rep1_downsampled_sensitivity.bw",   False),
    ("Genomic rep2\nsensitivity","MCF10A_genomic_rep2_downsampled_sensitivity.bw",   False),
    ("A_Inp rep1\nsensitivity",  "A_Inp_rep1_sensitivity.bw",                        False),
]


# ---------------------------------------------------------------------------
# Helper: fetch binned values from bigwig
# ---------------------------------------------------------------------------
def fetch_binned(bw_path: str, chrom: str, start: int, end: int,
                 binsize: int, flip: bool = False) -> np.ndarray:
    bw    = pyBigWig.open(bw_path)
    n     = (end - start) // binsize
    vals  = np.array(
        bw.stats(chrom, start, start + n * binsize,
                 type="mean", nBins=n, exact=False),
        dtype=float
    )
    bw.close()
    if flip:
        vals = -vals
    # Fill NaN with column mean
    finite = vals[np.isfinite(vals)]
    if len(finite) > 0:
        vals = np.where(np.isfinite(vals), vals, np.nanmean(finite))
    return vals


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--compartment", default="MCF10A_subcompartments_20000_T2T.bw")
    p.add_argument("--chrom",       default="chr2")
    p.add_argument("--start",       type=int, default=10_000_000)
    p.add_argument("--end",         type=int, default=14_000_000)
    p.add_argument("--binsize",     type=int, default=20_000,
                   help="Bin size in bp (default 20kb — matches compartment track resolution)")
    p.add_argument("--out",         default="track_comparison.png")
    return p.parse_args()


def main():
    args   = parse_args()
    chrom  = args.chrom
    start  = args.start
    end    = args.end
    bs     = args.binsize
    n_bins = (end - start) // bs
    mb     = np.array([start + i * bs for i in range(n_bins)]) / 1e6

    print(f"Region : {chrom}:{start:,}-{end:,}  ({(end-start)/1e6:.1f} Mb)")
    print(f"Bins   : {n_bins}  ({bs/1000:.0f} kb each)\n")

    # Load compartment reference — use raw values (positive=B, negative=A in this file)
    print(f"Loading compartment reference: {args.compartment}")
    comp = fetch_binned(args.compartment, chrom, start, end, bs, flip=False)
    print(f"  range: {comp.min():.3f} – {comp.max():.3f}  "
          f"(positive=B, negative=A — B-enriched tracks expect positive r)\n")

    # Load and correlate each track
    results = []
    track_vals = []

    print(f"{'Track':<28} {'Pearson r':>10} {'Spearman r':>11} {'p-value':>10}  direction")
    print("-" * 70)

    for label, fname, flip in TRACKS:
        try:
            vals = fetch_binned(fname, chrom, start, end, bs, flip=flip)
        except Exception as e:
            print(f"  {label:<26} ERROR: {e}")
            track_vals.append(None)
            results.append((label, np.nan, np.nan, np.nan))
            continue

        # Mask bins where either track has no data
        mask = np.isfinite(vals) & np.isfinite(comp)
        if mask.sum() < 10:
            print(f"  {label:<26} insufficient data")
            track_vals.append(vals)
            results.append((label, np.nan, np.nan, np.nan))
            continue

        pr, pp = pearsonr(comp[mask],  vals[mask])
        sr, sp = spearmanr(comp[mask], vals[mask])

        # Expected direction: positive r = track positively correlated with A compartment
        # For condense-seq: B enriched, so negative r with compartment score is expected
        # For sensitivity (flipped): same
        direction = "+" if pr > 0 else "-"
        sig = "***" if pp < 0.001 else ("**" if pp < 0.01 else ("*" if pp < 0.05 else "ns"))
        clean_label = label.replace("\n", " ")
        print(f"  {clean_label:<26} {pr:>10.4f} {sr:>11.4f} {pp:>10.2e}  {direction} {sig}")

        track_vals.append(vals)
        results.append((label, pr, sr, pp))

    # ── Plot ──────────────────────────────────────────────────────────────
    n_tracks = len(TRACKS)
    fig = plt.figure(figsize=(14, 3 + n_tracks * 1.4))
    gs  = gridspec.GridSpec(n_tracks + 2, 1, hspace=0.15)

    # Panel 0: compartment reference (positive=B, negative=A — raw file convention)
    ax0 = fig.add_subplot(gs[0])
    ax0.fill_between(mb, comp, where=(comp >= 0),
                     color="#2c3e50", alpha=0.9, label="B compartment (+)")
    ax0.fill_between(mb, comp, where=(comp < 0),
                     color="#95a5a6", alpha=0.9, label="A compartment (-)")
    ax0.axhline(0, color="black", lw=0.5)
    ax0.set_ylabel("Compartment\n(B=dark+, A=grey-)", fontsize=7)
    ax0.legend(fontsize=6, loc="upper right")
    ax0.set_title(
        f"Track comparison — {chrom}:{start//1_000_000}–{end//1_000_000} Mb",
        fontsize=10
    )
    ax0.set_xlim(mb[0], mb[-1])
    plt.setp(ax0.get_xticklabels(), visible=False)

    # One panel per track
    colors = ["#e74c3c", "#e67e22", "#f39c12", "#27ae60", "#2980b9", "#8e44ad"]
    for i, ((label, pr, sr, pp), vals) in enumerate(zip(results, track_vals)):
        ax = fig.add_subplot(gs[i + 1], sharex=ax0)
        if vals is not None:
            # Normalise for display
            v = vals.copy()
            v = (v - np.nanmean(v)) / (np.nanstd(v) + 1e-10)
            ax.fill_between(mb, v, where=(v >= 0), color=colors[i % len(colors)],
                            alpha=0.7)
            ax.fill_between(mb, v, where=(v < 0), color=colors[i % len(colors)],
                            alpha=0.3)
            ax.axhline(0, color="black", lw=0.4)

        sig = "***" if pp < 0.001 else ("**" if pp < 0.01 else ("*" if pp < 0.05 else "ns"))
        clean = label.replace("\n", " ")
        ax.set_ylabel(f"{clean}\nr={pr:.3f}{sig}", fontsize=6.5, labelpad=2)
        ax.set_xlim(mb[0], mb[-1])
        ax.tick_params(labelsize=6)
        if i < n_tracks - 1:
            plt.setp(ax.get_xticklabels(), visible=False)

    # Panel: bar chart of Pearson r values
    ax_bar = fig.add_subplot(gs[n_tracks + 1])
    labels_clean = [r[0].replace("\n", " ") for r in results]
    rs     = [r[1] for r in results]
    pvals  = [r[3] for r in results]
    bar_colors = [colors[i % len(colors)] for i in range(len(rs))]
    bars = ax_bar.bar(range(len(rs)), rs, color=bar_colors, alpha=0.8, edgecolor="black", lw=0.5)
    ax_bar.axhline(0, color="black", lw=0.8)
    ax_bar.axhline(0.3,  color="grey", lw=0.5, ls="--", alpha=0.5)
    ax_bar.axhline(-0.3, color="grey", lw=0.5, ls="--", alpha=0.5)
    ax_bar.set_xticks(range(len(rs)))
    ax_bar.set_xticklabels(labels_clean, fontsize=7, rotation=15, ha="right")
    ax_bar.set_ylabel("Pearson r\nvs compartment\n(pos=B-enriched)", fontsize=8)
    ax_bar.set_title("Correlation with MCF10A compartment (positive=B, negative=A)", fontsize=9)

    # Significance stars on bars
    for i, (r_val, pv) in enumerate(zip(rs, pvals)):
        if np.isfinite(r_val) and np.isfinite(pv):
            sig = "***" if pv < 0.001 else ("**" if pv < 0.01 else ("*" if pv < 0.05 else ""))
            if sig:
                y = r_val + (0.02 if r_val >= 0 else -0.04)
                ax_bar.text(i, y, sig, ha="center", fontsize=8)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nSaved: {args.out}")


if __name__ == "__main__":
    main()
