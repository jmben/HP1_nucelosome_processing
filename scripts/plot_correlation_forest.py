#!/usr/bin/env python3
"""
plot_correlation_forest.py  —  Forest plot of genome-wide track correlations
=============================================================================
Usage:
    python plot_correlation_forest.py --out correlation_forest.png
=============================================================================
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ---------------------------------------------------------------------------
# Data — paste updated values here if rerunning at different bin sizes
# ---------------------------------------------------------------------------
DATA = [
    # (label,                              pearson_r,  ci_lo,    ci_hi,   spearman_r, category)
    ("GM condense-seq\n(paper, GM12878)",  -0.5840,   -0.6076,  -0.5593,  -0.6150,  "condense-seq"),
    ("20U supernatant rep1\n(MCF10A)",     -0.5419,   -0.5673,  -0.5155,  -0.6246,  "condense-seq"),
    ("Genomic rep2\nsensitivity",          -0.4300,   -0.4594,  -0.3997,  -0.4575,  "sensitivity"),
    ("Genomic rep1\nsensitivity",          -0.2619,   -0.2956,  -0.2274,  -0.3068,  "sensitivity"),
    ("A_Inp rep1\nsensitivity",            -0.1400,   -0.1762,  -0.1035,  -0.3132,  "sensitivity"),
]

CATEGORY_COLORS = {
    "condense-seq": "#e67e22",   # orange
    "sensitivity":  "#2980b9",   # blue
}

CATEGORY_LABELS = {
    "condense-seq": "Condense-seq",
    "sensitivity":  "MNase sensitivity",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--out",      default="correlation_forest.png")
    p.add_argument("--binsize",  default="1 Mb", help="Bin size label for title")
    p.add_argument("--width",    type=float, default=8)
    p.add_argument("--height",   type=float, default=5)
    p.add_argument("--dpi",      type=int,   default=200)
    return p.parse_args()


def main():
    args = parse_args()
    n_tracks = len(DATA)

    fig, ax = plt.subplots(figsize=(args.width, args.height))

    y_positions = np.arange(n_tracks)[::-1]  # top to bottom

    for i, (label, pr, lo, hi, sr, cat) in enumerate(DATA):
        y   = y_positions[i]
        col = CATEGORY_COLORS[cat]

        # CI line
        ax.plot([lo, hi], [y, y], color=col, lw=2.0, solid_capstyle="round", zorder=2)

        # Pearson r dot (filled)
        ax.scatter(pr, y, color=col, s=80, zorder=3, edgecolors="white", lw=0.8)

        # Spearman r tick mark
        ax.scatter(sr, y, color=col, s=25, marker="D", zorder=3,
                   alpha=0.6, edgecolors="none")

        # Value annotation
        ax.text(lo - 0.01, y, f"{pr:.3f}", ha="right", va="center",
                fontsize=7.5, color=col)

    # Reference lines
    ax.axvline(0, color="black", lw=0.8, ls="-", zorder=1)
    for x in [-0.2, -0.4, -0.6]:
        ax.axvline(x, color="grey", lw=0.4, ls="--", alpha=0.5, zorder=0)

    # Y-axis labels
    ax.set_yticks(y_positions)
    ax.set_yticklabels([d[0] for d in DATA], fontsize=8.5)

    # Colour-code y tick labels by category
    for tick, (_, _, _, _, _, cat) in zip(ax.get_yticklabels(), DATA):
        tick.set_color(CATEGORY_COLORS[cat])

    # Axis formatting
    ax.set_xlabel("Pearson r (genome-wide, 1 Mb bins)\nNegative = B-compartment enriched",
                  fontsize=9)
    ax.set_xlim(-0.70, 0.05)
    ax.set_ylim(-0.7, n_tracks - 0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(left=False)

    # Category shading
    # Find y ranges per category
    cats_seen = {}
    for i, (_, _, _, _, _, cat) in enumerate(DATA):
        y = y_positions[i]
        if cat not in cats_seen:
            cats_seen[cat] = [y, y]
        else:
            cats_seen[cat][0] = min(cats_seen[cat][0], y)
            cats_seen[cat][1] = max(cats_seen[cat][1], y)

    for cat, (ymin, ymax) in cats_seen.items():
        ax.axhspan(ymin - 0.4, ymax + 0.4,
                   color=CATEGORY_COLORS[cat], alpha=0.06, zorder=0)

    # Legend
    legend_elements = [
        mpatches.Patch(color=CATEGORY_COLORS[cat], alpha=0.7, label=CATEGORY_LABELS[cat])
        for cat in ["condense-seq", "sensitivity"]
    ]
    # Add marker legend
    from matplotlib.lines import Line2D
    legend_elements += [
        Line2D([0], [0], marker="o", color="grey", lw=0,
               markersize=7, label="Pearson r (95% CI)"),
        Line2D([0], [0], marker="D", color="grey", lw=0,
               markersize=5, alpha=0.6, label="Spearman r"),
    ]
    ax.legend(handles=legend_elements, fontsize=7.5, loc="lower right",
              framealpha=0.9, edgecolor="grey")

    ax.set_title(
        f"Genome-wide correlation with MCF10A Hi-C compartment score\n"
        f"({args.binsize} bins, n ≈ 2,866, all chromosomes)",
        fontsize=9.5, pad=10
    )

    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"✓  Saved → {args.out}")


if __name__ == "__main__":
    main()
