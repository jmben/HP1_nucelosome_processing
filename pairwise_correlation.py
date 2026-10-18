#!/usr/bin/env python3
"""
pairwise_correlation.py
=============================================================================
Computes pairwise Pearson r between all tracks at 1 Mb resolution,
then produces:
  1. Correlation heatmap (all pairs)
  2. Clustered heatmap (hierarchical clustering)
  3. Forest-style pairwise plot for key comparisons

Usage:
    python pairwise_correlation.py \
        --out_dir correlation_results/
=============================================================================
"""

import argparse
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D
    from matplotlib.colors import TwoSlopeNorm
except ImportError:
    raise ImportError("pip install matplotlib")

try:
    import pyBigWig
except ImportError:
    raise ImportError("pip install pyBigWig")


# ---------------------------------------------------------------------------
# Track definitions
# ---------------------------------------------------------------------------
TRACKS = [
    # (label,                         filename,                                           category)
    ("GM condense-seq\n(paper)",       "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T.bw", "Condense-seq"),
    ("20U supernatant rep1",           "20U_supernatant_rep1_1kb_inverted.bw",            "Condense-seq"),
    ("HiC",                            "MCF10A_subcompartments_20000_T2T.bw",             "Condense-seq"),
    ("Genomic rep2\nsensitivity",      "MCF10A_genomic_rep2_downsampled_sensitivity.bw",  "MNase"),
]

CATEGORY_COLORS = {
    "Condense-seq": "#E69F00",
    "MNase":        "#0072B2",
    "Chromatin":    "#009E73",
}

CHROMS_USE = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
BINSIZE    = 1_000_000


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def fetch_all_tracks(binsize: int) -> dict:
    """
    Returns dict: label -> np.array of binned values (genome-wide, NaN-masked).
    Only bins where ALL tracks have data are retained.
    """
    print("Loading tracks ...")
    raw = {}
    chrom_lengths = None

    for label, fname, cat in TRACKS:
        if not Path(fname).exists():
            print(f"  WARNING: {fname} not found, skipping {label!r}")
            continue
        bw     = pyBigWig.open(fname)
        chroms = {c: l for c, l in bw.chroms().items()
                  if c in CHROMS_USE}
        if chrom_lengths is None:
            chrom_lengths = chroms

        vals = []
        for chrom in sorted(chrom_lengths):
            n = chrom_lengths[chrom] // binsize
            if n < 3:
                continue
            try:
                v = np.array(bw.stats(chrom, 0, n * binsize,
                                      type="mean", nBins=n, exact=False),
                             dtype=float)
            except Exception:
                v = np.full(n, np.nan)
            vals.append(v)
        bw.close()

        raw[label] = np.concatenate(vals)
        finite = np.isfinite(raw[label]).sum()
        print(f"  {label.replace(chr(10),' '):<30} {finite:,} finite bins")

    # Common mask — only keep bins where all tracks have data
    labels = list(raw.keys())
    stack  = np.stack([raw[l] for l in labels], axis=0)
    mask   = np.all(np.isfinite(stack), axis=0)
    print(f"\nBins with complete data across all tracks: {mask.sum():,}")

    return {l: stack[i][mask] for i, l in enumerate(labels)}


# ---------------------------------------------------------------------------
# Pairwise correlation matrix
# ---------------------------------------------------------------------------
def compute_pairwise(data: dict) -> tuple:
    """
    Returns (labels, pearson_r_matrix, spearman_r_matrix, ci_lo_matrix, ci_hi_matrix)
    """
    from scipy.stats import norm

    labels = list(data.keys())
    n      = len(labels)
    pr_mat = np.zeros((n, n))
    sr_mat = np.zeros((n, n))
    lo_mat = np.zeros((n, n))
    hi_mat = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                pr_mat[i, j] = sr_mat[i, j] = 1.0
                lo_mat[i, j] = hi_mat[i, j] = 1.0
                continue
            pr, _ = pearsonr(data[labels[i]], data[labels[j]])
            sr, _ = spearmanr(data[labels[i]], data[labels[j]])
            pr_mat[i, j] = pr
            sr_mat[i, j] = sr
            # Fisher z CI (95%)
            nobs = len(data[labels[i]])
            z    = np.arctanh(pr)
            se   = 1 / np.sqrt(nobs - 3)
            zc   = norm.ppf(0.975)
            lo_mat[i, j] = np.tanh(z - zc * se)
            hi_mat[i, j] = np.tanh(z + zc * se)

    return labels, pr_mat, sr_mat, lo_mat, hi_mat


# ---------------------------------------------------------------------------
# Plot 1: Clustered correlation heatmap
# ---------------------------------------------------------------------------
def plot_heatmap(labels: list, pr_mat: np.ndarray, out_path: Path):
    """Hierarchically clustered Pearson r heatmap."""

    # Cluster using 1-r as distance
    dist_condensed = squareform(1 - pr_mat, checks=False)
    Z     = linkage(dist_condensed, method="average")
    order = leaves_list(Z)

    mat_ordered  = pr_mat[np.ix_(order, order)]
    labels_clean = [labels[i].replace("\n", " ") for i in order]
    cats_ordered = [TRACKS[[l for l, _, _ in TRACKS].index(labels[i])][2]
                    for i in order]

    fig, ax = plt.subplots(figsize=(13, 11))
    norm    = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im      = ax.imshow(mat_ordered, cmap="RdBu_r", norm=norm,
                        aspect="auto", interpolation="nearest")

    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(labels_clean, rotation=45, ha="right", fontsize=12)
    ax.set_yticklabels(labels_clean, fontsize=12)

    # Colour tick labels by category
    for tick, cat in zip(ax.get_xticklabels(), cats_ordered):
        tick.set_color(CATEGORY_COLORS[cat])
    for tick, cat in zip(ax.get_yticklabels(), cats_ordered):
        tick.set_color(CATEGORY_COLORS[cat])

    # Annotate cells
    for i in range(len(order)):
        for j in range(len(order)):
            v  = mat_ordered[i, j]
            tc = "white" if abs(v) > 0.6 else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=14, color=tc)

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("Pearson r", fontsize=14)
    cbar.ax.tick_params(labelsize=11)

    # Category legend
    patches = [mpatches.Patch(color=c, label=k)
               for k, c in CATEGORY_COLORS.items()]
    ax.legend(handles=patches, fontsize=11, loc="upper right",
              bbox_to_anchor=(1.18, 1.0), framealpha=0.9)

    ax.set_title("Pairwise Pearson r between tracks\n"
                 "(1 Mb bins, genome-wide, hierarchically clustered)",
                 fontsize=18, pad=12)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 2: Forest plot of key pairwise comparisons
# ---------------------------------------------------------------------------
def plot_pairwise_forest(labels: list, pr_mat: np.ndarray,
                         lo_mat: np.ndarray, hi_mat: np.ndarray,
                         sr_mat: np.ndarray, out_path: Path):
    """
    Forest plot showing r for each track vs GM condense-seq (paper reference)
    and vs Genomic rep2 sensitivity, to show where each track sits relative
    to the two key references.
    """
    # Reference tracks to compare against
    refs = ["GM condense-seq\n(paper)", "Genomic rep2\nsensitivity"]
    ref_colors = {"GM condense-seq\n(paper)":  "#E69F00",
                  "Genomic rep2\nsensitivity":  "#0072B2"}

    # All other tracks (excluding the references themselves)
    others = [l for l in labels if l not in refs]

    fig, axes = plt.subplots(1, len(refs),
                             figsize=(8 * len(refs), max(8, len(others) * 0.9 + 3)),
                             sharey=True)

    DARK     = "#2c3e50"
    CI_COLOR = "#aaaaaa"

    for ax, ref in zip(axes, refs):
        ref_idx  = labels.index(ref)
        col      = ref_colors[ref]
        y_pos    = np.arange(len(others))[::-1]

        for i, other in enumerate(others):
            other_idx = labels.index(other)
            pr = pr_mat[ref_idx, other_idx]
            lo = lo_mat[ref_idx, other_idx]
            hi = hi_mat[ref_idx, other_idx]
            sr = sr_mat[ref_idx, other_idx]
            y  = y_pos[i]

            # CI line
            ax.plot([lo, hi], [y, y], color=CI_COLOR, lw=2.5,
                    solid_capstyle="round", zorder=2)
            # Pearson dot
            ax.scatter(pr, y, color=col, s=120, zorder=4,
                       edgecolors="white", lw=1.2)
            # Spearman diamond
            ax.scatter(sr, y - 0.28, color=col, s=55, marker="D",
                       zorder=3, facecolors="none", edgecolors=col, lw=1.2)

            # r value label
            label_x = min(lo, pr) - 0.03 if pr < 0 else max(hi, pr) + 0.03
            ha      = "right" if pr < 0 else "left"
            ax.text(label_x, y + 0.07, f"r = {pr:.3f}",
                    ha=ha, va="bottom", fontsize=12, color=DARK,
                    fontweight="bold")
            ax.text(label_x, y - 0.28, f"ρ = {sr:.3f}",
                    ha=ha, va="center", fontsize=10, color=CI_COLOR,
                    style="italic")

        ax.axvline(0, color="black", lw=1.0, zorder=1)
        for x in [-0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8]:
            ax.axvline(x, color="grey", lw=0.4, ls="--", alpha=0.3, zorder=0)

        # Alternating shading
        for i in range(len(others)):
            if i % 2 == 0:
                ax.axhspan(y_pos[i] - 0.5, y_pos[i] + 0.5,
                           color="grey", alpha=0.06, zorder=0)

        ax.set_yticks(y_pos)
        ax.set_yticklabels([l.replace("\n", " ") for l in others],
                           fontsize=13, color=DARK)
        ax.tick_params(left=False)
        ax.set_xlim(-1.0, 1.0)
        ax.set_xticks([-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8])
        ax.set_xticklabels([f"{x:.1f}" for x in
                            [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]],
                           fontsize=12)
        ax.set_xlabel("Pearson r", fontsize=13)
        ax.set_ylim(-0.9, len(others) - 0.1)
        ax.set_title(f"Correlation with\n{ref.replace(chr(10), ' ')}",
                     fontsize=14, color=col, fontweight="bold", pad=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        legend_els = [
            Line2D([0], [0], marker="o", color=col, lw=0,
                   markersize=10, markeredgecolor="white",
                   label="Pearson r (95% CI)"),
            Line2D([0], [0], marker="D", color="none", lw=0,
                   markersize=8, markeredgecolor=col, label="Spearman ρ"),
        ]
        ax.legend(handles=legend_els, fontsize=11, loc="lower right",
                  framealpha=0.9, edgecolor="grey")

    fig.suptitle("Pairwise track correlations (1 Mb bins, genome-wide)",
                 fontsize=15, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------------------
def print_table(labels: list, pr_mat: np.ndarray,
                sr_mat: np.ndarray, out_path: Path):
    clean = [l.replace("\n", " ") for l in labels]
    n     = len(labels)

    lines = ["Pairwise Pearson r matrix (lower triangle)\n",
             f"{'Track':<32}" + "".join(f"  {c[:8]:>8}" for c in clean),
             "-" * (32 + 10 * n)]

    for i in range(n):
        row = f"{clean[i]:<32}"
        for j in range(n):
            if j > i:
                row += f"  {'':>8}"
            else:
                row += f"  {pr_mat[i,j]:>8.4f}"
        lines.append(row)

    text = "\n".join(lines)
    print(text)
    with open(out_path, "w") as f:
        f.write(text + "\n")
    print(f"\n  Table saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--binsize",  type=int,   default=1_000_000)
    p.add_argument("--out_dir",  default="correlation_results/")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load all tracks
    data = fetch_all_tracks(args.binsize)

    if len(data) < 2:
        print("Need at least 2 tracks — check file paths.")
        return

    # Compute pairwise correlations
    print("\nComputing pairwise correlations ...")
    labels, pr_mat, sr_mat, lo_mat, hi_mat = compute_pairwise(data)

    # Print and save table
    print_table(labels, pr_mat, sr_mat, out_dir / "pairwise_correlation_table.txt")

    # Plots
    print("\nGenerating plots ...")
    plot_heatmap(labels, pr_mat,
                 out_dir / "pairwise_heatmap.png")
    plot_pairwise_forest(labels, pr_mat, lo_mat, hi_mat, sr_mat,
                         out_dir / "pairwise_forest.png")

    print(f"\n✓  Done. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
