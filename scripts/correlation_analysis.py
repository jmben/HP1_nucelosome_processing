#!/usr/bin/env python3
"""
correlation_analysis.py  —  Per-chromosome correlations + scatter plots
=============================================================================
Computes Pearson r between each track and the Hi-C compartment score
at 1 Mb resolution, broken down per chromosome, then produces:
  1. Per-chromosome heatmap
  2. Genome-wide forest plot with 95% CI
  3. Scatter plots (track vs compartment) for key tracks

Usage:
    python correlation_analysis.py \
        --compartment MCF10A_subcompartments_20000_T2T.bw \
        --out_dir     correlation_results/
=============================================================================
"""

import argparse
import numpy as np
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D
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
    # (label,                         filename,                                           category,       flip)
    # Compartment file convention: negative=B, positive=A
    # B-enriched tracks (high in B): flip=False → will show positive r (B covaries with B)
    # A-enriched tracks (high in A): flip=False → will show negative r (A anti-covaries with B)
    # Note: sensitivity files are high in A (open chromatin) → negative r expected
    # Note: H3K9me3, H3K27me3, HP1 are high in B → positive r expected
    # Note: H3K27ac, ATAC are high in A → negative r expected
    ("GM condense-seq\n(paper)",       "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T.bw", "Condense-seq", False),
    ("20U supernatant rep1",           "20U_supernatant_rep1_1kb_inverted.bw",             "Condense-seq", False),
    ("A_Inp rep1\nsensitivity",        "A_Inp_rep1_sensitivity.bw",                        "MNase",        False),
    ("Genomic rep1\nsensitivity",      "MCF10A_genomic_rep1_downsampled_sensitivity.bw",   "MNase",        False),
    ("Genomic rep2\nsensitivity",      "MCF10A_genomic_rep2_downsampled_sensitivity.bw",   "MNase",        False),
    ("ATAC-seq",                       "MCF10A_ATACseq_rep1_SRR12006157_coverage.bw",     "Chromatin",    False),
    ("H3K27ac",                        "MCF10A_H3K27ac_mean_RKPM.bw",                     "Chromatin",    False),
    ("H3K27me3",                       "MCF10A_H3K27me3_mean_RKPM.bw",                    "Chromatin",    False),
    ("H3K9me3",                        "MCF10A_H3K9me3_mean_RKPM.bw",                     "Chromatin",    False),
    ("HP1α ChIP",                      "WT_HP1a_signal.bw",                               "Chromatin",    False),  # fold enrichment from broadPeak col7
    ("HP1β ChIP",                      "WT_HP1b_signal.bw",                               "Chromatin",    False),  # fold enrichment from broadPeak col7
]

# Expected sign of Pearson r (compartment: negative=B, positive=A):
# Positive r = track enriched in B (condense-seq, H3K9me3, H3K27me3, HP1)
# Negative r = track enriched in A (ATAC, H3K27ac, sensitivity tracks)

CATEGORY_COLORS = {
    "Condense-seq": "#e67e22",
    "MNase":        "#2980b9",
    "Chromatin":    "#27ae60",
}

CHROMS_USE = [f"chr{i}" for i in list(range(1, 23))] + ["chrX"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def fetch_binned(bw_path: str, chrom: str, start: int, end: int,
                 binsize: int, flip: bool = False) -> np.ndarray:
    bw   = pyBigWig.open(bw_path)
    n    = (end - start) // binsize
    vals = np.array(
        bw.stats(chrom, start, start + n * binsize,
                 type="mean", nBins=n, exact=False),
        dtype=float)
    bw.close()
    if flip:
        vals = -vals
    finite = vals[np.isfinite(vals)]
    if len(finite) > 0:
        vals = np.where(np.isfinite(vals), vals, np.nanmean(finite))
    return vals


def pearson_ci(r: float, n: int, alpha: float = 0.05) -> tuple:
    from scipy.stats import norm
    if abs(r) >= 1.0 or n <= 3:
        return r, r
    z     = np.arctanh(r)
    se    = 1.0 / np.sqrt(n - 3)
    zc    = norm.ppf(1 - alpha / 2)
    return np.tanh(z - zc * se), np.tanh(z + zc * se)


def get_chrom_lengths(comp_file: str) -> dict:
    bw = pyBigWig.open(comp_file)
    cl = {c: l for c, l in bw.chroms().items()
          if c in CHROMS_USE}
    bw.close()
    return cl


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def compute_correlations(comp_file: str, binsize: int) -> dict:
    """
    Returns dict: track_label -> {
        'genome': (r, lo, hi, sr, n),
        'chroms': {chrom: (r, n)}
    }
    """
    chroms = get_chrom_lengths(comp_file)
    results = {}

    for label, fname, cat, flip in TRACKS:
        print(f"  {label.replace(chr(10), ' ')} ...")
        chrom_rs = {}
        all_comp, all_track = [], []

        # Check file exists
        if not Path(fname).exists():
            print(f"    WARNING: {fname} not found, skipping")
            results[label] = None
            continue

        bw_comp  = pyBigWig.open(comp_file)
        bw_track = pyBigWig.open(fname)

        for chrom, chrom_len in sorted(chroms.items()):
            n_bins = chrom_len // binsize
            if n_bins < 5:
                continue
            try:
                cv = np.array(bw_comp.stats(chrom, 0, n_bins*binsize,
                              type="mean", nBins=n_bins, exact=False), dtype=float)
                tv = np.array(bw_track.stats(chrom, 0, n_bins*binsize,
                              type="mean", nBins=n_bins, exact=False), dtype=float)
            except Exception:
                continue
            if flip:
                tv = -tv
            mask = np.isfinite(cv) & np.isfinite(tv)
            if mask.sum() < 5:
                continue
            pr, _ = pearsonr(cv[mask], tv[mask])
            chrom_rs[chrom] = (pr, int(mask.sum()))
            all_comp.extend(cv[mask].tolist())
            all_track.extend(tv[mask].tolist())

        bw_comp.close()
        bw_track.close()

        if len(all_comp) > 10:
            c, t = np.array(all_comp), np.array(all_track)
            pr, _  = pearsonr(c, t)
            sr, _  = spearmanr(c, t)
            lo, hi = pearson_ci(pr, len(c))
            results[label] = {
                "genome": (pr, lo, hi, sr, len(c)),
                "chroms": chrom_rs,
                "category": cat,
                "fname": fname,
                "all_comp": c,
                "all_track": t,
            }
        else:
            results[label] = None

    return results


# ---------------------------------------------------------------------------
# Plot 1: Per-chromosome heatmap
# ---------------------------------------------------------------------------
def plot_chrom_heatmap(results: dict, out_path: Path):
    labels   = [l for l, _, _, _ in TRACKS if results.get(l) is not None]
    chroms   = CHROMS_USE
    data_mat = np.full((len(labels), len(chroms)), np.nan)

    for i, label in enumerate(labels):
        res = results[label]
        if res is None:
            continue
        for j, chrom in enumerate(chroms):
            if chrom in res["chroms"]:
                data_mat[i, j] = res["chroms"][chrom][0]

    fig, ax = plt.subplots(figsize=(14, max(4, len(labels) * 0.55 + 1.5)))
    norm    = TwoSlopeNorm(vmin=-0.8, vcenter=0, vmax=0.8)
    im      = ax.imshow(data_mat, cmap="RdBu_r", norm=norm,
                        aspect="auto", interpolation="nearest")

    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=7)
    ax.set_yticks(range(len(labels)))
    clean_labels = [l.replace("\n", " ") for l in labels]
    ax.set_yticklabels(clean_labels, fontsize=8)

    # Colour y-tick labels by category
    for tick, label in zip(ax.get_yticklabels(), labels):
        cat = results[label]["category"]
        tick.set_color(CATEGORY_COLORS[cat])

    # Annotate cells with r value
    for i in range(len(labels)):
        for j in range(len(chroms)):
            v = data_mat[i, j]
            if np.isfinite(v):
                tc = "white" if abs(v) > 0.5 else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=5, color=tc)

    plt.colorbar(im, ax=ax, label="Pearson r", shrink=0.6, pad=0.01)
    ax.set_title("Per-chromosome Pearson r vs MCF10A Hi-C compartment score (1 Mb bins)",
                 fontsize=10, pad=10)

    # Category legend
    patches = [mpatches.Patch(color=c, label=k)
               for k, c in CATEGORY_COLORS.items()]
    ax.legend(handles=patches, fontsize=7, loc="upper right",
              bbox_to_anchor=(1.12, 1.0))

    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 2: Genome-wide forest plot — monochromatic, x from -1 to 1
# ---------------------------------------------------------------------------
def plot_forest(results: dict, out_path: Path):
    valid    = [(l, results[l]) for l, _, _, _ in TRACKS
                if results.get(l) is not None]
    n_tracks = len(valid)

    # Taller figure with more row spacing for larger text
    fig, ax  = plt.subplots(figsize=(12, max(8, n_tracks * 1.1 + 2.5)))
    y_pos    = np.arange(n_tracks)[::-1]

    DARK     = "#2c3e50"
    CI_COLOR = "#7f8c8d"

    ROW_H    = 1.0   # vertical spacing between rows

    for i, (label, res) in enumerate(valid):
        y                  = y_pos[i] * ROW_H
        pr, lo, hi, sr, n  = res["genome"]

        # CI line
        ax.plot([lo, hi], [y, y], color=CI_COLOR, lw=3.0,
                solid_capstyle="round", zorder=2)

        # Pearson r — filled circle
        ax.scatter(pr, y, color=DARK, s=130, zorder=4,
                   edgecolors="white", lw=1.2)

        # Spearman r — open diamond, offset vertically
        ax.scatter(sr, y - 0.32, color=DARK, s=60, marker="D",
                   zorder=3, facecolors="none", edgecolors=DARK, lw=1.2)

        # Pearson label — left of CI bar
        label_x = min(lo, pr) - 0.03
        ax.text(label_x, y + 0.07, f"r = {pr:.3f}",
                ha="right", va="bottom", fontsize=18, color=DARK,
                fontweight="bold")

        # Spearman label — below, aligned with diamond
        ax.text(label_x, y - 0.32, f"ρ = {sr:.3f}",
                ha="right", va="center", fontsize=16, color=CI_COLOR,
                style="italic")

        # 95% CI text — right of CI bar
        ax.text(max(hi, pr) + 0.02, y,
                f"[{lo:.2f}, {hi:.2f}]",
                ha="left", va="center", fontsize=16, color=CI_COLOR)

    # Reference lines
    ax.axvline(0, color="black", lw=1.2, zorder=1)
    for x in [-0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8]:
        ax.axvline(x, color="grey", lw=0.5, ls="--", alpha=0.35, zorder=0)

    # Y-axis track labels
    ax.set_yticks([y * ROW_H for y in y_pos])
    ax.set_yticklabels([l.replace("\n", " ") for l, _ in valid],
                       fontsize=20, color=DARK)
    ax.tick_params(left=False, labelsize=14)

    # X-axis
    ax.set_xlim(-1.0, 1.0)
    ax.set_xticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels([f"{x:.1f}" for x in
                        [-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0]],
                       fontsize=20)
    ax.set_xlabel("Pearson r vs MCF10A Hi-C compartment score (1 Mb bins)\n"
                  "Negative = B-compartment enriched  |  Positive = A-compartment enriched",
                  fontsize=13)

    ax.set_ylim(-0.9 * ROW_H, (n_tracks - 0.1) * ROW_H)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Alternating row shading
    for i in range(n_tracks):
        if i % 2 == 0:
            ax.axhspan((y_pos[i] - 0.5) * ROW_H, (y_pos[i] + 0.5) * ROW_H,
                       color="grey", alpha=0.06, zorder=0)

    # Legend
    legend_els = [
        Line2D([0], [0], marker="o", color=DARK, lw=0,
               markersize=11, markeredgecolor="white", label="Pearson r (95% CI)"),
        Line2D([0], [0], marker="D", color="none", lw=0,
               markersize=9, markeredgecolor=DARK, label="Spearman ρ"),
        Line2D([0], [0], color=CI_COLOR, lw=3.0, label="95% CI"),
    ]
    ax.legend(handles=legend_els, fontsize=18, loc="lower left",
              framealpha=0.9, edgecolor="grey")

    ax.set_title(
        "Genome-wide correlation with MCF10A Hi-C compartment\n"
        "(1 Mb bins, n ≈ 2,800, all autosomes + chrX)",
        fontsize=14, pad=14, color=DARK
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 3: Scatter plots with consistent y-axes per category
# ---------------------------------------------------------------------------
def plot_scatters(results: dict, out_path: Path):
    valid = [(l, results[l]) for l, _, _, _ in TRACKS
             if results.get(l) is not None]

    ncols = 4
    nrows = int(np.ceil(len(valid) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 3.8, nrows * 3.6))
    axes = axes.flatten()

    # Compute global percentile limits per category for consistent y-axes
    cat_ylims = {}
    for label, res in valid:
        cat = res["category"]
        t   = res["all_track"]
        lo  = np.percentile(t, 1)
        hi  = np.percentile(t, 99)
        if cat not in cat_ylims:
            cat_ylims[cat] = [lo, hi]
        else:
            cat_ylims[cat][0] = min(cat_ylims[cat][0], lo)
            cat_ylims[cat][1] = max(cat_ylims[cat][1], hi)

    # Fixed x-axis range — compartment score
    all_comp = np.concatenate([res["all_comp"] for _, res in valid])
    x_lo = np.percentile(all_comp, 1)
    x_hi = np.percentile(all_comp, 99)

    DARK = "#2c3e50"

    for i, (label, res) in enumerate(valid):
        ax  = axes[i]
        c   = res["all_comp"]
        t   = res["all_track"]
        pr, lo, hi, sr, n = res["genome"]
        cat = res["category"]

        # Clip to percentile range for display
        x_mask = (c >= x_lo) & (c <= x_hi)
        c_plot  = c[x_mask]
        t_plot  = t[x_mask]

        # Hexbin density
        ax.hexbin(c_plot, t_plot, gridsize=35, cmap="Blues",
                  mincnt=1, linewidths=0.1, alpha=0.85)

        # Regression line
        m, b_int = np.polyfit(c, t, 1)
        x_line   = np.linspace(x_lo, x_hi, 100)
        ax.plot(x_line, m * x_line + b_int, color=DARK, lw=1.8, zorder=5)

        # Consistent axes
        ax.set_xlim(x_lo, x_hi)
        ylo, yhi = cat_ylims[cat]
        pad = (yhi - ylo) * 0.08
        ax.set_ylim(ylo - pad, yhi + pad)

        ax.set_title(label.replace("\n", " "), fontsize=9.5,
                     color=CATEGORY_COLORS[cat], fontweight="bold", pad=4)
        ax.set_xlabel("Compartment score", fontsize=8)
        ax.set_ylabel("Track value", fontsize=8)
        ax.tick_params(labelsize=7.5)

        # Annotation box — Pearson on top line, Spearman below
        txt = f"r = {pr:.3f}  [{lo:.3f}, {hi:.3f}]\nρ = {sr:.3f}"
        ax.text(0.04, 0.97, txt,
                transform=ax.transAxes, fontsize=8, va="top", ha="left",
                color=DARK, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec="lightgrey", alpha=0.85))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Track vs MCF10A Hi-C compartment score (1 Mb bins, all chromosomes)",
                 fontsize=12, y=1.01, color=DARK)
    fig.tight_layout(h_pad=1.5, w_pad=1.0)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--compartment", default="MCF10A_subcompartments_20000_T2T.bw")
    p.add_argument("--binsize",     type=int, default=1_000_000)
    p.add_argument("--out_dir",     default="correlation_results/")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nBin size : {args.binsize//1_000_000} Mb")
    print(f"Output   : {out_dir}\n")

    print("Computing correlations ...")
    results = compute_correlations(args.compartment, args.binsize)

    # Print summary table
    print(f"\n{'Track':<30} {'Pearson r':>10} {'95% CI':>18} "
          f"{'Spearman r':>11} {'n':>6}")
    print("-" * 80)
    for label, _, cat, _ in TRACKS:
        res = results.get(label)
        if res is None:
            print(f"  {label.replace(chr(10),' '):<28}  (file not found)")
            continue
        pr, lo, hi, sr, n = res["genome"]
        print(f"  {label.replace(chr(10),' '):<28} {pr:>10.4f} "
              f" [{lo:>7.4f},{hi:>7.4f}]  {sr:>10.4f}  {n:>6,}")

    print("\nGenerating plots ...")
    plot_chrom_heatmap(results, out_dir / "per_chrom_heatmap.png")
    plot_forest(results,        out_dir / "genome_wide_forest.png")
    plot_scatters(results,      out_dir / "scatter_plots.png")

    print(f"\n✓  Done. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
