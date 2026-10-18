#!/usr/bin/env python3
"""
subcompartment_analysis.py  —  Track values stratified by subcompartment
=============================================================================
Assigns each genomic bin a subcompartment label (A0-A3, B0-B3) based on
majority overlap with the subcompartment bedgraph, then plots track value
distributions per subcompartment for each track.

Usage:
    python subcompartment_analysis.py \
        --bedgraph MCF10A_subcompartments_20000_T2T.bg \
        --binsize  1000000 \
        --out_dir  subcompartment_results/
=============================================================================
"""

import argparse
import numpy as np
from pathlib import Path
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
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
    ("GM condense-seq (paper)",   "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T.bw", False),
    ("20U supernatant rep1",      "20U_supernatant_rep1_1kb_inverted.bw",            False),
    ("A_Inp rep1 sensitivity",    "A_Inp_rep1_sensitivity.bw",                       False),
    ("Genomic rep1 sensitivity",  "MCF10A_genomic_rep1_downsampled_sensitivity.bw",  False),
    ("Genomic rep2 sensitivity",  "MCF10A_genomic_rep2_downsampled_sensitivity.bw",  False),
    ("ATAC-seq",                  "MCF10A_ATACseq_rep1_SRR12006157_coverage.bw",    False),
    ("H3K27ac",                   "MCF10A_H3K27ac_mean_RKPM.bw",                    False),
    ("H3K27me3",                  "MCF10A_H3K27me3_mean_RKPM.bw",                   False),
    ("H3K9me3",                   "MCF10A_H3K9me3_mean_RKPM.bw",                    False),
    ("HP1α ChIP",                 "WT_HP1a_signal.bw",                              False),
    ("HP1β ChIP",                 "WT_HP1b_signal.bw",                              False),
]

# Subcompartment order and colours (A = active, B = inactive)
SUBCOMPS   = ["A3", "A2", "A1", "A0", "B0", "B1", "B2", "B3"]
SUBCOMP_COLORS = {
    "A0": "#c0392b", "A1": "#e74c3c", "A2": "#f39c12", "A3": "#f1c40f",
    "B0": "#1a5276", "B1": "#2980b9", "B2": "#7fb3d3", "B3": "#aed6f1",
}

CHROMS_USE = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


# ---------------------------------------------------------------------------
# Load subcompartment bedgraph into per-bin labels
# ---------------------------------------------------------------------------
def load_subcompartments(bg_path: str, binsize: int) -> dict:
    """
    Returns dict: (chrom, bin_start) -> dominant_subcompartment_label
    Dominant = subcompartment with most bp overlap in that bin.
    """
    print(f"Loading subcompartments from {bg_path} ...")
    # bin_overlap[(chrom, bin_start)][label] = total bp overlap
    bin_overlap = defaultdict(lambda: defaultdict(int))

    with open(bg_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom, start, end, label = parts[0], int(parts[1]), int(parts[2]), parts[3]
            if chrom not in CHROMS_USE:
                continue

            # Find all bins this interval overlaps
            bin_start_first = (start // binsize) * binsize
            b = bin_start_first
            while b < end:
                overlap = min(end, b + binsize) - max(start, b)
                if overlap > 0:
                    bin_overlap[(chrom, b)][label] += overlap
                b += binsize

    # Assign dominant label per bin
    bin_labels = {}
    for (chrom, b), label_counts in bin_overlap.items():
        bin_labels[(chrom, b)] = max(label_counts, key=label_counts.get)

    label_counts_total = defaultdict(int)
    for label in bin_labels.values():
        label_counts_total[label] += 1
    print(f"  Total bins assigned: {len(bin_labels):,}")
    for sc in SUBCOMPS:
        print(f"    {sc}: {label_counts_total.get(sc, 0):,} bins")

    return bin_labels


# ---------------------------------------------------------------------------
# Fetch track values at binned positions
# ---------------------------------------------------------------------------
def fetch_track_by_subcomp(bw_path: str, bin_labels: dict,
                            binsize: int, flip: bool) -> dict:
    """
    Returns dict: subcompartment_label -> list of track values
    """
    bw     = pyBigWig.open(bw_path)
    chroms = set(c for c, _ in bin_labels.keys())
    values = defaultdict(list)

    for chrom in sorted(chroms):
        if chrom not in (bw.chroms() or {}):
            continue
        # Get all bins for this chrom
        chrom_bins = sorted(b for c, b in bin_labels if c == chrom)
        for b in chrom_bins:
            label = bin_labels[(chrom, b)]
            try:
                v = bw.stats(chrom, b, b + binsize, type="mean", nBins=1, exact=False)
                if v and v[0] is not None and np.isfinite(v[0]):
                    val = -v[0] if flip else v[0]
                    values[label].append(val)
            except Exception:
                continue

    bw.close()
    return {k: np.array(v) for k, v in values.items()}


# ---------------------------------------------------------------------------
# Plot 1: Violin plots per track
# ---------------------------------------------------------------------------
def plot_violins(all_data: dict, out_path: Path):
    """
    all_data: {track_label: {subcomp: np.array of values}}
    One row per track, one violin per subcompartment.
    """
    valid_tracks = [t for t in all_data if all_data[t]]
    n_tracks     = len(valid_tracks)
    n_sc         = len(SUBCOMPS)

    fig, axes = plt.subplots(n_tracks, 1,
                             figsize=(12, n_tracks * 2.2),
                             sharex=False)
    if n_tracks == 1:
        axes = [axes]

    for ax, track_label in zip(axes, valid_tracks):
        data_by_sc = all_data[track_label]
        plot_data  = [data_by_sc.get(sc, np.array([])) for sc in SUBCOMPS]
        plot_data  = [d if len(d) > 1 else np.array([np.nan]) for d in plot_data]

        # Z-score normalise across all subcompartments for comparability
        all_vals = np.concatenate([d for d in plot_data if np.isfinite(d).any()])
        mu, sd   = np.nanmean(all_vals), np.nanstd(all_vals)
        if sd > 0:
            plot_data = [(d - mu) / sd for d in plot_data]

        positions = np.arange(1, n_sc + 1)
        vp = ax.violinplot(
            [d[np.isfinite(d)] for d in plot_data],
            positions=positions,
            showmedians=True,
            showextrema=False,
            widths=0.7,
        )

        # Colour each violin by subcompartment
        for body, sc in zip(vp["bodies"], SUBCOMPS):
            body.set_facecolor(SUBCOMP_COLORS[sc])
            body.set_alpha(0.75)
            body.set_edgecolor("grey")
            body.set_linewidth(0.5)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.5)

        ax.set_ylabel(track_label.replace("\n", " "), fontsize=8, labelpad=4)
        ax.set_xticks(positions)
        ax.set_xticklabels(SUBCOMPS, fontsize=8)
        ax.axvline(4.5, color="black", lw=0.8, ls="--", alpha=0.5)  # A|B boundary
        ax.axhline(0, color="grey", lw=0.4, ls=":")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=7)

    axes[-1].set_xlabel("Subcompartment", fontsize=10)
    fig.suptitle("Track values by subcompartment (z-scored, 1 Mb bins)",
                 fontsize=11, y=1.002)
    fig.tight_layout(h_pad=0.5)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 2: Mean per subcompartment heatmap
# ---------------------------------------------------------------------------
def plot_subcomp_heatmap(all_data: dict, out_path: Path):
    valid_tracks = [t for t in all_data if all_data[t]]
    n_tracks     = len(valid_tracks)

    # Build matrix: tracks × subcompartments (z-scored means)
    mat = np.full((n_tracks, len(SUBCOMPS)), np.nan)
    for i, track_label in enumerate(valid_tracks):
        data_by_sc = all_data[track_label]
        vals = [data_by_sc.get(sc, np.array([])) for sc in SUBCOMPS]
        means = np.array([np.nanmean(v) if len(v) > 0 else np.nan for v in vals])
        # Z-score across subcompartments
        mu, sd = np.nanmean(means), np.nanstd(means)
        if sd > 0:
            mat[i] = (means - mu) / sd

    from scipy.stats import spearmanr
    sc_rank = {sc: i for i, sc in enumerate(SUBCOMPS)}

    fig, ax = plt.subplots(figsize=(10, max(4, n_tracks * 0.55 + 1.5)))
    vmax = np.nanpercentile(np.abs(mat), 98)
    im   = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                     aspect="auto", interpolation="nearest")

    # Build y-tick labels with Spearman r appended
    y_labels = []
    for track_label in valid_tracks:
        data_by_sc = all_data[track_label]
        all_vals, all_ranks = [], []
        for sc in SUBCOMPS:
            v = data_by_sc.get(sc, np.array([]))
            if len(v) > 0:
                all_vals.extend(v.tolist())
                all_ranks.extend([sc_rank[sc]] * len(v))
        if len(all_vals) > 10:
            sr, sp = spearmanr(all_ranks, all_vals)
            sig = "***" if sp < 0.001 else ("**" if sp < 0.01 else "*" if sp < 0.05 else "")
            y_labels.append(f"{track_label}  (ρ={sr:.3f}{sig})")
        else:
            y_labels.append(track_label)

    ax.set_xticks(range(len(SUBCOMPS)))
    ax.set_xticklabels(SUBCOMPS, fontsize=10)
    ax.set_yticks(range(n_tracks))
    ax.set_yticklabels(y_labels, fontsize=8.5)

    # Annotate with z-score
    for i in range(n_tracks):
        for j in range(len(SUBCOMPS)):
            v  = mat[i, j]
            if np.isfinite(v):
                tc = "white" if abs(v) > vmax * 0.6 else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7.5, color=tc)

    # A|B boundary line
    ax.axvline(3.5, color="black", lw=2.0)
    ax.text(1.5, -0.7, "A compartment", ha="center", fontsize=9,
            transform=ax.get_xaxis_transform())
    ax.text(5.5, -0.7, "B compartment", ha="center", fontsize=9,
            transform=ax.get_xaxis_transform())

    plt.colorbar(im, ax=ax, label="Z-score (mean track value)", shrink=0.6)
    ax.set_title("Mean track value per subcompartment (z-scored across subcompartments)",
                 fontsize=10, pad=12)

    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 3: Line plot — mean per subcompartment ordered A0→B3
# ---------------------------------------------------------------------------
def plot_subcomp_profiles(all_data: dict, out_path: Path):
    valid_tracks = [t for t in all_data if all_data[t]]

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(SUBCOMPS))

    # Group by assay type for colour
    track_colors = {}
    for label, fname, flip in TRACKS:
        if "condense" in fname.lower() or "supernatant" in fname.lower():
            track_colors[label] = "#e67e22"
        elif "sensitivity" in fname.lower():
            track_colors[label] = "#2980b9"
        else:
            track_colors[label] = "#27ae60"

    for track_label in valid_tracks:
        data_by_sc = all_data[track_label]
        means = np.array([np.nanmean(data_by_sc.get(sc, [np.nan]))
                          for sc in SUBCOMPS])
        # Z-score
        mu, sd = np.nanmean(means), np.nanstd(means)
        if sd > 0:
            means = (means - mu) / sd

        col = track_colors.get(track_label, "#555555")
        ax.plot(x, means, marker="o", lw=1.8, ms=5,
                color=col, alpha=0.8, label=track_label)

    ax.axvline(3.5, color="black", lw=1.0, ls="--", alpha=0.6)
    ax.axhline(0, color="grey", lw=0.5, ls=":")
    ax.set_xticks(x)
    ax.set_xticklabels(SUBCOMPS, fontsize=11)
    ax.set_ylabel("Z-scored mean track value", fontsize=10)
    ax.set_xlabel("Subcompartment (A3 = most active → B3 = most inactive)",
                  fontsize=9)
    ax.set_title("Track profiles across subcompartments", fontsize=11)

    # Legend outside
    ax.legend(fontsize=7.5, bbox_to_anchor=(1.01, 1), loc="upper left",
              framealpha=0.9, edgecolor="grey")

    # Category patches
    patches = [
        mpatches.Patch(color="#e67e22", label="Condense-seq"),
        mpatches.Patch(color="#2980b9", label="MNase sensitivity"),
        mpatches.Patch(color="#27ae60", label="Chromatin marks"),
    ]
    ax.text(1.5, ax.get_ylim()[0] * 0.9 if ax.get_ylim()[0] < 0 else -0.1,
            "← A compartment →", ha="center", fontsize=8, color="grey")
    ax.text(5.5, ax.get_ylim()[0] * 0.9 if ax.get_ylim()[0] < 0 else -0.1,
            "← B compartment →", ha="center", fontsize=8, color="grey")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bedgraph", default="MCF10A_subcompartments_20000_T2T.bg")
    p.add_argument("--binsize",  type=int, default=1_000_000)
    p.add_argument("--out_dir",  default="subcompartment_results/")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load subcompartment labels
    bin_labels = load_subcompartments(args.bedgraph, args.binsize)

    # Fetch track values per subcompartment
    all_data = {}
    for label, fname, flip in TRACKS:
        if not Path(fname).exists():
            print(f"  Skipping {label} — {fname} not found")
            all_data[label] = {}
            continue
        print(f"  Fetching: {label} ...")
        all_data[label] = fetch_track_by_subcomp(fname, bin_labels,
                                                  args.binsize, flip)

    # Print summary table with statistics
    from scipy.stats import spearmanr, kruskal

    # Subcompartment rank (A3=most active=0 ... B3=most inactive=7)
    sc_rank = {sc: i for i, sc in enumerate(SUBCOMPS)}

    print(f"\n{'Track':<30}", end="")
    for sc in SUBCOMPS:
        print(f"  {sc:>6}", end="")
    print(f"  {'Spearman r':>11}  {'KW p':>10}")
    print("-" * (30 + 9 * len(SUBCOMPS) + 25))

    stats_results = {}
    for label, _, _ in TRACKS:
        data = all_data.get(label, {})
        if not data:
            continue
        print(f"{label:<30}", end="")
        vals_all = np.concatenate([v for v in data.values() if len(v) > 0])
        mu, sd   = np.nanmean(vals_all), np.nanstd(vals_all)

        # Z-score means per subcompartment
        z_means = []
        for sc in SUBCOMPS:
            v = data.get(sc, np.array([]))
            m = np.nanmean(v) if len(v) > 0 else np.nan
            z = (m - mu) / sd if sd > 0 and np.isfinite(m) else np.nan
            z_means.append(z)
            print(f"  {z:>6.2f}", end="")

        # Spearman r between subcompartment rank and all individual bin values
        all_vals, all_ranks = [], []
        for sc in SUBCOMPS:
            v = data.get(sc, np.array([]))
            if len(v) > 0:
                all_vals.extend(v.tolist())
                all_ranks.extend([sc_rank[sc]] * len(v))

        sr, sp = spearmanr(all_ranks, all_vals)
        sig = "***" if sp < 0.001 else ("**" if sp < 0.01 else ("*" if sp < 0.05 else "ns"))
        print(f"  {sr:>10.4f}{sig}", end="")

        # Kruskal-Wallis test across subcompartments
        groups = [data.get(sc, np.array([])) for sc in SUBCOMPS]
        groups = [g for g in groups if len(g) > 1]
        if len(groups) >= 2:
            kw_stat, kw_p = kruskal(*groups)
            kw_sig = "***" if kw_p < 0.001 else ("**" if kw_p < 0.01 else ("*" if kw_p < 0.05 else "ns"))
            print(f"  {kw_p:>10.2e}{kw_sig}")
        else:
            print(f"  {'NA':>10}")

        stats_results[label] = {"spearman_r": sr, "spearman_p": sp,
                                 "z_means": z_means}

    # Save stats to text file
    stats_path = out_dir / "subcomp_statistics.txt"
    with open(stats_path, "w") as f:
        f.write("Subcompartment analysis statistics\n")
        f.write("=" * 60 + "\n\n")
        f.write("Spearman r: correlation between subcompartment rank\n")
        f.write("(A3=0 most active ... B3=7 most inactive) and\n")
        f.write("individual 1 Mb bin track values.\n")
        f.write("Positive r = track increases with heterochromatin.\n\n")
        f.write(f"{'Track':<30} {'Spearman r':>11} {'p-value':>12} {'sig':>5}\n")
        f.write("-" * 62 + "\n")
        for label in stats_results:
            sr  = stats_results[label]["spearman_r"]
            sp  = stats_results[label]["spearman_p"]
            sig = "***" if sp < 0.001 else ("**" if sp < 0.01 else ("*" if sp < 0.05 else "ns"))
            f.write(f"{label:<30} {sr:>11.4f} {sp:>12.2e} {sig:>5}\n")
    print(f"\n  Stats saved: {stats_path}")

    # Generate plots
    print("\nGenerating plots ...")
    plot_violins(all_data,           out_dir / "subcomp_violins.png")
    plot_subcomp_heatmap(all_data,   out_dir / "subcomp_heatmap.png")
    plot_subcomp_profiles(all_data,  out_dir / "subcomp_profiles.png")

    print(f"\n✓  Done. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
