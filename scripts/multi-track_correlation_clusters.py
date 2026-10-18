#!/usr/bin/env python3
"""
multi-track_correlation_clusters.py  (v3)
=========================================
Changes from v2:
  - TSS list replaced by cluster BED (deepTools format with cluster_N in col 12)
  - Per-cluster TSS profile plots (one panel per cluster)
  - Cluster-stratified correlation matrices
  - Genome-wide compartment-stratified correlation retained
"""

import os
import argparse
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# =============================================================================
# Config
# =============================================================================

DEFAULT_BIGWIG_DIR  = "."
DEFAULT_CLUSTER_BED = "HP1_HPTM_2kb_sen_matrix_og_sort_k7.bed"
DEFAULT_OUTDIR      = "./correlation_output_v3"
DEFAULT_BIN_SIZE    = 10_000
DEFAULT_TSS_WINDOW  = 1_000   # ±1kb (BED entries are already 2kb wide, centred on TSS)
DEFAULT_TSS_BINS    = 80
DEFAULT_COV_PCT     = 20
COMPARTMENT_TRACK   = "MCF10A_subcompartments_20000_T2T.bw"
COV_FILTER_TRACK    = "20U_supernatant_rep1_1kb_inverted.bw"

SELECTED_TRACKS = [
    "MCF10A_subcompartments_20000_T2T.bw",
    "A_Inp_rep1_distribution.bw",
    "A_Inp_rep1_sensitivity.bw",
    "MCF10A_ATACseq_rep1_SRR12006157_coverage.bw",
    "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T.bw",
    "20U_supernatant_rep1_1kb_inverted.bw",
    "20U_supernatant_rep2_1kb_inverted.bw",
]

SHORT_LABELS = {
    "MCF10A_subcompartments_20000_T2T":              "Hi-C compartments",
    "A_Inp_rep1_distribution":                       "Nuc distribution",
    "A_Inp_rep1_sensitivity":                        "Nuc sensitivity",
    "MCF10A_ATACseq_rep1_SRR12006157_coverage":      "ATAC-seq",
    "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T":   "GM NCP condensability",
    "20U_supernatant_rep1_1kb_inverted":             "Cond 20U supernatant rep1",
    "20U_supernatant_rep2_1kb_inverted":             "Cond 20U supernatant rep2",
}

TRACK_ORDER = [
    "MCF10A_subcompartments_20000_T2T",
    "A_Inp_rep1_distribution",
    "A_Inp_rep1_sensitivity",
    "MCF10A_ATACseq_rep1_SRR12006157_coverage",
    "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T",
    "20U_supernatant_rep1_1kb_inverted",
    "20U_supernatant_rep2_1kb_inverted",
]

TRACK_COLORS = {
    "MCF10A_subcompartments_20000_T2T":              "#7B2D8B",
    "A_Inp_rep1_distribution":                       "#2166AC",
    "A_Inp_rep1_sensitivity":                        "#6BAED6",
    "MCF10A_ATACseq_rep1_SRR12006157_coverage":      "#1A9850",
    "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T":   "#D6604D",
    "20U_supernatant_rep1_1kb_inverted":             "#B2182B",
    "20U_supernatant_rep2_1kb_inverted":             "#D6604D",
}

# Colours for cluster lines in per-cluster plots
CLUSTER_COLORS = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF", "#999999",
]


# =============================================================================
# I/O helpers
# =============================================================================

def load_chrom_sizes_from_bw(bw_path):
    bw = pyBigWig.open(bw_path)
    sizes = dict(bw.chroms())
    bw.close()
    return sizes


def fetch_bw_bins(bw_path, chrom, chrom_size, bin_size):
    n_bins = chrom_size // bin_size
    if n_bins == 0:
        return np.array([])
    try:
        bw = pyBigWig.open(bw_path)
        vals = bw.stats(chrom, 0, n_bins * bin_size, type="mean", nBins=n_bins)
        bw.close()
        return np.array(vals, dtype=np.float64)
    except Exception:
        return np.full(n_bins, np.nan)


# =============================================================================
# Cluster BED loader
# =============================================================================

def load_cluster_bed(bed_path, chrom_sizes):
    """
    Parses a deepTools-style BED with cluster labels in the last column.
    Centres each entry on (start+end)//2 and records strand from col 5.
    Returns:
        tss_all   : list of (chrom, centre, strand, cluster_label)
        clusters  : sorted list of unique cluster labels
    """
    tss_all  = []
    skipped  = 0

    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom  = parts[0]
            start  = int(parts[1])
            end    = int(parts[2])
            strand = parts[5] if len(parts) >= 6 else "+"
            # cluster label is the last column
            cluster = parts[-1].strip()

            centre = (start + end) // 2

            # validate against chrom sizes
            if chrom not in chrom_sizes:
                skipped += 1
                continue

            tss_all.append((chrom, centre, strand, cluster))

    clusters = sorted(set(t[3] for t in tss_all))
    print(f"  Loaded {len(tss_all):,} entries ({skipped} skipped) "
          f"across {len(clusters)} clusters: {clusters}", flush=True)
    return tss_all, clusters


# =============================================================================
# Genome-wide matrix
# =============================================================================

def build_genome_wide_matrix(bw_paths, chrom_sizes, bin_size):
    chroms = [c for c, s in chrom_sizes.items()
              if s > bin_size * 10
              and not any(x in c for x in ["_", "chrM", "chrEBV", "chrUn"])]
    print(f"  Using {len(chroms)} chromosomes", flush=True)

    all_data = {os.path.splitext(os.path.basename(p))[0]: [] for p in bw_paths}
    for chrom in chroms:
        chrom_size = chrom_sizes[chrom]
        print(f"    {chrom}", flush=True)
        for bw_path in bw_paths:
            stem = os.path.splitext(os.path.basename(bw_path))[0]
            vals = fetch_bw_bins(bw_path, chrom, chrom_size, bin_size)
            all_data[stem].extend(vals.tolist())

    return pd.DataFrame(all_data)


def apply_coverage_filter(df, coverage_col, percentile=20):
    if coverage_col not in df.columns:
        print(f"  WARNING: coverage filter col '{coverage_col}' not found, skipping")
        return df
    thresh = np.nanpercentile(df[coverage_col].dropna(), percentile)
    filtered = df[df[coverage_col] >= thresh]
    print(f"  Coverage filter (>={percentile}th pct, thresh={thresh:.3f}): "
          f"{len(filtered):,} / {len(df):,} bins kept", flush=True)
    return filtered


def stratify_by_compartment(df, compartment_col):
    if compartment_col not in df.columns:
        print(f"  WARNING: '{compartment_col}' not found")
        return df.copy(), df.copy()
    df_a = df[df[compartment_col] > 0].copy()
    df_b = df[df[compartment_col] < 0].copy()
    print(f"  A: {len(df_a):,} bins | B: {len(df_b):,} bins", flush=True)
    return df_a, df_b


# =============================================================================
# Correlation + heatmap plotting
# =============================================================================

def compute_correlation_matrix(df, method="spearman"):
    df_clean = df.dropna()
    print(f"  {len(df_clean):,} complete rows used", flush=True)
    return df_clean.corr(method=method)


def _draw_heatmap(ax, corr_df, cols, labels, title, norm):
    data = corr_df.reindex(index=cols, columns=cols).values
    n    = len(cols)
    im   = ax.imshow(data, cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title(title, fontsize=9)
    for i in range(n):
        for j in range(n):
            val = data[i, j]
            if not np.isnan(val):
                color = "white" if abs(val) > 0.6 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=5.5, color=color)
    return im


def plot_correlation_heatmap(corr_df, title, out_path, track_order=None):
    stems  = track_order or list(corr_df.columns)
    cols   = [s for s in stems if s in corr_df.columns]
    labels = [SHORT_LABELS.get(c, c) for c in cols]
    n      = len(cols)
    norm   = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    fig, ax = plt.subplots(figsize=(max(7, n * 0.65), max(6, n * 0.65)))
    im = _draw_heatmap(ax, corr_df, cols, labels, title, norm)
    plt.colorbar(im, ax=ax, label="Spearman ρ", shrink=0.7)
    plt.tight_layout()
    for ext in [".pdf", ".png"]:
        plt.savefig(out_path + ext, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}.pdf/png", flush=True)


def plot_three_heatmaps(corr_all, corr_a, corr_b, out_path, track_order=None):
    stems  = track_order or list(corr_all.columns)
    cols   = [s for s in stems if s in corr_all.columns]
    labels = [SHORT_LABELS.get(c, c) for c in cols]
    n      = len(cols)
    norm   = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    fig, axes = plt.subplots(1, 3, figsize=(n * 1.9, n * 0.72 + 1.5))
    for ax, corr_df, title in zip(axes,
                                   [corr_all, corr_a, corr_b],
                                   ["All bins", "A compartment", "B compartment"]):
        im = _draw_heatmap(ax, corr_df, cols, labels, title, norm)
    plt.colorbar(im, ax=axes[-1], label="Spearman ρ", shrink=0.8)
    plt.suptitle("Genome-wide correlation stratified by Hi-C compartment",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    for ext in [".pdf", ".png"]:
        plt.savefig(out_path + ext, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}.pdf/png", flush=True)


# =============================================================================
# TSS profile computation
# =============================================================================

def compute_profiles_for_tss_list(bw_paths, tss_list, window, n_bins, chrom_sizes):
    """
    tss_list: list of (chrom, centre, strand)  — no cluster column
    Returns: profiles dict, sems dict, per_tss dict
    """
    valid = [(c, p, s) for c, p, s in tss_list
             if c in chrom_sizes
             and p - window >= 0
             and p + window <= chrom_sizes[c]]

    profiles, sems, per_tss = {}, {}, {}
    for bw_path in bw_paths:
        stem   = os.path.splitext(os.path.basename(bw_path))[0]
        matrix = np.full((len(valid), n_bins), np.nan)
        bw     = pyBigWig.open(bw_path)
        for i, (chrom, pos, strand) in enumerate(valid):
            try:
                vals = bw.stats(chrom, pos - window, pos + window,
                                type="mean", nBins=n_bins)
                row = np.array(vals, dtype=np.float64)
                if strand == "-":
                    row = row[::-1]
                matrix[i] = row
            except Exception:
                pass
        bw.close()
        n_valid        = np.sum(~np.isnan(matrix), axis=0)
        profiles[stem] = np.nanmean(matrix, axis=0)
        sems[stem]     = np.nanstd(matrix, axis=0) / np.sqrt(np.maximum(n_valid, 1))
        per_tss[stem]  = matrix

    return profiles, sems, per_tss, len(valid)


# =============================================================================
# Per-cluster profile plots
# =============================================================================

def plot_per_cluster_profiles(all_profiles, all_sems, all_counts,
                               clusters, window, n_bins, out_path):
    """
    all_profiles : dict  cluster -> {stem -> mean_array}
    all_sems     : dict  cluster -> {stem -> sem_array}
    all_counts   : dict  cluster -> n_tss
    One column per track, one row per cluster.
    """
    x       = np.linspace(-window, window, n_bins) / 1000
    ordered = [s for s in TRACK_ORDER if s in list(all_profiles[clusters[0]].keys())]
    n_tracks   = len(ordered)
    n_clusters = len(clusters)

    fig, axes = plt.subplots(n_clusters, n_tracks,
                             figsize=(n_tracks * 2.8, n_clusters * 2.2),
                             sharex=True)

    # ensure 2D array even if 1 cluster or 1 track
    if n_clusters == 1:
        axes = axes[np.newaxis, :]
    if n_tracks == 1:
        axes = axes[:, np.newaxis]

    for row, cluster in enumerate(clusters):
        profiles = all_profiles[cluster]
        sems     = all_sems[cluster]
        n_tss    = all_counts[cluster]

        for col, stem in enumerate(ordered):
            ax    = axes[row, col]
            mean  = profiles.get(stem, np.zeros(n_bins))
            sem   = sems.get(stem, np.zeros(n_bins))
            color = TRACK_COLORS.get(stem, "#888888")

            ax.plot(x, mean, color=color, linewidth=1.5)
            ax.fill_between(x, mean - sem, mean + sem, color=color, alpha=0.2)
            ax.axvline(0, color="black", linestyle="--", linewidth=0.6)
            ax.axhline(0, color="grey",  linestyle=":",  linewidth=0.5)

            # Column header (track name) — top row only
            if row == 0:
                ax.set_title(SHORT_LABELS.get(stem, stem), fontsize=8,
                             pad=4, wrap=True)

            # Row label (cluster + n) — left column only
            if col == 0:
                ax.set_ylabel(f"{cluster}\n(n={n_tss:,})", fontsize=8)

            ax.tick_params(labelsize=6)
            ax.set_xlim(-window / 1000, window / 1000)

    # shared x label
    for col in range(n_tracks):
        axes[-1, col].set_xlabel("Distance from TSS (kb)", fontsize=7)

    plt.suptitle(f"Signal profiles by promoter cluster (±{window//1000}kb)",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    for ext in [".pdf", ".png"]:
        plt.savefig(out_path + ext, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}.pdf/png", flush=True)


def plot_overlay_by_cluster(all_profiles, clusters, window, n_bins,
                             stem, out_path):
    """For a single track, overlay all cluster mean profiles on one axis."""
    x = np.linspace(-window, window, n_bins) / 1000
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, cluster in enumerate(clusters):
        mean  = all_profiles[cluster].get(stem, np.zeros(n_bins))
        color = CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
        ax.plot(x, mean, color=color, linewidth=1.8, label=cluster)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(0, color="grey",  linestyle=":",  linewidth=0.6)
    ax.set_xlabel("Distance from TSS (kb)", fontsize=11)
    ax.set_ylabel("Mean signal", fontsize=11)
    ax.set_title(f"{SHORT_LABELS.get(stem, stem)} — by cluster", fontsize=11)
    ax.legend(fontsize=8, framealpha=0.9)
    plt.tight_layout()
    for ext in [".pdf", ".png"]:
        plt.savefig(out_path + ext, dpi=200, bbox_inches="tight")
    plt.close()


# =============================================================================
# Main
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bw-dir",           default=DEFAULT_BIGWIG_DIR)
    p.add_argument("--cluster-bed",      default=DEFAULT_CLUSTER_BED)
    p.add_argument("--outdir",           default=DEFAULT_OUTDIR)
    p.add_argument("--bin-size",         type=int, default=DEFAULT_BIN_SIZE)
    p.add_argument("--tss-window",       type=int, default=DEFAULT_TSS_WINDOW)
    p.add_argument("--tss-bins",         type=int, default=DEFAULT_TSS_BINS)
    p.add_argument("--cov-pct",          type=int, default=DEFAULT_COV_PCT)
    p.add_argument("--corr-method",      default="spearman",
                   choices=["spearman", "pearson"])
    p.add_argument("--skip-genome-wide", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "correlation_tables"), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "cluster_profiles"), exist_ok=True)

    # ── Build bigWig list ─────────────────────────────────────────────────────
    all_bws = []
    for fname in SELECTED_TRACKS:
        path = os.path.join(args.bw_dir, fname)
        if os.path.exists(path):
            all_bws.append(path)
        else:
            print(f"  WARNING: not found, skipping: {fname}")
    print(f"Using {len(all_bws)} tracks", flush=True)

    chrom_sizes = load_chrom_sizes_from_bw(all_bws[0])
    chrom_sizes = {c: s for c, s in chrom_sizes.items()
                   if not any(x in c for x in ["_", "chrM", "chrEBV", "chrUn"])}

    compartment_stem = os.path.splitext(COMPARTMENT_TRACK)[0]
    cov_filter_stem  = os.path.splitext(COV_FILTER_TRACK)[0]

    # ── PART 1: Genome-wide compartment-stratified correlation ────────────────
    if not args.skip_genome_wide:
        print("\n[PART 1] Genome-wide correlation...", flush=True)
        gw_df = build_genome_wide_matrix(all_bws, chrom_sizes, args.bin_size)

        if args.cov_pct > 0:
            gw_df = apply_coverage_filter(gw_df, cov_filter_stem, args.cov_pct)

        gw_df.to_csv(os.path.join(args.outdir, "correlation_tables",
                                  "genome_wide_raw.tsv.gz"),
                     sep="\t", index=False, compression="gzip")

        corr_all = compute_correlation_matrix(gw_df, args.corr_method)
        corr_all.to_csv(os.path.join(args.outdir, "correlation_tables",
                                     "corr_all.tsv"), sep="\t")

        df_a, df_b = stratify_by_compartment(gw_df, compartment_stem)
        corr_a = compute_correlation_matrix(df_a, args.corr_method)
        corr_b = compute_correlation_matrix(df_b, args.corr_method)
        corr_a.to_csv(os.path.join(args.outdir, "correlation_tables",
                                   "corr_A_compartment.tsv"), sep="\t")
        corr_b.to_csv(os.path.join(args.outdir, "correlation_tables",
                                   "corr_B_compartment.tsv"), sep="\t")

        for corr, tag in [(corr_all, "all"),
                          (corr_a,   "A_compartment"),
                          (corr_b,   "B_compartment")]:
            plot_correlation_heatmap(
                corr,
                title=f"Genome-wide — {tag.replace('_',' ')} "
                      f"({args.corr_method}, {args.bin_size//1000}kb)",
                out_path=os.path.join(args.outdir, f"genome_wide_corr_{tag}"),
                track_order=TRACK_ORDER
            )
        plot_three_heatmaps(
            corr_all, corr_a, corr_b,
            out_path=os.path.join(args.outdir, "genome_wide_corr_stratified"),
            track_order=TRACK_ORDER
        )

    # ── PART 2: Per-cluster TSS profiles ─────────────────────────────────────
    print("\n[PART 2] Loading cluster BED...", flush=True)
    tss_all, clusters = load_cluster_bed(args.cluster_bed, chrom_sizes)

    all_profiles = {}
    all_sems     = {}
    all_counts   = {}
    all_per_tss  = {}

    for cluster in clusters:
        tss_subset = [(c, p, s) for c, p, s, cl in tss_all if cl == cluster]
        print(f"\n  Cluster {cluster}: {len(tss_subset):,} TSSs", flush=True)

        profiles, sems, per_tss, n_valid = compute_profiles_for_tss_list(
            all_bws, tss_subset, args.tss_window, args.tss_bins, chrom_sizes)

        all_profiles[cluster] = profiles
        all_sems[cluster]     = sems
        all_counts[cluster]   = n_valid
        all_per_tss[cluster]  = per_tss

    # Grid plot: clusters × tracks
    print("\n  Plotting per-cluster grid...", flush=True)
    plot_per_cluster_profiles(
        all_profiles, all_sems, all_counts, clusters,
        args.tss_window, args.tss_bins,
        out_path=os.path.join(args.outdir, "tss_profiles_by_cluster")
    )

    # Per-track overlay plots (all clusters on one axis)
    print("  Plotting per-track cluster overlays...", flush=True)
    overlay_dir = os.path.join(args.outdir, "cluster_profiles")
    for stem in TRACK_ORDER:
        if stem in all_profiles[clusters[0]]:
            plot_overlay_by_cluster(
                all_profiles, clusters, args.tss_window, args.tss_bins,
                stem=stem,
                out_path=os.path.join(overlay_dir, f"overlay_{stem}")
            )

    # ── PART 3: Per-cluster correlation matrices ──────────────────────────────
    print("\n[PART 3] Per-cluster correlation matrices...", flush=True)
    for cluster in clusters:
        tss_df = pd.DataFrame({
            stem: np.nanmean(mat, axis=1)
            for stem, mat in all_per_tss[cluster].items()
        })
        tss_corr = compute_correlation_matrix(tss_df, args.corr_method)
        tss_corr.to_csv(os.path.join(args.outdir, "correlation_tables",
                                     f"corr_{cluster}.tsv"), sep="\t")
        plot_correlation_heatmap(
            tss_corr,
            title=f"TSS correlation — {cluster} (n={all_counts[cluster]:,}, "
                  f"±{args.tss_window//1000}kb)",
            out_path=os.path.join(args.outdir, f"tss_corr_{cluster}"),
            track_order=TRACK_ORDER
        )

    print(f"\n  Done! Outputs in: {args.outdir}")
    print("    tss_profiles_by_cluster.pdf   — grid: clusters × tracks")
    print("    cluster_profiles/overlay_*.pdf — per-track, all clusters overlaid")
    print("    tss_corr_cluster_N.pdf         — correlation per cluster")
    print("    genome_wide_corr_stratified.pdf")


if __name__ == "__main__":
    main()
