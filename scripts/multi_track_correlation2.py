#!/usr/bin/env python3
"""
multi_track_correlation2.py  (v2)
==================================
1. Genome-wide correlation matrix across all bigWig tracks at a chosen bin size
2. TSS-centered profile plots (mean signal ± SEM across all TSSs)
3. TSS-region correlation matrix (same tracks, restricted to TSS ± window)

Outputs:
  - genome_wide_correlation.pdf/png
  - tss_correlation.pdf/png
  - tss_profiles.pdf/png
  - correlation_tables/  (TSV files of the raw matrices)
  - tss_matrix.tab       (raw per-TSS values for custom analysis)
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from concurrent.futures import ThreadPoolExecutor

# =============================================================================
# Configuration — edit these paths if not passing via command line
# =============================================================================

DEFAULT_BIGWIG_DIR = "."
DEFAULT_TSS_BED    = "T2T_TSS_1kb_reduced_clean.bed"
DEFAULT_OUTDIR     = "./correlation_output"
DEFAULT_BIN_SIZE   = 10_000   # bp for genome-wide binning
DEFAULT_TSS_WINDOW = 5_000    # bp each side of TSS centre
DEFAULT_TSS_BINS   = 100      # number of bins across TSS window

# ── Reduced track selection — only these files will be used ──────────────────
SELECTED_TRACKS = [
    "MCF10A_subcompartments_20000_T2T.bw",
    "A_Inp_rep1_distribution.bw",
    "A_Inp_rep1_sensitivity.bw",
    "MCF10A_ATACseq_rep1_SRR12006157_coverage.bw",
    "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T.bw",
    "20U_supernatant_rep1_1kb_inverted.bw",
    "20U_supernatant_rep2_1kb_inverted.bw",
]

# Track groupings (used for colour-coding TSS profiles)
TRACK_GROUPS = {
    "Hi-C compartments":               ["MCF10A_subcompartments_20000_T2T.bw"],
    "Nucleosome distribution/sensitivity": [
        "A_Inp_rep1_distribution.bw",
        "A_Inp_rep1_sensitivity.bw",
    ],
    "ATAC-seq":                        ["MCF10A_ATACseq_rep1_SRR12006157_coverage.bw"],
    "Nucleosome (external)":           ["GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T.bw"],
    "Condensability (20U supernatant rep1)": ["20U_supernatant_rep1_1kb_inverted.bw"],
    "Condensability (20U supernatant rep2)": ["20U_supernatant_rep2_1kb_inverted.bw"],
}

# Short display labels (filename stem -> short label)
SHORT_LABELS = {
    "MCF10A_subcompartments_20000_T2T":              "Hi-C compartments",
    "A_Inp_rep1_distribution":                       "Nuc distribution",
    "A_Inp_rep1_sensitivity":                        "Nuc sensitivity",
    "MCF10A_ATACseq_rep1_SRR12006157_coverage":      "ATAC-seq",
    "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T":   "GM NCP condensability",
    "20U_supernatant_rep1_1kb_inverted":             "Cond 20U supernatant rep1",
    "20U_supernatant_rep2_1kb_inverted":             "Cond 20U supernatant rep2",
}

# Preferred track order for plots (stems, without .bw)
TRACK_ORDER = [
    "MCF10A_subcompartments_20000_T2T",
    "A_Inp_rep1_distribution",
    "A_Inp_rep1_sensitivity",
    "MCF10A_ATACseq_rep1_SRR12006157_coverage",
    "GSE252941_GM_NCP_sp_1rep_deep_1kb_score_T2T",
    "20U_supernatant_rep1_1kb_inverted",
    "20U_supernatant_rep2_1kb_inverted",
]


# =============================================================================
# I/O helpers
# =============================================================================

def load_chrom_sizes_from_bw(bw_path):
    bw = pyBigWig.open(bw_path)
    sizes = bw.chroms()
    bw.close()
    return sizes


def fetch_bw_values(bw_path, chrom, start, end, bins=None):
    """Fetch mean values from a bigWig, optionally binned."""
    try:
        bw = pyBigWig.open(bw_path)
        if bins:
            vals = bw.stats(chrom, start, end, type="mean", nBins=bins)
        else:
            vals = bw.stats(chrom, start, end, type="mean", nBins=(end - start))
        bw.close()
        return np.array(vals, dtype=np.float64)
    except Exception:
        if bins:
            return np.full(bins, np.nan)
        return np.full(end - start, np.nan)


# =============================================================================
# PART 1: Genome-wide binned correlation
# =============================================================================

def build_genome_wide_matrix(bw_paths, chrom_sizes, bin_size, min_chroms=None):
    """
    For each chromosome, tile into bin_size windows and fetch mean signal
    from each bigWig. Returns a DataFrame (rows=bins, cols=tracks).
    """
    chroms = [c for c, s in chrom_sizes.items()
              if s > bin_size * 10  # skip tiny contigs
              and not any(x in c for x in ["_", "chrM", "chrEBV"])]

    if min_chroms:
        chroms = chroms[:min_chroms]

    print(f"  Using {len(chroms)} chromosomes for genome-wide analysis", flush=True)

    all_rows = {os.path.splitext(os.path.basename(p))[0]: [] for p in bw_paths}

    for chrom in chroms:
        chrom_len = chrom_sizes[chrom]
        bins_starts = np.arange(0, chrom_len, bin_size)
        print(f"    {chrom}: {len(bins_starts)} bins", flush=True)

        for bw_path in bw_paths:
            stem = os.path.splitext(os.path.basename(bw_path))[0]
            try:
                bw = pyBigWig.open(bw_path)
                vals = bw.stats(chrom, 0, chrom_len, type="mean",
                                nBins=len(bins_starts))
                bw.close()
                all_rows[stem].extend(vals if vals else [np.nan] * len(bins_starts))
            except Exception as e:
                all_rows[stem].extend([np.nan] * len(bins_starts))

    df = pd.DataFrame(all_rows)
    return df


def compute_correlation_matrix(df, method="spearman"):
    """Compute pairwise correlation matrix, dropping bins with any NaN."""
    df_clean = df.dropna()
    print(f"  {len(df_clean):,} / {len(df):,} bins retained after dropping NaNs", flush=True)
    if method == "spearman":
        corr = df_clean.corr(method="spearman")
    else:
        corr = df_clean.corr(method="pearson")
    return corr


def plot_correlation_heatmap(corr_df, title, out_path, track_order=None):
    """Plot a publication-quality correlation heatmap."""
    # Reorder columns/rows
    stems = [os.path.splitext(os.path.basename(t))[0] if t.endswith(".bw") else t
             for t in (track_order or [])]
    cols = [s for s in stems if s in corr_df.columns] + \
           [c for c in corr_df.columns if c not in stems]
    corr_ordered = corr_df.loc[cols, cols]

    labels = [SHORT_LABELS.get(c, c) for c in cols]

    n = len(cols)
    fig_size = max(8, n * 0.55)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.85))

    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im = ax.imshow(corr_ordered.values, cmap="RdBu_r", norm=norm, aspect="auto")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)

    # Annotate cells
    for i in range(n):
        for j in range(n):
            val = corr_ordered.values[i, j]
            color = "white" if abs(val) > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=6, color=color)

    plt.colorbar(im, ax=ax, label="Spearman ρ", shrink=0.7)
    ax.set_title(title, fontsize=12, pad=15)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}", flush=True)


# =============================================================================
# PART 2: TSS profiles
# =============================================================================

def load_tss_bed(bed_path, window):
    """
    Load TSS BED file. Expects at least 3 columns (chrom, start, end).
    If 6 columns, uses strand (col 5) to orient the window.
    Returns list of (chrom, centre, strand).
    """
    tss_list = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start, end = int(parts[1]), int(parts[2])
            centre = (start + end) // 2
            strand = parts[5] if len(parts) >= 6 else "+"
            tss_list.append((chrom, centre, strand))
    print(f"  Loaded {len(tss_list):,} TSSs from {bed_path}", flush=True)
    return tss_list


def compute_tss_profiles(bw_paths, tss_list, window, n_bins, chrom_sizes):
    """
    For each bigWig, compute the mean profile across all TSSs.
    Returns dict: stem -> (mean_profile, sem_profile) arrays of length n_bins.
    Also returns per-TSS matrix for correlation.
    """
    profiles  = {}   # stem -> mean array
    sems      = {}   # stem -> SEM array
    per_tss   = {}   # stem -> 2D array (n_tss, n_bins)

    valid_tss = [(c, pos, strand) for c, pos, strand in tss_list
                 if c in chrom_sizes
                 and pos - window >= 0
                 and pos + window <= chrom_sizes[c]]
    print(f"  {len(valid_tss):,} TSSs within chromosome bounds", flush=True)

    for bw_path in bw_paths:
        stem = os.path.splitext(os.path.basename(bw_path))[0]
        print(f"    Profiling {stem}...", flush=True)
        matrix = np.full((len(valid_tss), n_bins), np.nan)

        bw = pyBigWig.open(bw_path)
        for i, (chrom, pos, strand) in enumerate(valid_tss):
            s = pos - window
            e = pos + window
            try:
                vals = bw.stats(chrom, s, e, type="mean", nBins=n_bins)
                row = np.array(vals, dtype=np.float64)
                if strand == "-":
                    row = row[::-1]
                matrix[i] = row
            except Exception:
                pass
        bw.close()

        # Mean and SEM ignoring NaN
        mean_prof = np.nanmean(matrix, axis=0)
        sem_prof  = np.nanstd(matrix, axis=0) / np.sqrt(np.sum(~np.isnan(matrix), axis=0))
        profiles[stem]  = mean_prof
        sems[stem]      = sem_prof
        per_tss[stem]   = matrix

    return profiles, sems, per_tss, valid_tss


def plot_tss_profiles(profiles, sems, window, n_bins, out_path, track_order=None):
    """Plot TSS profile lines with shaded SEM, grouped by track type."""
    x = np.linspace(-window, window, n_bins)

    # Group colours
    group_colors = {
        "Hi-C compartments":                    "#7B2D8B",
        "Nucleosome distribution/sensitivity":  "#2166AC",
        "ATAC-seq":                             "#1A9850",
        "Nucleosome (external)":                "#D6604D",
        "Condensability (20U supernatant rep1)":"#B2182B",
        "Condensability (20U supernatant rep2)":"#D6604D",
    }

    # Build ordered list
    stems_ordered = [os.path.splitext(os.path.basename(t))[0]
                     for t in (track_order or [])]
    stems_ordered = [s for s in stems_ordered if s in profiles] + \
                    [s for s in profiles if s not in stems_ordered]

    # Map stem -> group colour
    stem_color = {}
    for group, fnames in TRACK_GROUPS.items():
        for fname in fnames:
            stem = os.path.splitext(fname)[0]
            stem_color[stem] = group_colors.get(group, "#555555")

    fig, ax = plt.subplots(figsize=(10, 6))

    for stem in stems_ordered:
        if stem not in profiles:
            continue
        mean = profiles[stem]
        sem  = sems[stem]
        label = SHORT_LABELS.get(stem, stem)
        color = stem_color.get(stem, "#888888")

        # Normalise to z-score for comparable plotting
        std = np.nanstd(mean)
        if std > 0:
            mean_z = (mean - np.nanmean(mean)) / std
            sem_z  = sem / std
        else:
            mean_z = mean
            sem_z  = sem

        ax.plot(x / 1000, mean_z, label=label, color=color, linewidth=1.5)
        ax.fill_between(x / 1000, mean_z - sem_z, mean_z + sem_z,
                        color=color, alpha=0.15)

    ax.axvline(0, color="black", linestyle="--", linewidth=0.8, label="TSS")
    ax.set_xlabel("Distance from TSS (kb)", fontsize=11)
    ax.set_ylabel("Z-score (mean ± SEM)", fontsize=11)
    ax.set_title("Signal profiles around TSSs", fontsize=13)
    ax.legend(fontsize=7, ncol=2, loc="upper right",
              framealpha=0.9, bbox_to_anchor=(1.35, 1))
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_path}", flush=True)


def build_tss_mean_df(per_tss):
    """
    Collapse per-TSS matrices to per-TSS mean values (one value per TSS per track).
    Returns DataFrame (rows=TSSs, cols=tracks) for correlation at TSS regions.
    """
    data = {stem: np.nanmean(mat, axis=1) for stem, mat in per_tss.items()}
    return pd.DataFrame(data)


# =============================================================================
# Main
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bw-dir",    default=DEFAULT_BIGWIG_DIR)
    p.add_argument("--tss-bed",   default=DEFAULT_TSS_BED)
    p.add_argument("--outdir",    default=DEFAULT_OUTDIR)
    p.add_argument("--bin-size",  type=int, default=DEFAULT_BIN_SIZE)
    p.add_argument("--tss-window",type=int, default=DEFAULT_TSS_WINDOW)
    p.add_argument("--tss-bins",  type=int, default=DEFAULT_TSS_BINS)
    p.add_argument("--corr-method", default="spearman", choices=["spearman","pearson"])
    p.add_argument("--skip-genome-wide", action="store_true",
                   help="Skip genome-wide correlation (faster, TSS only)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "correlation_tables"), exist_ok=True)

    # Use only the selected tracks, in the defined order
    all_bws = []
    for fname in SELECTED_TRACKS:
        path = os.path.join(args.bw_dir, fname)
        if os.path.exists(path):
            all_bws.append(path)
        else:
            print(f"  WARNING: selected track not found, skipping: {fname}", flush=True)
    print(f"Using {len(all_bws)} selected bigWig files:", flush=True)
    for b in all_bws:
        print(f"  {os.path.basename(b)}")

    # Get chrom sizes from first bigWig
    chrom_sizes = load_chrom_sizes_from_bw(all_bws[0])
    # Filter to main chromosomes only
    chrom_sizes = {c: s for c, s in chrom_sizes.items()
                   if not any(x in c for x in ["_", "chrM", "chrEBV", "chrUn"])}

    # ── PART 1: Genome-wide correlation ──────────────────────────────────────
    if not args.skip_genome_wide:
        print("\n[PART 1] Building genome-wide binned matrix...", flush=True)
        gw_df = build_genome_wide_matrix(all_bws, chrom_sizes, args.bin_size)

        # Save raw matrix
        gw_tsv = os.path.join(args.outdir, "correlation_tables", "genome_wide_raw.tsv.gz")
        gw_df.to_csv(gw_tsv, sep="\t", index=False, compression="gzip")
        print(f"  Raw matrix saved: {gw_tsv}")

        print(f"\n  Computing {args.corr_method} correlation...", flush=True)
        gw_corr = compute_correlation_matrix(gw_df, method=args.corr_method)
        gw_corr.to_csv(os.path.join(args.outdir, "correlation_tables",
                                    "genome_wide_correlation.tsv"), sep="\t")

        plot_correlation_heatmap(
            gw_corr,
            title=f"Genome-wide correlation ({args.corr_method}, {args.bin_size//1000}kb bins)",
            out_path=os.path.join(args.outdir, "genome_wide_correlation.pdf"),
            track_order=TRACK_ORDER
        )
        plot_correlation_heatmap(
            gw_corr,
            title=f"Genome-wide correlation ({args.corr_method}, {args.bin_size//1000}kb bins)",
            out_path=os.path.join(args.outdir, "genome_wide_correlation.png"),
            track_order=TRACK_ORDER
        )

    # ── PART 2: TSS profiles and correlation ─────────────────────────────────
    print("\n[PART 2] Computing TSS profiles...", flush=True)
    tss_list = load_tss_bed(args.tss_bed, args.tss_window)

    profiles, sems, per_tss, valid_tss = compute_tss_profiles(
        all_bws, tss_list, args.tss_window, args.tss_bins, chrom_sizes
    )

    # Profile plot
    plot_tss_profiles(
        profiles, sems, args.tss_window, args.tss_bins,
        out_path=os.path.join(args.outdir, "tss_profiles.pdf"),
        track_order=TRACK_ORDER
    )
    plot_tss_profiles(
        profiles, sems, args.tss_window, args.tss_bins,
        out_path=os.path.join(args.outdir, "tss_profiles.png"),
        track_order=TRACK_ORDER
    )

    # TSS-level correlation (one mean value per TSS per track)
    print("\n  Computing TSS-level correlation...", flush=True)
    tss_df = build_tss_mean_df(per_tss)
    tss_df.to_csv(os.path.join(args.outdir, "correlation_tables",
                               "tss_per_tss_values.tsv.gz"),
                  sep="\t", index=False, compression="gzip")

    tss_corr = compute_correlation_matrix(tss_df, method=args.corr_method)
    tss_corr.to_csv(os.path.join(args.outdir, "correlation_tables",
                                 "tss_correlation.tsv"), sep="\t")

    plot_correlation_heatmap(
        tss_corr,
        title=f"TSS-region correlation ({args.corr_method}, ±{args.tss_window//1000}kb)",
        out_path=os.path.join(args.outdir, "tss_correlation.pdf"),
        track_order=TRACK_ORDER
    )
    plot_correlation_heatmap(
        tss_corr,
        title=f"TSS-region correlation ({args.corr_method}, ±{args.tss_window//1000}kb)",
        out_path=os.path.join(args.outdir, "tss_correlation.png"),
        track_order=TRACK_ORDER
    )

    print("\n================================================")
    print(f" Done! All outputs in: {args.outdir}")
    print("   genome_wide_correlation.pdf/png")
    print("   tss_profiles.pdf/png")
    print("   tss_correlation.pdf/png")
    print("   correlation_tables/")
    print("================================================")


if __name__ == "__main__":
    main()
