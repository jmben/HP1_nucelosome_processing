#!/usr/bin/env python3
"""
replot_cluster_grid.py
======================
Regenerates the per-cluster TSS profile grid plot without recomputing
any bigWig data. Reads directly from the bigWig files and cluster BED.

Usage:
    python3 replot_cluster_grid.py \
        --bw-dir /data2/jbenoit/condense-seq/condensability_output/bigwig \
        --cluster-bed /pauper/jbenoit/FINAL_HP1_FIGS/HP1_HPTM_2kb_sen_matrix_og_sort_k7.bed \
        --outdir ./correlation_output_v3 \
        --tss-window 1000 \
        --tss-bins 80
"""

import os
import argparse
import numpy as np
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# Config — edit these if needed
# =============================================================================

DEFAULT_BIGWIG_DIR  = "."
DEFAULT_CLUSTER_BED = "HP1_HPTM_2kb_sen_matrix_og_sort_k7.bed"
DEFAULT_OUTDIR      = "./correlation_output_v3"
DEFAULT_TSS_WINDOW  = 1_000
DEFAULT_TSS_BINS    = 80

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


# =============================================================================
# Helpers
# =============================================================================

def load_chrom_sizes_from_bw(bw_path):
    bw = pyBigWig.open(bw_path)
    sizes = dict(bw.chroms())
    bw.close()
    return sizes


def load_cluster_bed(bed_path, chrom_sizes):
    tss_all, skipped = [], 0
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom   = parts[0]
            centre  = (int(parts[1]) + int(parts[2])) // 2
            strand  = parts[5] if len(parts) >= 6 else "+"
            cluster = parts[-1].strip()
            if chrom not in chrom_sizes:
                skipped += 1
                continue
            tss_all.append((chrom, centre, strand, cluster))
    clusters = sorted(set(t[3] for t in tss_all))
    print(f"  Loaded {len(tss_all):,} entries ({skipped} skipped), "
          f"{len(clusters)} clusters: {clusters}", flush=True)
    return tss_all, clusters


def compute_profiles(bw_paths, tss_list, window, n_bins, chrom_sizes):
    """tss_list: list of (chrom, centre, strand)"""
    valid = [(c, p, s) for c, p, s in tss_list
             if c in chrom_sizes
             and p - window >= 0
             and p + window <= chrom_sizes[c]]

    profiles, sems = {}, {}
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

    return profiles, sems, len(valid)


# =============================================================================
# Plot
# =============================================================================

def plot_grid(all_profiles, all_sems, all_counts, clusters,
              ordered, window, n_bins, out_path):
    x          = np.linspace(-window, window, n_bins) / 1000
    n_tracks   = len(ordered)
    n_clusters = len(clusters)

    fig, axes = plt.subplots(n_clusters, n_tracks,
                             figsize=(n_tracks * 2.8, n_clusters * 2.2),
                             sharex=True)
    if n_clusters == 1:
        axes = axes[np.newaxis, :]
    if n_tracks == 1:
        axes = axes[:, np.newaxis]

    # Compute shared y limits per column
    col_ylims = {}
    for col, stem in enumerate(ordered):
        all_vals = []
        for cluster in clusters:
            mean = all_profiles[cluster].get(stem, np.zeros(n_bins))
            sem  = all_sems[cluster].get(stem, np.zeros(n_bins))
            all_vals.extend((mean - sem).tolist())
            all_vals.extend((mean + sem).tolist())
        all_vals = [v for v in all_vals if not np.isnan(v)]
        if all_vals:
            ymin, ymax = min(all_vals), max(all_vals)
            pad = (ymax - ymin) * 0.08
            col_ylims[col] = (ymin - pad, ymax + pad)
        else:
            col_ylims[col] = (0, 1)

    for row, cluster in enumerate(clusters):
        for col, stem in enumerate(ordered):
            ax    = axes[row, col]
            mean  = all_profiles[cluster].get(stem, np.zeros(n_bins))
            sem   = all_sems[cluster].get(stem, np.zeros(n_bins))
            color = TRACK_COLORS.get(stem, "#888888")

            ax.plot(x, mean, color=color, linewidth=1.5)
            ax.fill_between(x, mean - sem, mean + sem, color=color, alpha=0.2)
            ax.axvline(0, color="black", linestyle="--", linewidth=0.6)
            ax.axhline(0, color="grey",  linestyle=":",  linewidth=0.5)
            ax.set_ylim(col_ylims[col])
            ax.set_xlim(-window / 1000, window / 1000)
            ax.tick_params(labelsize=6)

            if row == 0:
                ax.set_title(SHORT_LABELS.get(stem, stem), fontsize=8, pad=4)
            if col == 0:
                ax.set_ylabel(f"{cluster}\n(n={all_counts[cluster]:,})", fontsize=8)

    for col in range(n_tracks):
        axes[-1, col].set_xlabel("Distance from TSS (kb)", fontsize=7)

    plt.suptitle(f"Signal profiles by promoter cluster (±{window//1000}kb)",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    for ext in [".pdf", ".png"]:
        plt.savefig(out_path + ext, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}.pdf/png", flush=True)


# =============================================================================
# Main
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bw-dir",      default=DEFAULT_BIGWIG_DIR)
    p.add_argument("--cluster-bed", default=DEFAULT_CLUSTER_BED)
    p.add_argument("--outdir",      default=DEFAULT_OUTDIR)
    p.add_argument("--tss-window",  type=int, default=DEFAULT_TSS_WINDOW)
    p.add_argument("--tss-bins",    type=int, default=DEFAULT_TSS_BINS)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Build bigWig list
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

    print("\nLoading cluster BED...", flush=True)
    tss_all, clusters = load_cluster_bed(args.cluster_bed, chrom_sizes)

    ordered = [s for s in TRACK_ORDER
               if any(os.path.splitext(os.path.basename(b))[0] == s for b in all_bws)]

    all_profiles, all_sems, all_counts = {}, {}, {}
    for cluster in clusters:
        tss_subset = [(c, p, s) for c, p, s, cl in tss_all if cl == cluster]
        print(f"\n  {cluster}: {len(tss_subset):,} TSSs", flush=True)
        profiles, sems, n_valid = compute_profiles(
            all_bws, tss_subset, args.tss_window, args.tss_bins, chrom_sizes)
        all_profiles[cluster] = profiles
        all_sems[cluster]     = sems
        all_counts[cluster]   = n_valid

    print("\nPlotting...", flush=True)
    plot_grid(all_profiles, all_sems, all_counts, clusters,
              ordered, args.tss_window, args.tss_bins,
              out_path=os.path.join(args.outdir, "tss_profiles_by_cluster"))


if __name__ == "__main__":
    main()
