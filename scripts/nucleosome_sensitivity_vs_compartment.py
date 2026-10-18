#!/usr/bin/env python3
"""
nucleosome_sensitivity_vs_compartment.py
=============================================================================
Assigns each nucleosome sensitivity call to a Hi-C subcompartment based on
the TSS position of its associated gene, then tests whether sensitivity
(Sensitive/Resistant/Neutral) is enriched in specific compartments.

Input files:
  - TSS_10A_*_sensitivity_overview.txt
    Columns: nuc_chr nuc_start nuc_end sensitivity gene_chr gene_start gene_end
             gene_name score strand
  - MCF10A_subcompartments_20000_T2T.bg
    Columns: chrom start end subcompartment score strand

Usage:
    python nucleosome_sensitivity_vs_compartment.py \
        --bedgraph MCF10A_subcompartments_20000_T2T.bg \
        --files    TSS_10A_*txt \
        --out_dir  nucleosome_results/
=============================================================================
"""

import argparse
import glob
import numpy as np
from pathlib import Path
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.gridspec import GridSpec
except ImportError:
    raise ImportError("pip install matplotlib")

try:
    from scipy.stats import chi2_contingency, fisher_exact
except ImportError:
    raise ImportError("pip install scipy")


# Subcompartment order: most active → most inactive
SUBCOMPS = ["A3", "A2", "A1", "A0", "B0", "B1", "B2", "B3"]
SENS_CALLS = ["Resistant", "Neutral", "Sensitive"]  # bottom to top in stacked bar

SUBCOMP_COLORS = {
    "A3": "#c0392b", "A2": "#e74c3c", "A1": "#f39c12", "A0": "#f1c40f",
    "B0": "#1a5276", "B1": "#2980b9", "B2": "#7fb3d3", "B3": "#aed6f1",
}
SENS_COLORS = {
    "Sensitive": "#E69F00",   # gold
    "Neutral":   "#D5D8DC",   # light grey
    "Resistant": "#1B3A6B",   # navy
}

# Sample groupings
SAMPLE_GROUPS = {
    "Scr":  ["TSS_10A_Scr_clone1_sensitivity_overview.txt",
             "TSS_10A_Scr_clone2_sensitivity_overview.txt"],
    "HP1a": ["TSS_10A_HP1a_clone5_sensitivity_overview.txt",
             "TSS_10A_HP1a_clone6_sensitivity_overview.txt"],
    "HP1b": ["TSS_10A_HP1b_clone15_sensitivity_overview.txt",
             "TSS_10A_HP1b_cloneD_sensitivity_overview.txt"],
}


# ---------------------------------------------------------------------------
# Build subcompartment interval tree (simple sorted list for overlap)
# ---------------------------------------------------------------------------
def load_subcompartments(bg_path: str) -> dict:
    """Returns dict: chrom -> sorted list of (start, end, label)"""
    comp = defaultdict(list)
    with open(bg_path) as f:
        for line in f:
            p = line.strip().split()
            if len(p) < 4:
                continue
            chrom, start, end, label = p[0], int(p[1]), int(p[2]), p[3]
            comp[chrom].append((start, end, label))
    for chrom in comp:
        comp[chrom].sort()
    print(f"  Loaded {sum(len(v) for v in comp.values()):,} subcompartment intervals")
    return dict(comp)


def query_subcompartment(comp: dict, chrom: str, pos: int) -> str:
    """Binary search for subcompartment label at position."""
    intervals = comp.get(chrom, [])
    if not intervals:
        return None
    lo, hi = 0, len(intervals) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e, label = intervals[mid]
        if pos < s:
            hi = mid - 1
        elif pos >= e:
            lo = mid + 1
        else:
            return label
    return None


# ---------------------------------------------------------------------------
# Load sensitivity file and assign subcompartments
# ---------------------------------------------------------------------------
def load_sensitivity_file(fpath: str, comp: dict) -> list:
    """
    Returns list of dicts with keys:
      nuc_chr, nuc_start, nuc_end, sensitivity,
      gene_chr, gene_start, gene_end, gene_name,
      tss_pos, subcompartment
    """
    records = []
    with open(fpath) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 10:
                p = line.strip().split()
            if len(p) < 10:
                continue
            nuc_chr   = p[0]
            nuc_start = int(p[1])
            nuc_end   = int(p[2])
            sens      = p[3]
            gene_chr  = p[4]
            gene_start = int(p[5])
            gene_end   = int(p[6])
            gene_name  = p[7]
            strand    = p[9] if len(p) > 9 else "."

            if sens not in SENS_CALLS:
                continue

            # TSS = gene_start if +, gene_end if -
            tss = gene_start if strand == "+" else gene_end
            sc  = query_subcompartment(comp, gene_chr, tss)
            if sc is None:
                continue

            records.append({
                "nuc_chr":   nuc_chr,
                "nuc_start": nuc_start,
                "nuc_end":   nuc_end,
                "sensitivity": sens,
                "gene_name": gene_name,
                "tss_pos":   tss,
                "subcompartment": sc,
            })
    return records


# ---------------------------------------------------------------------------
# Aggregate: count sensitivity calls per subcompartment
# ---------------------------------------------------------------------------
def count_by_subcomp(records: list) -> dict:
    """Returns dict: subcomp -> {Sensitive: n, Neutral: n, Resistant: n}"""
    counts = {sc: {s: 0 for s in SENS_CALLS} for sc in SUBCOMPS}
    for r in records:
        sc   = r["subcompartment"]
        sens = r["sensitivity"]
        if sc in counts:
            counts[sc][sens] += 1
    return counts


# ---------------------------------------------------------------------------
# Statistical test: chi-squared per subcompartment vs rest
# ---------------------------------------------------------------------------
def test_enrichment(counts: dict) -> dict:
    """
    For each subcompartment, test whether sensitivity distribution
    differs from the genome-wide average using chi-squared.
    Returns dict: subcomp -> (chi2, p, sensitive_pct, resistant_pct)
    """
    # Genome-wide totals
    total = {s: sum(counts[sc][s] for sc in SUBCOMPS) for s in SENS_CALLS}
    grand_total = sum(total.values())

    results = {}
    for sc in SUBCOMPS:
        sc_total = sum(counts[sc].values())
        if sc_total == 0:
            continue
        # 2×3 contingency: [sc vs rest] × [Sensitive, Neutral, Resistant]
        rest  = {s: total[s] - counts[sc][s] for s in SENS_CALLS}
        table = [[counts[sc][s] for s in SENS_CALLS],
                 [rest[s]       for s in SENS_CALLS]]
        chi2, p, dof, expected = chi2_contingency(table)
        sens_pct = 100 * counts[sc]["Sensitive"] / sc_total
        res_pct  = 100 * counts[sc]["Resistant"] / sc_total
        results[sc] = (chi2, p, sens_pct, res_pct, sc_total)
    return results


# ---------------------------------------------------------------------------
# Plot 1b: Grouped bar — Sensitive vs Resistant counts per subcompartment
# ---------------------------------------------------------------------------
def plot_grouped_bars(counts_dict: dict, out_path: Path):
    """
    Side-by-side bars: Sensitive (gold) and Resistant (navy) per subcompartment.
    Neutral shown as thin grey bar for context.
    Significance stars above each pair vs genome-wide average.
    One panel per group (Scr, HP1a, HP1b).
    """
    from scipy.stats import chi2_contingency

    n_groups = len(counts_dict)
    fig, axes = plt.subplots(1, n_groups,
                             figsize=(5.5 * n_groups, 5.5),
                             sharey=False)
    if n_groups == 1:
        axes = [axes]

    GOLD  = "#D4A017"   # Sensitive
    NAVY  = "#1B3A6B"   # Resistant
    LGREY = "#D5D8DC"   # Neutral

    BAR_W = 0.28
    x     = np.arange(len(SUBCOMPS))

    for ax, (group, counts) in zip(axes, counts_dict.items()):

        # Genome-wide totals for significance test
        total_s = sum(counts[sc]["Sensitive"] for sc in SUBCOMPS)
        total_r = sum(counts[sc]["Resistant"] for sc in SUBCOMPS)
        total_n = sum(counts[sc]["Neutral"]   for sc in SUBCOMPS)
        grand   = total_s + total_r + total_n

        sens_vals = np.array([counts[sc]["Sensitive"] for sc in SUBCOMPS], float)
        res_vals  = np.array([counts[sc]["Resistant"] for sc in SUBCOMPS], float)
        neut_vals = np.array([counts[sc]["Neutral"]   for sc in SUBCOMPS], float)
        totals    = sens_vals + res_vals + neut_vals

        # Draw bars
        b_neut = ax.bar(x,          neut_vals, BAR_W * 3,
                        color=LGREY, alpha=0.35, label="Neutral",
                        edgecolor="none", zorder=1)
        b_sens = ax.bar(x - BAR_W/2, sens_vals, BAR_W,
                        color=GOLD, label="Sensitive",
                        edgecolor="white", linewidth=0.4, zorder=3)
        b_res  = ax.bar(x + BAR_W/2, res_vals,  BAR_W,
                        color=NAVY, label="Resistant",
                        edgecolor="white", linewidth=0.4, zorder=3)

        # Significance stars above each subcompartment pair
        ymax = max(max(sens_vals), max(res_vals)) * 1.05
        for i, sc in enumerate(SUBCOMPS):
            sc_total = totals[i]
            if sc_total < 10:
                continue
            rest_s = total_s - counts[sc]["Sensitive"]
            rest_r = total_r - counts[sc]["Resistant"]
            rest_n = total_n - counts[sc]["Neutral"]
            table = [
                [counts[sc]["Sensitive"], counts[sc]["Resistant"], counts[sc]["Neutral"]],
                [rest_s, rest_r, rest_n]
            ]
            try:
                _, p, _, _ = chi2_contingency(table)
            except Exception:
                p = 1.0
            sig = "***" if p < 0.001 else ("**" if p < 0.01 else
                  ("*" if p < 0.05 else ""))
            if sig:
                bar_top = max(sens_vals[i], res_vals[i]) * 1.03
                ax.text(x[i], bar_top, sig, ha="center", va="bottom",
                        fontsize=9, color="#2c3e50", fontweight="bold")

        # A|B divider
        ax.axvline(3.5, color="black", lw=1.0, ls="--", alpha=0.5)
        ax.text(1.5, ymax * 1.12, "← A compartment →",
                ha="center", fontsize=8, color="#922b21")
        ax.text(5.5, ymax * 1.12, "← B compartment →",
                ha="center", fontsize=8, color="#1a5276")

        ax.set_xticks(x)
        ax.set_xticklabels(SUBCOMPS, fontsize=9)
        ax.set_xlabel("Subcompartment\n(A3 = most active → B3 = most inactive)",
                      fontsize=9)
        ax.set_ylabel("Nucleosome count", fontsize=10)
        ax.set_title(group, fontsize=12, fontweight="bold")
        ax.set_ylim(0, ymax * 1.18)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if ax == axes[0]:
            ax.legend(fontsize=8.5, loc="upper right",
                      handles=[b_sens, b_res, b_neut])

    fig.suptitle(
        "Nucleosome sensitivity at gene promoters by subcompartment\n"
        "Significance (* p<0.05, ** p<0.01, *** p<0.001) vs genome-wide average",
        fontsize=11, y=1.01
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")
def plot_stacked_bars(counts_dict: dict, title: str, out_path: Path):
    """
    counts_dict: {group_label: {subcomp: {sens: count}}}
    """
    n_groups = len(counts_dict)
    fig, axes = plt.subplots(1, n_groups,
                             figsize=(5 * n_groups, 5),
                             sharey=True)
    if n_groups == 1:
        axes = [axes]

    for ax, (group, counts) in zip(axes, counts_dict.items()):
        x      = np.arange(len(SUBCOMPS))
        totals = [sum(counts[sc].values()) for sc in SUBCOMPS]

        bottoms = np.zeros(len(SUBCOMPS))
        for sens in SENS_CALLS:
            proportions = np.array([
                100 * counts[sc][sens] / max(totals[i], 1)
                for i, sc in enumerate(SUBCOMPS)
            ])
            ax.bar(x, proportions, bottom=bottoms,
                   color=SENS_COLORS[sens], label=sens,
                   edgecolor="white", linewidth=0.3)
            bottoms += proportions

        ax.set_xticks(x)
        ax.set_xticklabels(SUBCOMPS, fontsize=9, rotation=0)
        ax.set_title(group, fontsize=11, fontweight="bold")
        ax.set_xlabel("Subcompartment", fontsize=9)
        ax.axvline(3.5, color="black", lw=1.0, ls="--", alpha=0.5)

        # Annotate n per subcompartment
        for i, n in enumerate(totals):
            ax.text(i, 102, f"n={n//1000:.0f}k" if n >= 1000 else f"n={n}",
                    ha="center", va="bottom", fontsize=6.5, color="grey")

        if ax == axes[0]:
            ax.set_ylabel("Nucleosomes (%)", fontsize=10)
            ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(title, fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 2: Sensitive fraction across subcompartments (line plot)
# ---------------------------------------------------------------------------
def plot_sensitive_fraction(counts_dict: dict, out_path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    from scipy.stats import spearmanr
    sc_rank = {sc: i for i, sc in enumerate(SUBCOMPS)}

    x = np.arange(len(SUBCOMPS))

    # Panel 1: Sensitive %
    ax = axes[0]
    for group, counts in counts_dict.items():
        col = "#e74c3c" if "Scr" in group else (
              "#2980b9" if "HP1a" in group else "#27ae60")
        sens_pct = np.array([
            100 * counts[sc]["Sensitive"] / max(sum(counts[sc].values()), 1)
            for sc in SUBCOMPS
        ])
        # Spearman r across subcompartment ranks
        all_vals, all_ranks = [], []
        for sc in SUBCOMPS:
            v = counts[sc]["Sensitive"]
            t = max(sum(counts[sc].values()), 1)
            all_vals.append(100 * v / t)
            all_ranks.append(sc_rank[sc])
        sr, _ = spearmanr(all_ranks, all_vals)
        ax.plot(x, sens_pct, marker="o", lw=2, ms=5,
                color=col, alpha=0.8, label=f"{group} (ρ={sr:.3f})")
    ax.set_xticks(x)
    ax.set_xticklabels(SUBCOMPS, fontsize=9)
    ax.set_ylabel("% Sensitive nucleosomes", fontsize=10)
    ax.set_title("Sensitive fraction per subcompartment", fontsize=11)
    ax.axvline(3.5, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel 2: Resistant %
    ax = axes[1]
    for group, counts in counts_dict.items():
        col = "#e74c3c" if "Scr" in group else (
              "#2980b9" if "HP1a" in group else "#27ae60")
        res_pct = np.array([
            100 * counts[sc]["Resistant"] / max(sum(counts[sc].values()), 1)
            for sc in SUBCOMPS
        ])
        all_vals2 = [100 * counts[sc]["Resistant"] /
                     max(sum(counts[sc].values()), 1) for sc in SUBCOMPS]
        sr2, _ = spearmanr(list(sc_rank.values()), all_vals2)
        ax.plot(x, res_pct, marker="o", lw=2, ms=5,
                color=col, alpha=0.8, label=f"{group} (ρ={sr2:.3f})")
    ax.set_xticks(x)
    ax.set_xticklabels(SUBCOMPS, fontsize=9)
    ax.set_ylabel("% Resistant nucleosomes", fontsize=10)
    ax.set_title("Resistant fraction per subcompartment", fontsize=11)
    ax.axvline(3.5, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("Nucleosome sensitivity at gene promoters by subcompartment",
                 fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 3: Heatmap of z-scored proportions
# ---------------------------------------------------------------------------
def plot_heatmap(counts_dict: dict, out_path: Path):
    groups  = list(counts_dict.keys())
    n_g     = len(groups)
    sens_order = SENS_CALLS

    fig, axes = plt.subplots(1, len(sens_order),
                             figsize=(5 * len(sens_order), max(3, n_g * 0.6 + 1.5)))

    for ax, sens in zip(axes, sens_order):
        mat = np.zeros((n_g, len(SUBCOMPS)))
        for i, group in enumerate(groups):
            counts = counts_dict[group]
            for j, sc in enumerate(SUBCOMPS):
                total = sum(counts[sc].values())
                mat[i, j] = 100 * counts[sc][sens] / max(total, 1)

        # Z-score across subcompartments per group
        for i in range(n_g):
            mu, sd = mat[i].mean(), mat[i].std()
            if sd > 0:
                mat[i] = (mat[i] - mu) / sd

        vmax = max(1.5, np.nanpercentile(np.abs(mat), 95))
        im = ax.imshow(mat, cmap="RdBu_r" if sens == "Sensitive" else "RdBu",
                       vmin=-vmax, vmax=vmax, aspect="auto")

        ax.set_xticks(range(len(SUBCOMPS)))
        ax.set_xticklabels(SUBCOMPS, fontsize=8)
        ax.set_yticks(range(n_g))
        ax.set_yticklabels(groups if ax == axes[0] else [], fontsize=9)
        ax.set_title(f"% {sens}\n(z-scored)", fontsize=10,
                     color=SENS_COLORS[sens], fontweight="bold")
        ax.axvline(3.5, color="black", lw=1.5)

        for i in range(n_g):
            for j in range(len(SUBCOMPS)):
                tc = "white" if abs(mat[i, j]) > vmax * 0.6 else "black"
                ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center",
                        fontsize=7, color=tc)

        plt.colorbar(im, ax=ax, shrink=0.6, label="Z-score")

    fig.suptitle("Nucleosome sensitivity proportions by subcompartment (z-scored)",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Gene-level promoter classification
# ---------------------------------------------------------------------------
PROMOTER_CLASSES = [
    "Sens/Sens",   # sensitive upstream AND downstream
    "Sens/Res",    # sensitive upstream, resistant downstream (active NDR pattern)
    "Res/Res",     # resistant upstream AND downstream
    "Res/Sens",    # resistant upstream, sensitive downstream
]

PCLASS_COLORS = {
    "Sens/Sens": "#27ae60",   # green
    "Sens/Res":  "#D4A017",   # gold
    "Res/Res":   "#1B3A6B",   # navy
    "Res/Sens":  "#8e44ad",   # purple
}


def sr_score(group: list) -> float:
    """
    Compute (S-R)/(S+R) for a list of sensitivity calls, ignoring Neutral.
    Returns NaN if no S or R calls present.
    Range: -1 (all resistant) to +1 (all sensitive).
    """
    s = sum(1 for r in group if r["sensitivity"] == "Sensitive")
    r = sum(1 for r in group if r["sensitivity"] == "Resistant")
    if s + r == 0:
        return np.nan
    return (s - r) / (s + r)


def classify_promoters(records: list,
                       min_nucs: int = 5,
                       threshold: float = 0.0) -> dict:
    """
    Classify each gene's promoter using continuous S/R score (ignoring Neutral).

    Score = (Sensitive - Resistant) / (Sensitive + Resistant) per region.
    Positive score = more sensitive, negative = more resistant.

    Upstream = closer to TSS (strand-aware):
      + strand: lower genomic coordinates
      - strand: higher genomic coordinates (reversed)

    Classification (threshold=0 by default):
      Sens/Sens: upstream > threshold AND downstream > threshold
      Sens/Res:  upstream > threshold AND downstream < -threshold
      Res/Res:   upstream < -threshold AND downstream < -threshold
      Res/Sens:  upstream < -threshold AND downstream > threshold

    Genes where upstream OR downstream has no S/R calls are excluded.

    Parameters
    ----------
    min_nucs  : minimum total nucleosomes per gene
    threshold : score threshold (0 = any imbalance counts)
    """
    from collections import defaultdict
    gene_records = defaultdict(list)
    for r in records:
        key = (r["gene_name"], r["tss_pos"])
        gene_records[key].append(r)

    gene_classes = {}
    for key, recs in gene_records.items():
        if len(recs) < min_nucs:
            continue

        # Consistent subcompartment
        sc_counts = defaultdict(int)
        for r in recs:
            sc_counts[r["subcompartment"]] += 1
        sc = max(sc_counts, key=sc_counts.get)

        strand = recs[0].get("strand", "+")

        # Sort by nucleosome midpoint, orient 5'→3'
        recs_sorted = sorted(recs,
                             key=lambda r: (r["nuc_start"] + r["nuc_end"]) / 2)
        if strand == "-":
            recs_sorted = recs_sorted[::-1]

        n = len(recs_sorted)
        t = max(1, n // 3)
        up   = recs_sorted[:t]
        down = recs_sorted[-t:]

        us = sr_score(up)
        ds = sr_score(down)

        # Skip if no scoreable S/R calls in either region
        if np.isnan(us) or np.isnan(ds):
            continue

        # Classify
        if us > threshold and ds > threshold:
            cls = "Sens/Sens"
        elif us > threshold and ds < -threshold:
            cls = "Sens/Res"
        elif us < -threshold and ds < -threshold:
            cls = "Res/Res"
        elif us < -threshold and ds > threshold:
            cls = "Res/Sens"
        else:
            # Ambiguous — only occurs when threshold > 0
            continue

        gene_classes[key] = {
            "class":          cls,
            "subcompartment": sc,
            "n_nucs":         n,
            "strand":         strand,
            "gene_name":      recs[0]["gene_name"],
            "up_score":       us,
            "down_score":     ds,
        }

    return gene_classes


def count_promoter_classes(gene_classes: dict) -> dict:
    """Returns dict: subcomp -> {class: count}"""
    from collections import defaultdict
    counts = {sc: defaultdict(int) for sc in SUBCOMPS}
    for key, info in gene_classes.items():
        sc  = info["subcompartment"]
        cls = info["class"]
        if sc in counts:
            counts[sc][cls] += 1
    return counts


def plot_promoter_classes(class_counts_by_group: dict, out_path: Path):
    """
    Grouped bar plot of promoter class counts per subcompartment.
    One panel per group (Scr, HP1a, HP1b).
    Bars: Sensitive (green), Mixed (orange), Resistant (navy), Other (grey).
    Significance stars from chi-squared vs genome average.
    """
    from scipy.stats import chi2_contingency

    CLASSES_PLOT = ["Sens/Sens", "Sens/Res", "Res/Res", "Res/Sens"]
    BAR_W  = 0.18
    OFFSET = np.array([-1.5, -0.5, 0.5, 1.5]) * BAR_W
    x      = np.arange(len(SUBCOMPS))

    n_groups = len(class_counts_by_group)
    fig, axes = plt.subplots(1, n_groups,
                             figsize=(6 * n_groups, 5.5),
                             sharey=False)
    if n_groups == 1:
        axes = [axes]

    for ax, (group, counts) in zip(axes, class_counts_by_group.items()):

        # Genome-wide totals
        grand = {cls: sum(counts[sc][cls] for sc in SUBCOMPS)
                 for cls in CLASSES_PLOT}
        grand_total = sum(grand.values())

        bars_drawn = []
        for k, cls in enumerate(CLASSES_PLOT):
            vals = np.array([counts[sc][cls] for sc in SUBCOMPS], float)
            col  = PCLASS_COLORS[cls]
            b    = ax.bar(x + OFFSET[k], vals, BAR_W,
                          color=col, label=cls.replace("\n", " "),
                          edgecolor="white", linewidth=0.4, zorder=3)
            bars_drawn.append((cls, vals))

        # Significance stars above each subcompartment
        ymax = max(
            max(counts[sc][cls] for cls in CLASSES_PLOT)
            for sc in SUBCOMPS
        ) * 1.05
        for i, sc in enumerate(SUBCOMPS):
            sc_total = sum(counts[sc][cls] for cls in CLASSES_PLOT)
            if sc_total < 5:
                continue
            rest = {cls: grand[cls] - counts[sc][cls] for cls in CLASSES_PLOT}
            table = [
                [counts[sc][cls] for cls in CLASSES_PLOT],
                [rest[cls]       for cls in CLASSES_PLOT],
            ]
            try:
                _, p, _, _ = chi2_contingency(table)
            except Exception:
                p = 1.0
            sig = "***" if p < 0.001 else ("**" if p < 0.01 else
                  ("*" if p < 0.05 else ""))
            if sig:
                ax.text(x[i], ymax * 1.03, sig,
                        ha="center", va="bottom",
                        fontsize=9, color="#2c3e50", fontweight="bold")

        ax.axvline(3.5, color="black", lw=1.0, ls="--", alpha=0.5)
        ax.text(1.5, ymax * 1.14, "← A compartment →",
                ha="center", fontsize=8, color="#922b21")
        ax.text(5.5, ymax * 1.14, "← B compartment →",
                ha="center", fontsize=8, color="#1a5276")

        ax.set_xticks(x)
        ax.set_xticklabels(SUBCOMPS, fontsize=9)
        ax.set_xlabel("Subcompartment\n(A3 = most active → B3 = most inactive)",
                      fontsize=9)
        ax.set_ylabel("Gene count", fontsize=10)
        ax.set_title(group, fontsize=12, fontweight="bold")
        ax.set_ylim(0, ymax * 1.2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if ax == axes[0]:
            handles = [mpatches.Patch(color=PCLASS_COLORS[c],
                                      label=c.replace("\n", " "))
                       for c in CLASSES_PLOT]
            ax.legend(handles=handles, fontsize=8.5, loc="upper right")

    fig.suptitle(
        "Promoter sensitivity class by subcompartment (gene-level)\n"
        "Significance (* p<0.05, ** p<0.01, *** p<0.001) vs genome-wide average",
        fontsize=11, y=1.01
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bedgraph", default="MCF10A_subcompartments_20000_T2T.bg")
    p.add_argument("--files",    nargs="+",
                   default=sorted(glob.glob("TSS_10A_*txt")))
    p.add_argument("--out_dir",  default="nucleosome_results/")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading subcompartments ...")
    comp = load_subcompartments(args.bedgraph)

    # Process each file
    all_records  = {}   # fname -> records
    group_counts = {}   # group_label -> merged counts

    for fpath in sorted(args.files):
        fname = Path(fpath).name
        print(f"Loading {fname} ...")
        records = load_sensitivity_file(fpath, comp)
        all_records[fname] = records
        counts = count_by_subcomp(records)

        # Determine group
        group = "Other"
        for g, files in SAMPLE_GROUPS.items():
            if any(fname == Path(f).name for f in files):
                group = g
                break

        # Merge into group counts
        sample_label = fname.replace("TSS_10A_", "").replace(
            "_sensitivity_overview.txt", "")
        group_counts[sample_label] = counts

        print(f"  {len(records):,} nucleosomes assigned to subcompartments")
        print(f"  Sensitivity breakdown:", end="")
        for s in SENS_CALLS:
            n = sum(counts[sc][s] for sc in SUBCOMPS)
            print(f"  {s}={n:,}", end="")
        print()

    # Merged counts per group
    merged_counts = {}
    for group_name, files in SAMPLE_GROUPS.items():
        mc = {sc: {s: 0 for s in SENS_CALLS} for sc in SUBCOMPS}
        for fpath in files:
            fname = Path(fpath).name
            if fname in all_records:
                c = count_by_subcomp(all_records[fname])
                for sc in SUBCOMPS:
                    for s in SENS_CALLS:
                        mc[sc][s] += c[sc][s]
        merged_counts[group_name] = mc

    # Statistical tests on merged data
    print(f"\n{'='*70}")
    print("Chi-squared test: sensitivity distribution vs genome average")
    print(f"{'='*70}")
    print(f"\n{'Group':<10} {'Subcomp':<8} {'%Sens':>7} {'%Res':>7} "
          f"{'chi2':>8} {'p':>12} {'sig':>5} {'n':>8}")
    print("-" * 65)

    stats_out = []
    for group_name, counts in merged_counts.items():
        stats = test_enrichment(counts)
        for sc in SUBCOMPS:
            if sc not in stats:
                continue
            chi2, p, sp, rp, n = stats[sc]
            sig = "***" if p < 0.001 else ("**" if p < 0.01 else
                  ("*" if p < 0.05 else "ns"))
            print(f"{group_name:<10} {sc:<8} {sp:>7.1f} {rp:>7.1f} "
                  f"{chi2:>8.1f} {p:>12.2e} {sig:>5} {n:>8,}")
            stats_out.append((group_name, sc, sp, rp, chi2, p, sig, n))
        print()

    # Save stats
    with open(out_dir / "sensitivity_stats.txt", "w") as f:
        f.write("Group\tSubcompartment\t%Sensitive\t%Resistant\t"
                "chi2\tp-value\tsig\tn\n")
        for row in stats_out:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"  Stats saved: {out_dir/'sensitivity_stats.txt'}")

    # Plots
    print("\nGenerating plots ...")
    plot_grouped_bars(merged_counts,
                      out_dir / "sensitivity_grouped_bars.png")
    plot_stacked_bars(merged_counts,
                      "Nucleosome sensitivity at gene promoters by subcompartment",
                      out_dir / "sensitivity_stacked_bars.png")
    plot_sensitive_fraction(merged_counts,
                            out_dir / "sensitivity_fractions.png")
    plot_heatmap(group_counts,
                 out_dir / "sensitivity_heatmap.png")

    # Gene-level promoter classification
    print("\nClassifying promoters by positional sensitivity pattern ...")
    gene_class_by_group = {}
    for group_name, files in SAMPLE_GROUPS.items():
        # Merge records from both clones in group
        merged_records = []
        for fpath in files:
            fname = Path(fpath).name
            if fname in all_records:
                merged_records.extend(all_records[fname])
        gene_classes = classify_promoters(merged_records)
        gene_class_by_group[group_name] = gene_classes
        # Count classes per subcompartment
        class_counts = count_promoter_classes(gene_classes)
        total_genes  = sum(sum(class_counts[sc].values())
                           for sc in SUBCOMPS)
        print(f"  {group_name}: {total_genes:,} genes classified")
        class_totals = defaultdict(int)
        for sc in SUBCOMPS:
            for cls, n in class_counts[sc].items():
                class_totals[cls] += n
        for cls, n in sorted(class_totals.items(),
                             key=lambda x: -x[1]):
            print(f"    {cls}: {n:,} ({100*n/total_genes:.1f}%)")

    # Plot gene-level grouped bars
    gene_class_counts = {g: count_promoter_classes(gc)
                         for g, gc in gene_class_by_group.items()}
    plot_promoter_classes(gene_class_counts,
                          out_dir / "promoter_classes_grouped_bars.png")

    print(f"\n✓  Done. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
