#!/usr/bin/env python3
"""
promoter_class_proportions.py
=============================================================================
Plots promoter sensitivity class proportions (%) per subcompartment,
normalised within each subcompartment so bars sum to 100%.
Companion to nucleosome_sensitivity_vs_compartment.py — reuses the same
classifier logic but focuses on relative enrichment rather than raw counts.

Usage:
    python promoter_class_proportions.py \
        --bedgraph MCF10A_subcompartments_20000_T2T.bg \
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
except ImportError:
    raise ImportError("pip install matplotlib")

try:
    from scipy.stats import chi2_contingency
except ImportError:
    raise ImportError("pip install scipy")


# ---------------------------------------------------------------------------
# Constants — must match nucleosome_sensitivity_vs_compartment.py
# ---------------------------------------------------------------------------
SUBCOMPS = ["A3", "A2", "A1", "A0", "B0", "B1", "B2", "B3"]

PROMOTER_CLASSES = ["Sens/Sens", "Sens/Res", "Res/Res", "Res/Sens"]

# Colorblind-friendly palette (Wong 2011)
PCLASS_COLORS = {
    "Sens/Sens": "#009E73",   # teal green
    "Sens/Res":  "#E69F00",   # orange/gold
    "Res/Res":   "#0072B2",   # blue
    "Res/Sens":  "#CC79A7",   # pink/mauve
}

# Colorblind-friendly group colors
GROUP_COLORS = {
    "Scr":  "#000000",   # black
    "HP1a": "#E69F00",   # orange
    "HP1b": "#56B4E9",   # sky blue
}

SAMPLE_GROUPS = {
    "Scr":  ["TSS_10A_Scr_clone1_sensitivity_overview.txt",
             "TSS_10A_Scr_clone2_sensitivity_overview.txt"],
    "HP1a": ["TSS_10A_HP1a_clone5_sensitivity_overview.txt",
             "TSS_10A_HP1a_clone6_sensitivity_overview.txt"],
    "HP1b": ["TSS_10A_HP1b_clone15_sensitivity_overview.txt",
             "TSS_10A_HP1b_cloneD_sensitivity_overview.txt"],
}


# ---------------------------------------------------------------------------
# Subcompartment lookup
# ---------------------------------------------------------------------------
def load_subcompartments(bg_path: str) -> dict:
    comp = defaultdict(list)
    with open(bg_path) as f:
        for line in f:
            p = line.strip().split()
            if len(p) < 4:
                continue
            comp[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for chrom in comp:
        comp[chrom].sort()
    return dict(comp)


def query_subcompartment(comp: dict, chrom: str, pos: int) -> str:
    intervals = comp.get(chrom, [])
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
# Load sensitivity file
# ---------------------------------------------------------------------------
def load_sensitivity_file(fpath: str, comp: dict) -> list:
    records = []
    with open(fpath) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 10:
                p = line.strip().split()
            if len(p) < 10:
                continue
            sens   = p[3]
            strand = p[9]
            if sens not in ("Sensitive", "Neutral", "Resistant"):
                continue
            gene_chr   = p[4]
            gene_start = int(p[5])
            gene_end   = int(p[6])
            tss = gene_start if strand == "+" else gene_end
            sc  = query_subcompartment(comp, gene_chr, tss)
            if sc is None:
                continue
            records.append({
                "nuc_start":     int(p[1]),
                "nuc_end":       int(p[2]),
                "sensitivity":   sens,
                "gene_name":     p[7],
                "tss_pos":       tss,
                "subcompartment": sc,
                "strand":        strand,
            })
    return records


# ---------------------------------------------------------------------------
# Classifier (continuous S/R score, ignoring Neutral)
# ---------------------------------------------------------------------------
def sr_score(group: list) -> float:
    s = sum(1 for r in group if r["sensitivity"] == "Sensitive")
    r = sum(1 for r in group if r["sensitivity"] == "Resistant")
    if s + r == 0:
        return np.nan
    return (s - r) / (s + r)


def classify_promoters(records: list,
                       min_nucs: int = 5,
                       threshold: float = 0.0) -> dict:
    gene_records = defaultdict(list)
    for r in records:
        gene_records[(r["gene_name"], r["tss_pos"])].append(r)

    gene_classes = {}
    for key, recs in gene_records.items():
        if len(recs) < min_nucs:
            continue
        sc_counts = defaultdict(int)
        for r in recs:
            sc_counts[r["subcompartment"]] += 1
        sc     = max(sc_counts, key=sc_counts.get)
        strand = recs[0]["strand"]

        recs_sorted = sorted(recs,
                             key=lambda r: (r["nuc_start"] + r["nuc_end"]) / 2)
        if strand == "-":
            recs_sorted = recs_sorted[::-1]

        n  = len(recs_sorted)
        t  = max(1, n // 3)
        us = sr_score(recs_sorted[:t])
        ds = sr_score(recs_sorted[-t:])

        if np.isnan(us) or np.isnan(ds):
            continue

        if us > threshold and ds > threshold:
            cls = "Sens/Sens"
        elif us > threshold and ds < -threshold:
            cls = "Sens/Res"
        elif us < -threshold and ds < -threshold:
            cls = "Res/Res"
        elif us < -threshold and ds > threshold:
            cls = "Res/Sens"
        else:
            continue

        gene_classes[key] = {
            "class": cls, "subcompartment": sc,
            "n_nucs": n, "strand": strand,
            "gene_name": recs[0]["gene_name"],
        }
    return gene_classes


def count_promoter_classes(gene_classes: dict) -> dict:
    counts = {sc: defaultdict(int) for sc in SUBCOMPS}
    for key, info in gene_classes.items():
        sc  = info["subcompartment"]
        cls = info["class"]
        if sc in counts:
            counts[sc][cls] += 1
    return counts


# ---------------------------------------------------------------------------
# Plot 1: Proportional stacked bars (100%) per subcompartment
# ---------------------------------------------------------------------------
def plot_proportional_stacked(class_counts_by_group: dict, out_path: Path):
    """
    Each bar = 100%, coloured by promoter class proportion.
    Colorblind-friendly palette. Legend outside to the right.
    Significance stars above bars. n annotated below x-axis.
    """
    n_groups = len(class_counts_by_group)
    fig, axes = plt.subplots(1, n_groups,
                             figsize=(4.5 * n_groups + 2.0, 5.5),
                             sharey=True)
    if n_groups == 1:
        axes = [axes]

    for ax_idx, (ax, (group, counts)) in enumerate(
            zip(axes, class_counts_by_group.items())):

        grand  = {cls: sum(counts[sc][cls] for sc in SUBCOMPS)
                  for cls in PROMOTER_CLASSES}
        x      = np.arange(len(SUBCOMPS))
        totals = np.array([sum(counts[sc][cls] for cls in PROMOTER_CLASSES)
                           for sc in SUBCOMPS], float)
        bottoms = np.zeros(len(SUBCOMPS))

        for cls in PROMOTER_CLASSES:
            vals = np.array([counts[sc][cls] for sc in SUBCOMPS], float)
            pct  = np.where(totals > 0, 100 * vals / totals, 0)
            bars = ax.bar(x, pct, bottom=bottoms,
                          color=PCLASS_COLORS[cls],
                          label=cls,
                          edgecolor="white", linewidth=0.5, zorder=3)
            # Label inside bar segment if wide enough
            for i, (p, b) in enumerate(zip(pct, bottoms)):
                if p > 8:
                    ax.text(x[i], b + p / 2, f"{p:.0f}%",
                            ha="center", va="center",
                            fontsize=6.5, color="white", fontweight="bold")
            bottoms += pct

        # Significance stars above bars
        for i, sc in enumerate(SUBCOMPS):
            sc_total = totals[i]
            if sc_total < 5:
                continue
            rest  = {cls: grand[cls] - counts[sc][cls] for cls in PROMOTER_CLASSES}
            table = [[counts[sc][cls] for cls in PROMOTER_CLASSES],
                     [rest[cls]       for cls in PROMOTER_CLASSES]]
            try:
                _, p, _, _ = chi2_contingency(table)
            except Exception:
                p = 1.0
            sig = "***" if p < 0.001 else ("**" if p < 0.01 else
                  ("*" if p < 0.05 else ""))
            if sig:
                ax.text(x[i], 101.5, sig, ha="center", va="bottom",
                        fontsize=8.5, color="#333333", fontweight="bold")

        # A|B boundary
        ax.axvline(3.5, color="#333333", lw=1.2, ls="--", alpha=0.7)

        ax.set_xticks(x)
        ax.set_xticklabels(SUBCOMPS, fontsize=9.5)
        ax.set_ylim(0, 112)
        ax.set_title(group, fontsize=13, fontweight="bold", pad=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(ax_idx == 0)
        ax.tick_params(left=(ax_idx == 0))

        if ax_idx == 0:
            ax.set_ylabel("Genes (%)", fontsize=11)

        # n per subcompartment below x-axis — outside plot area
        for i, n in enumerate(totals):
            ax.text(x[i], -7, f"n={int(n):,}",
                    ha="center", va="top", fontsize=6.5,
                    color="#666666", rotation=35)

        # A / B labels
        ax.text(1.5, 108, "A compartment", ha="center", fontsize=8,
                color="#922b21", style="italic")
        ax.text(5.5, 108, "B compartment", ha="center", fontsize=8,
                color="#1a5276", style="italic")

    # Shared legend to the right of all panels
    handles = [mpatches.Patch(color=PCLASS_COLORS[c], label=c)
               for c in PROMOTER_CLASSES]
    fig.legend(handles=handles, fontsize=9.5, title="Promoter class",
               title_fontsize=10, loc="center right",
               bbox_to_anchor=(1.0, 0.5), framealpha=0.9,
               edgecolor="#cccccc")

    fig.suptitle(
        "Promoter sensitivity class proportions by subcompartment\n"
        "Significance (* p<0.05  ** p<0.01  *** p<0.001) vs genome-wide average",
        fontsize=11, y=1.02
    )
    fig.tight_layout(rect=[0, 0, 0.88, 1])
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Plot 2: Line plot of each class proportion across subcompartments
# ---------------------------------------------------------------------------
def plot_class_profiles(class_counts_by_group: dict, out_path: Path):
    """
    One panel per promoter class, one line per group.
    Y-axis = % of genes in that subcompartment with that class.
    """
    n_classes = len(PROMOTER_CLASSES)
    fig, axes = plt.subplots(1, n_classes,
                             figsize=(4.5 * n_classes, 4.5),
                             sharey=False)

    GROUP_COLORS_LOCAL = GROUP_COLORS  # use module-level colorblind palette
    x = np.arange(len(SUBCOMPS))

    for ax, cls in zip(axes, PROMOTER_CLASSES):
        for group, counts in class_counts_by_group.items():
            totals = np.array([sum(counts[sc][c] for c in PROMOTER_CLASSES)
                               for sc in SUBCOMPS], float)
            vals   = np.array([counts[sc][cls] for sc in SUBCOMPS], float)
            pct    = np.where(totals > 0, 100 * vals / totals, 0)
            col    = GROUP_COLORS_LOCAL.get(group, "grey")
            ls     = "-" if group == "Scr" else ("--" if group == "HP1a" else ":")
            ax.plot(x, pct, marker="o", lw=2, ms=5,
                    color=col, alpha=0.9, label=group, linestyle=ls)

        ax.set_xticks(x)
        ax.set_xticklabels(SUBCOMPS, fontsize=8.5)
        ax.set_title(cls, fontsize=11, fontweight="bold",
                     color=PCLASS_COLORS[cls])
        ax.set_ylabel("% of genes in subcompartment", fontsize=9)
        ax.axvline(3.5, color="black", lw=0.8, ls="--", alpha=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=8)

    fig.suptitle(
        "Promoter class % within each subcompartment across sample groups",
        fontsize=11, y=1.01
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Statistical test: A compartment vs B compartment per class per group
# ---------------------------------------------------------------------------
def test_ab_enrichment(class_counts_by_group: dict, out_path: Path):
    """
    For each group and promoter class, test whether the proportion differs
    between A compartment (A0-A3) and B compartment (B0-B3) using
    chi-squared. Reports odds ratio and p-value.
    """
    A_COMPS = ["A3", "A2", "A1", "A0"]
    B_COMPS = ["B0", "B1", "B2", "B3"]

    print(f"\n{'='*72}")
    print("A vs B compartment enrichment per promoter class")
    print(f"{'='*72}")
    print(f"\n{'Group':<8} {'Class':<12} {'A%':>6} {'B%':>6} "
          f"{'OR':>7} {'chi2':>8} {'p':>12} {'sig':>5}")
    print("-" * 65)

    rows = []
    for group, counts in class_counts_by_group.items():
        for cls in PROMOTER_CLASSES:
            a_cls   = sum(counts[sc][cls] for sc in A_COMPS)
            a_other = sum(counts[sc][c]   for sc in A_COMPS
                          for c in PROMOTER_CLASSES if c != cls)
            b_cls   = sum(counts[sc][cls] for sc in B_COMPS)
            b_other = sum(counts[sc][c]   for sc in B_COMPS
                          for c in PROMOTER_CLASSES if c != cls)

            a_total = a_cls + a_other
            b_total = b_cls + b_other
            a_pct   = 100 * a_cls / a_total if a_total > 0 else np.nan
            b_pct   = 100 * b_cls / b_total if b_total > 0 else np.nan

            # Odds ratio: (a_cls/a_other) / (b_cls/b_other)
            if a_other > 0 and b_cls > 0 and b_other > 0:
                OR = (a_cls / a_other) / (b_cls / b_other)
            else:
                OR = np.nan

            table = [[a_cls, a_other], [b_cls, b_other]]
            try:
                chi2, p, _, _ = chi2_contingency(table)
            except Exception:
                chi2, p = np.nan, 1.0

            sig = "***" if p < 0.001 else ("**" if p < 0.01 else
                  ("*" if p < 0.05 else "ns"))
            print(f"{group:<8} {cls:<12} {a_pct:>6.1f} {b_pct:>6.1f} "
                  f"{OR:>7.3f} {chi2:>8.1f} {p:>12.2e} {sig:>5}")
            rows.append((group, cls, a_pct, b_pct, OR, chi2, p, sig))
        print()

    # Save to file
    with open(out_path, "w") as f:
        f.write("Group\tClass\tA_pct\tB_pct\tOddsRatio\tchi2\tp_value\tsig\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"  Stats saved: {out_path}")
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bedgraph", default="MCF10A_subcompartments_20000_T2T.bg")
    p.add_argument("--out_dir",  default="nucleosome_results/")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading subcompartments ...")
    comp = load_subcompartments(args.bedgraph)

    # Load and classify all files
    all_records = {}
    for group_name, files in SAMPLE_GROUPS.items():
        merged = []
        for fpath in files:
            if not Path(fpath).exists():
                print(f"  WARNING: {fpath} not found, skipping")
                continue
            print(f"  Loading {Path(fpath).name} ...")
            merged.extend(load_sensitivity_file(fpath, comp))
        all_records[group_name] = merged

    # Classify promoters per group
    class_counts_by_group = {}
    print("\nClassifying promoters ...")
    for group_name, records in all_records.items():
        gene_classes = classify_promoters(records)
        counts       = count_promoter_classes(gene_classes)
        class_counts_by_group[group_name] = counts

        total = sum(sum(counts[sc].values()) for sc in SUBCOMPS)
        print(f"\n  {group_name}: {total:,} genes classified")

        # Print proportion table
        print(f"  {'Subcomp':<8}", end="")
        for cls in PROMOTER_CLASSES:
            print(f"  {cls:>10}", end="")
        print(f"  {'n':>6}")
        print("  " + "-" * (8 + 13 * len(PROMOTER_CLASSES) + 8))

        for sc in SUBCOMPS:
            sc_total = sum(counts[sc][cls] for cls in PROMOTER_CLASSES)
            print(f"  {sc:<8}", end="")
            for cls in PROMOTER_CLASSES:
                pct = 100 * counts[sc][cls] / sc_total if sc_total > 0 else 0
                print(f"  {pct:>9.1f}%", end="")
            print(f"  {sc_total:>6,}")

    print("\nGenerating plots ...")
    plot_proportional_stacked(class_counts_by_group,
                              out_dir / "promoter_class_proportions_stacked.png")
    plot_class_profiles(class_counts_by_group,
                        out_dir / "promoter_class_profiles.png")
    test_ab_enrichment(class_counts_by_group,
                       out_dir / "promoter_class_AB_stats.txt")

    print(f"\n✓  Done. Outputs in {out_dir}")


if __name__ == "__main__":
    main()
