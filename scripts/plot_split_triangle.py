#!/usr/bin/env python3
"""
plot_split_triangle.py
=============================================================================
Plots two contact matrices as a split triangle figure:
  - Upper triangle: track A (e.g. condense-seq / 20U supernatant rep1)
  - Lower triangle: track B (e.g. MNase sensitivity / rep2)
  - Diagonal: shared genomic axis
  - Both raw contact and O/E versions produced

Usage:
    python plot_split_triangle.py \
        --upper_contact  results_supernatant_rep1_chr2_25kb_gpu/contact_matrix_avg.npy \
        --upper_phi      data/phi_supernatant_rep1_chr2_4mb_25kb.npz \
        --upper_label    "Condense-seq (20U supernatant rep1)" \
        --lower_contact  results_rep2_chr2_25kb_gpu/contact_matrix_avg.npy \
        --lower_phi      data/phi_rep2_chr2_4mb_25kb.npz \
        --lower_label    "Genomic rep2 sensitivity" \
        --out_dir        results_summary/
=============================================================================
"""

import argparse
import numpy as np
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import Normalize
    from matplotlib import colormaps
except ImportError:
    raise ImportError("pip install matplotlib")


# ---------------------------------------------------------------------------
# O/E matrix
# ---------------------------------------------------------------------------
def observed_over_expected(contact: np.ndarray) -> np.ndarray:
    n  = contact.shape[0]
    oe = np.zeros_like(contact, dtype=float)
    for d in range(1, n):
        diag = np.array([contact[i, i+d] for i in range(n-d)])
        mean_d = diag.mean()
        norm   = diag / mean_d if mean_d > 0 else np.ones_like(diag)
        for i in range(n-d):
            oe[i, i+d] = norm[i]
            oe[i+d, i] = norm[i]
    np.fill_diagonal(oe, np.nan)
    return oe


# ---------------------------------------------------------------------------
# Core split triangle plot
# ---------------------------------------------------------------------------
def plot_split_triangle(upper_mat: np.ndarray,
                        lower_mat: np.ndarray,
                        bin_starts: np.ndarray,
                        binsize: int,
                        upper_label: str,
                        lower_label: str,
                        out_path: Path,
                        cmap: str = "Reds",
                        vmax: float = None,
                        colorbar_label: str = "Contact probability",
                        title: str = ""):
    """
    Split triangle contact map.
    Upper triangle = upper_mat, lower triangle = lower_mat.
    Both are shown on the same colour scale.
    """
    n  = upper_mat.shape[0]
    mb = bin_starts / 1e6

    # Build combined display matrix
    combined = np.full((n, n), np.nan)

    # Upper triangle (row < col) from upper_mat
    for i in range(n):
        for j in range(i+1, n):
            combined[i, j] = upper_mat[i, j]

    # Lower triangle (row > col) from lower_mat
    for i in range(n):
        for j in range(i):
            combined[i, j] = lower_mat[i, j]

    # Shared colour scale
    if vmax is None:
        upper_vals = upper_mat[np.triu_indices(n, k=1)]
        lower_vals = lower_mat[np.tril_indices(n, k=-1)]
        all_vals   = np.concatenate([upper_vals, lower_vals])
        all_vals   = all_vals[np.isfinite(all_vals) & (all_vals > 0)]
        vmax       = np.percentile(all_vals, 99) if len(all_vals) > 0 else 1.0

    fig, ax = plt.subplots(figsize=(8, 7.5))

    if cmap == "RdBu_r":
        im = ax.imshow(combined, cmap=cmap, origin="upper",
                       vmin=-vmax, vmax=vmax,
                       interpolation="nearest", aspect="equal")
    else:
        im = ax.imshow(combined, cmap=cmap, origin="upper",
                       vmin=0, vmax=vmax,
                       interpolation="nearest", aspect="equal")

    # Diagonal line
    ax.plot([0, n-1], [0, n-1], color="#333333", lw=0.8, zorder=5)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.55, pad=0.02)
    cbar.set_label(colorbar_label, fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    # Genomic coordinate ticks
    mb_start = mb[0]
    mb_end   = mb[-1] + binsize / 1e6
    tick_mb  = np.linspace(mb_start, mb_end, 7)
    tick_idx = [(np.abs(mb - t)).argmin() for t in tick_mb]
    tick_idx = sorted(set(np.clip(tick_idx, 0, n-1)))

    ax.set_xticks(tick_idx)
    ax.set_yticks(tick_idx)
    ax.set_xticklabels([f"{mb[i]:.1f}" for i in tick_idx],
                       fontsize=9, rotation=45, ha="right")
    ax.set_yticklabels([f"{mb[i]:.1f}" for i in tick_idx], fontsize=9)
    ax.set_xlabel("Genomic position (Mb)", fontsize=11)
    ax.set_ylabel("Genomic position (Mb)", fontsize=11)

    # Labels for each triangle
    # Upper right label
    ax.text(0.97, 0.03, upper_label,
            transform=ax.transAxes, fontsize=10, fontweight="bold",
            ha="right", va="bottom", color="#c0392b",
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec="#c0392b", alpha=0.85))
    # Lower left label
    ax.text(0.03, 0.97, lower_label,
            transform=ax.transAxes, fontsize=10, fontweight="bold",
            ha="left", va="top", color="#2980b9",
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec="#2980b9", alpha=0.85))

    if title:
        ax.set_title(title, fontsize=12, pad=10)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--upper_contact", required=True)
    p.add_argument("--upper_phi",     required=True)
    p.add_argument("--upper_label",   default="Condense-seq (20U supernatant rep1)")
    p.add_argument("--lower_contact", required=True)
    p.add_argument("--lower_phi",     required=True)
    p.add_argument("--lower_label",   default="Genomic rep2 sensitivity")
    p.add_argument("--out_dir",       default="results_summary/")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load matrices
    print("Loading contact matrices ...")
    upper_contact = np.load(args.upper_contact).astype(float)
    lower_contact = np.load(args.lower_contact).astype(float)

    assert upper_contact.shape == lower_contact.shape, \
        f"Matrix shapes must match: {upper_contact.shape} vs {lower_contact.shape}"

    # Load genomic coordinates from phi file
    phi_data   = np.load(args.upper_phi, allow_pickle=True)
    bin_starts = phi_data["bin_starts"]
    binsize    = int(phi_data["binsize"].flat[0])
    chrom      = str(phi_data["chrom"][0])
    r_start    = int(phi_data["region_start"].flat[0])
    r_end      = int(phi_data["region_end"].flat[0])

    region_str = f"{chrom}:{r_start//1_000_000:.0f}–{r_end//1_000_000:.0f} Mb"
    print(f"Region: {region_str}  ({upper_contact.shape[0]} beads @ {binsize//1000} kb)")

    # Symmetrise
    def symmetrise(m):
        return np.triu(m) + np.triu(m, 1).T

    upper_contact = symmetrise(upper_contact)
    lower_contact = symmetrise(lower_contact)

    # --- Plot 1: Raw contact probability ---
    print("\nPlotting raw contact split triangle ...")
    plot_split_triangle(
        upper_mat     = upper_contact,
        lower_mat     = lower_contact,
        bin_starts    = bin_starts,
        binsize       = binsize,
        upper_label   = args.upper_label,
        lower_label   = args.lower_label,
        out_path      = out_dir / "split_triangle_contact.png",
        cmap          = "Reds",
        colorbar_label= "Contact probability",
        title         = f"Contact probability — {region_str}",
    )

    # --- Plot 2: O/E ---
    print("Computing O/E matrices ...")
    upper_oe = observed_over_expected(upper_contact)
    lower_oe = observed_over_expected(lower_contact)

    # Shared vmax for OE
    u_vals  = upper_oe[np.isfinite(upper_oe)]
    l_vals  = lower_oe[np.isfinite(lower_oe)]
    vmax_oe = np.percentile(np.abs(np.concatenate([u_vals, l_vals])), 97)

    print("Plotting O/E split triangle ...")
    plot_split_triangle(
        upper_mat     = upper_oe,
        lower_mat     = lower_oe,
        bin_starts    = bin_starts,
        binsize       = binsize,
        upper_label   = args.upper_label,
        lower_label   = args.lower_label,
        out_path      = out_dir / "split_triangle_OE.png",
        cmap          = "RdBu_r",
        vmax          = vmax_oe,
        colorbar_label= "O/E ratio",
        title         = f"O/E contact map — {region_str}",
    )

    print(f"\nDone. Outputs in {out_dir}")
    print(f"   split_triangle_contact.png")
    print(f"   split_triangle_OE.png")


if __name__ == "__main__":
    main()
