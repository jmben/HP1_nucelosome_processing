#!/usr/bin/env python3
"""
04_analyse_results.py  —  Contact matrix + A/B compartment analysis
=============================================================================
Loads the averaged contact probability matrix from simulation and:

  1. Plots the contact probability map (log-scaled)
  2. Computes the observed/expected (O/E) matrix
  3. Runs PCA on the O/E correlation matrix to extract compartment PC1
     (positive = A compartment, negative = B compartment)
  4. Saves a compartment score bedGraph and summary plots

Usage:
    python 04_analyse_results.py \
        --contact  results/contact_matrix_avg.npy \
        --phi      data/phi_vector.npz \
        --out_dir  results/
=============================================================================
"""

import argparse
import numpy as np
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")   # headless
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
except ImportError:
    raise ImportError("pip install matplotlib")

try:
    from scipy.stats import pearsonr
    from scipy.linalg import eigh
except ImportError:
    raise ImportError("pip install scipy")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Analyse simulation contact matrix.")
    p.add_argument("--contact",  default="results/contact_matrix_avg.npy")
    p.add_argument("--phi",      default="data/phi_vector.npz")
    p.add_argument("--out_dir",  default="results/")
    p.add_argument("--contact_threshold", type=float, default=1e-6,
                   help="Set contacts below this value to NaN for log-plotting")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Helper: observed/expected matrix
# ---------------------------------------------------------------------------
def observed_over_expected(contact: np.ndarray) -> np.ndarray:
    """
    O/E normalisation: divide each diagonal (genomic distance) by its
    mean value, removing the distance-dependent decay common to all
    polymer contact maps.

    Returns O/E matrix (NaN on diagonal).
    """
    n   = contact.shape[0]
    oe  = np.full_like(contact, np.nan, dtype=float)

    for d in range(1, n):
        diag_vals = np.array([contact[i, i + d] for i in range(n - d)])
        mean_d    = diag_vals.mean()
        if mean_d > 0:
            norm = diag_vals / mean_d
        else:
            norm = np.ones_like(diag_vals)
        for i in range(n - d):
            oe[i, i + d] = norm[i]
            oe[i + d, i] = norm[i]

    return oe


# ---------------------------------------------------------------------------
# A/B compartments via PCA on O/E correlation matrix
# ---------------------------------------------------------------------------
def ab_compartments_pca(oe: np.ndarray) -> tuple:
    """
    Standard Hi-C compartment calling:
    1. Build Pearson correlation matrix of O/E
    2. Eigendecompose → first eigenvector = compartment score
    3. Sign convention: positive bins associated with higher condensability
       (active/A compartment in the paper's framework where high
       condensability → B compartment — see note below) are signed
       consistently with the phi vector.

    Note on sign:
        In the paper, HIGH condensability → B compartment (heterochromatin).
        We flip the eigenvector sign so that PC1 > 0 = low condensability = A.

    Returns
    -------
    pc1 : np.ndarray, shape (N,)
        Compartment score (positive = A, negative = B)
    eigenvalues : np.ndarray
        All eigenvalues
    """
    # Replace NaN with 0 for correlation
    oe_clean = np.nan_to_num(oe, nan=0.0)

    # Pearson correlation matrix (row-wise)
    # Centre each row
    row_means = oe_clean.mean(axis=1, keepdims=True)
    row_stds  = oe_clean.std(axis=1,  keepdims=True)
    row_stds[row_stds == 0] = 1.0
    oe_z = (oe_clean - row_means) / row_stds

    corr = (oe_z @ oe_z.T) / oe_clean.shape[1]
    np.fill_diagonal(corr, 1.0)

    # Eigendecomposition (eigh for symmetric matrices, returns ascending order)
    eigenvalues, eigenvectors = eigh(corr)
    # Largest eigenvalue is last
    pc1 = eigenvectors[:, -1]

    return pc1, eigenvalues, corr


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_contact_map(contact: np.ndarray, out_path: Path,
                     threshold: float, title: str = "Contact probability map"):
    n = contact.shape[0]
    c = contact.copy().astype(float)
    c[c < threshold] = np.nan

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(
        np.log10(c),
        cmap="RdYlBu_r",
        origin="upper",
        interpolation="nearest",
        aspect="auto"
    )
    plt.colorbar(im, ax=ax, label="log₁₀(contact probability)")
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Bead index")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_oe_map(oe: np.ndarray, out_path: Path,
                title: str = "O/E contact map"):
    fig, ax = plt.subplots(figsize=(7, 6))
    vmax = np.nanpercentile(np.abs(oe), 98)
    im = ax.imshow(
        oe,
        cmap="RdBu_r",
        origin="upper",
        vmin=-vmax, vmax=vmax,
        interpolation="nearest",
        aspect="auto"
    )
    plt.colorbar(im, ax=ax, label="O/E ratio")
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Bead index")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_compartment_track(pc1: np.ndarray, phi: np.ndarray,
                            bin_starts: np.ndarray,
                            chrom: str, out_path: Path):
    # Use a gridspec: top two panels share genomic x-axis; bottom is independent scatter
    fig = plt.figure(figsize=(12, 9))
    gs  = fig.add_gridspec(3, 1, hspace=0.45)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    ax2 = fig.add_subplot(gs[2])          # independent — scatter plot

    mb = bin_starts / 1e6   # Mb coordinates for genomic tracks

    # ── Panel 1: Compartment score (PC1) ──────────────────────────────
    ax0.fill_between(mb, pc1, where=(pc1 >= 0),
                     color="#e74c3c", alpha=0.8, label="A compartment (PC1>0)")
    ax0.fill_between(mb, pc1, where=(pc1 < 0),
                     color="#3498db", alpha=0.8, label="B compartment (PC1<0)")
    ax0.axhline(0, color="black", lw=0.8, ls="--")
    ax0.set_ylabel("Compartment\nPC1 score")
    ax0.legend(fontsize=8, loc="upper right")
    ax0.set_title(
        f"A/B Compartments — {chrom}:{bin_starts[0]/1e6:.2f}–{bin_starts[-1]/1e6:.2f} Mb",
        fontsize=11
    )
    plt.setp(ax0.get_xticklabels(), visible=False)

    # ── Panel 2: φ (condensability) genomic track ─────────────────────
    ax1.plot(mb, phi, color="#2c3e50", lw=0.8)
    ax1.axhline(phi.mean(), color="grey", lw=0.6, ls="--",
                label=f"mean φ = {phi.mean():.3f}")
    ax1.legend(fontsize=7, loc="upper right")
    ax1.set_ylabel("φ\n(condensability)")
    ax1.set_xlabel("Genomic position (Mb)")
    ax1.set_xlim(mb[0], mb[-1])

    # ── Panel 3: Scatter PC1 vs φ (independent x-axis) ────────────────
    r, pval = pearsonr(pc1, phi)
    ax2.scatter(phi, pc1, s=12, alpha=0.7, c="#8e44ad", edgecolors="none")
    ax2.axhline(0, color="grey", lw=0.5, ls="--")
    ax2.axvline(phi.mean(), color="grey", lw=0.5, ls="--")
    # Regression line
    m, b = np.polyfit(phi, pc1, 1)
    x_line = np.linspace(phi.min(), phi.max(), 100)
    ax2.plot(x_line, m * x_line + b, color="#c0392b", lw=1.5, ls="-",
             label=f"slope={m:.3f}")
    ax2.legend(fontsize=8)
    ax2.set_xlabel("φ (condensability)")
    ax2.set_ylabel("Compartment PC1")
    ax2.set_title(f"PC1 vs φ  (r = {r:.3f}, p = {pval:.2e})", fontsize=10)
    ax2.set_xlim(phi.min() - 0.02, phi.max() + 0.02)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_correlation_matrix(corr: np.ndarray, out_path: Path):
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1,
                   origin="upper", aspect="auto", interpolation="nearest")
    plt.colorbar(im, ax=ax, label="Pearson r")
    ax.set_title("O/E Pearson correlation matrix", fontsize=12)
    ax.set_xlabel("Bead index")
    ax.set_ylabel("Bead index")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# BedGraph export
# ---------------------------------------------------------------------------
def save_bedgraph(pc1: np.ndarray, bin_starts: np.ndarray,
                  binsize: int, chrom: str, out_path: Path):
    with open(out_path, "w") as f:
        f.write(f'track type=bedGraph name="AB_compartment_PC1" '
                f'description="A/B compartment score from polymer sim"\n')
        for i, (start, score) in enumerate(zip(bin_starts, pc1)):
            end = start + binsize
            f.write(f"{chrom}\t{start}\t{end}\t{score:.6f}\n")
    print(f"  Saved: {out_path}")


def save_condensability_bedgraph(phi: np.ndarray, raw: np.ndarray,
                                  bin_starts: np.ndarray,
                                  binsize: int, chrom: str, out_path: Path):
    with open(out_path, "w") as f:
        f.write(f'track type=bedGraph name="condensability_raw" '
                f'description="Raw condensability scores (binned)"\n')
        for i, (start, score) in enumerate(zip(bin_starts, raw)):
            end = start + binsize
            f.write(f"{chrom}\t{start}\t{end}\t{score:.6f}\n")
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n[Step 1] Loading contact matrix ...")
    contact = np.load(args.contact).astype(float)
    n = contact.shape[0]
    print(f"  Shape: {contact.shape}  (N={n} beads)")
    print(f"  Mean contact prob: {contact.mean():.4f}")

    print("\n[Step 2] Loading φ vector ...")
    phi_data   = np.load(args.phi, allow_pickle=True)
    phi        = phi_data["phi"]
    raw_scores = phi_data["raw_scores"]
    bin_starts = phi_data["bin_starts"]
    binsize    = int(phi_data["binsize"].flat[0])
    chrom      = str(phi_data["chrom"][0])
    r_start    = int(phi_data["region_start"].flat[0])
    r_end      = int(phi_data["region_end"].flat[0])

    print(f"  Region : {chrom}:{r_start:,}-{r_end:,}")
    print(f"  Binsize: {binsize/1e3:.0f} kb")

    # ---- Contact map ------------------------------------------------
    print("\n[Step 3] Plotting contact map ...")
    plot_contact_map(
        contact,
        out_dir / "contact_map.png",
        args.contact_threshold,
        title=f"Contact probability map — {chrom}:{r_start//1000}–{r_end//1000} kb"
    )

    # ---- O/E --------------------------------------------------------
    print("\n[Step 4] Computing O/E matrix ...")
    oe = observed_over_expected(contact)
    plot_oe_map(oe, out_dir / "contact_map_OE.png",
                title=f"O/E contact map — {chrom}:{r_start//1000}–{r_end//1000} kb")
    np.save(out_dir / "oe_matrix.npy", oe)
    print(f"  Saved: {out_dir/'oe_matrix.npy'}")

    # ---- Compartments -----------------------------------------------
    print("\n[Step 5] Computing A/B compartments (PCA on O/E correlation) ...")
    pc1, eigenvalues, corr = ab_compartments_pca(oe)

    # Sign convention: A compartment (low condensability) = positive PC1
    # Correlate with raw scores: if positive correlation, flip
    r_sign, _ = pearsonr(pc1, raw_scores)
    if r_sign > 0:
        pc1 = -pc1
        print(f"  PC1 sign flipped (correlation with condensability was {r_sign:.3f})")
    else:
        print(f"  PC1 sign kept (correlation with condensability: {r_sign:.3f})")

    # Variance explained
    var_explained = eigenvalues[-1] / eigenvalues.sum() * 100
    print(f"  Variance explained by PC1: {var_explained:.1f}%")

    # A/B fractions
    n_A = (pc1 > 0).sum()
    n_B = (pc1 < 0).sum()
    print(f"  A compartment beads: {n_A} ({100*n_A/n:.1f}%)")
    print(f"  B compartment beads: {n_B} ({100*n_B/n:.1f}%)")

    # Plots
    plot_correlation_matrix(corr, out_dir / "oe_correlation_matrix.png")
    plot_compartment_track(pc1, phi, bin_starts, chrom,
                           out_dir / "compartment_track.png")

    # ---- Save outputs -----------------------------------------------
    print("\n[Step 6] Saving outputs ...")
    np.save(out_dir / "compartment_pc1.npy", pc1)
    save_bedgraph(pc1, bin_starts, binsize, chrom,
                  out_dir / "compartment_scores.bedgraph")
    save_condensability_bedgraph(phi, raw_scores, bin_starts, binsize, chrom,
                                  out_dir / "condensability_binned.bedgraph")

    # Summary text
    summary_path = out_dir / "summary.txt"
    with open(summary_path, "w") as f:
        f.write("=== Chromatin Polymer Simulation Summary ===\n\n")
        f.write(f"Region       : {chrom}:{r_start:,}-{r_end:,}\n")
        f.write(f"Beads        : {n}  (binsize = {binsize/1e3:.0f} kb)\n")
        f.write(f"PC1 variance : {var_explained:.1f}%\n")
        f.write(f"A compartment: {n_A} beads ({100*n_A/n:.1f}%)\n")
        f.write(f"B compartment: {n_B} beads ({100*n_B/n:.1f}%)\n")
        f.write(f"PC1 vs φ corr: r = {pearsonr(pc1, phi)[0]:.4f}\n")
        f.write(f"PC1 vs raw   : r = {pearsonr(pc1, raw_scores)[0]:.4f}\n")
    print(f"  Saved: {summary_path}")

    print("\n✓  Analysis complete.  Outputs:")
    for p in sorted(out_dir.glob("*")):
        print(f"    {p}")
    print()


if __name__ == "__main__":
    main()
