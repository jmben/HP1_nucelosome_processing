#!/usr/bin/env python3
"""
02_build_energy_matrix.py  —  Compute ε_ij interaction matrix
=============================================================================
Implements the pairwise interaction energy from the paper (Supp. Note 2):

    ε_ij = -β ( φ_i · φ_j  -  φ̄² )

where
    β    = 0.6   (paper default)
    φ̄    = 0.6   (centering constant, same as used in 01_bin_bigwig.py)
    φ_i  = normalised condensability score of bead i

The full contact potential depth is:
    prefactor_ij = ε₀ / |j-i|^α  +  ε_ij
with ε₀ = -0.5, α = -0.5 (i.e. |j-i|^0.5 in denominator, but note
the paper writes α = -0.5 meaning the exponent is negative so the term
grows with separation — see note in code).

Saves an (N×N) matrix to data/energy_matrix.npz for the simulator.

Usage:
    python 02_build_energy_matrix.py \
        --phi data/phi_vector.npz \
        --out  data/energy_matrix.npz \
        [--beta 0.6] [--eps0 -0.5] [--alpha -0.5]
=============================================================================
"""

import argparse
import numpy as np
from pathlib import Path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Build ε_ij interaction matrix.")
    p.add_argument("--phi",          default="data/phi_vector.npz", help="φ vector from step 1")
    p.add_argument("--out",          default="data/energy_matrix.npz")
    p.add_argument("--beta",         type=float, default=0.6,   help="β scaling factor (paper: 0.6)")
    p.add_argument("--eps0",         type=float, default=-0.5,  help="ε₀ for separation term (paper: -0.5)")
    p.add_argument("--alpha",        type=float, default=-0.5,  help="α exponent for separation (paper: -0.5)")
    p.add_argument("--no_sep_term",  action="store_true",
                   help="Drop the genomic-separation term (ε₀/|j-i|^α) and use "
                        "only the sequence-dependent ε_ij term.  Recommended when "
                        "the separation term dominates and drowns sequence signal.")
    p.add_argument("--beta_scale",   type=float, default=1.0,
                   help="Multiply β by this factor to amplify sequence-dependent "
                        "signal relative to other forces.  E.g. --beta_scale 30 "
                        "gives effective β=18 when --beta=0.6.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Matrix construction
# ---------------------------------------------------------------------------
def build_epsilon_ij(phi: np.ndarray, beta: float, phi_bar: float) -> np.ndarray:
    """
    Sequence-dependent interaction term:
        ε_ij = -β ( φ_i φ_j  -  φ̄² )

    Implemented as an outer product for efficiency.
    Shape: (N, N)
    """
    outer = np.outer(phi, phi)          # φ_i * φ_j  for all i,j
    eps_ij = -beta * (outer - phi_bar ** 2)
    return eps_ij


def build_separation_term(n: int, eps0: float, alpha: float) -> np.ndarray:
    """
    Genomic separation-dependent term:
        T_ij = ε₀ / |j - i|^α

    Note: with ε₀ = -0.5 and α = -0.5, the term is:
        T_ij = -0.5 / |j-i|^(-0.5) = -0.5 * |j-i|^0.5

    i.e. more negative (stronger attraction) for beads further apart —
    this counteracts the polymer's tendency to form local contacts and
    reproduces the power-law Hi-C contact frequency decay after
    heavy coarse-graining (paper Supp. Note 2).
    The diagonal (i == j) is set to 0.
    """
    sep = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = abs(j - i)
            val = eps0 / (d ** alpha)   # d^alpha where alpha=-0.5 → d^(-(-0.5)) = d^0.5
            sep[i, j] = val
            sep[j, i] = val
    return sep


def build_separation_term_fast(n: int, eps0: float, alpha: float) -> np.ndarray:
    """Vectorised version of build_separation_term (much faster for large N)."""
    idx  = np.arange(n)
    dist = np.abs(idx[:, None] - idx[None, :]).astype(float)
    np.fill_diagonal(dist, 1.0)          # avoid div-by-zero; diagonal zeroed below
    with np.errstate(divide="ignore", invalid="ignore"):
        sep = eps0 / (dist ** alpha)
    np.fill_diagonal(sep, 0.0)
    return sep


def build_full_prefactor(phi: np.ndarray, beta: float, phi_bar: float,
                         eps0: float, alpha: float,
                         no_sep_term: bool = False) -> tuple:
    """
    Full prefactor for the contact potential:
        prefactor_ij = T_ij  +  ε_ij / 2

    where T_ij = ε₀/|j-i|^α  (separation-dependent correction)
    and   ε_ij = -β(φ_i φ_j - φ̄²)  (sequence-dependent term)

    If no_sep_term=True, T_ij is set to zero and the prefactor is
    purely sequence-dependent.  Recommended when the separation term
    dominates the sequence signal (ratio > 5×).
    """
    n      = len(phi)
    eps_ij = build_epsilon_ij(phi, beta, phi_bar)

    if no_sep_term:
        sep_term = np.zeros((n, n))
        print("  Separation term: DISABLED (using ε_ij only)")
    else:
        sep_term = build_separation_term_fast(n, eps0, alpha)

    prefactor = sep_term + eps_ij / 2.0
    np.fill_diagonal(prefactor, 0.0)
    return prefactor, eps_ij, sep_term


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    print("\n[Step 1] Loading φ vector ...")
    data = np.load(args.phi, allow_pickle=True)
    phi      = data["phi"]
    phi_bar  = float(data["phi_bar"].flat[0])
    n_beads  = len(phi)
    binsize  = int(data["binsize"].flat[0])
    chrom    = str(data["chrom"][0])
    r_start  = int(data["region_start"].flat[0])
    r_end    = int(data["region_end"].flat[0])

    print(f"  Beads   : {n_beads}")
    print(f"  Region  : {chrom}:{r_start:,}-{r_end:,}")
    print(f"  φ̄       : {phi_bar}")
    effective_beta = args.beta * args.beta_scale
    print(f"  β       : {args.beta}  ×  scale {args.beta_scale}  =  {effective_beta:.2f} (effective)")
    print(f"  ε₀      : {args.eps0}   α: {args.alpha}")
    print(f"  sep term: {'DISABLED' if args.no_sep_term else 'enabled'}")

    print("\n[Step 2] Building interaction matrices ...")
    prefactor, eps_ij, sep_term = build_full_prefactor(
        phi, effective_beta, phi_bar, args.eps0, args.alpha,
        no_sep_term=args.no_sep_term
    )

    print(f"  ε_ij    : min={eps_ij.min():.4f}  max={eps_ij.max():.4f}")
    nz = sep_term[sep_term != 0]
    if len(nz):
        print(f"  sep term: min={nz.min():.4f}  max={nz.max():.4f}")
    else:
        print(f"  sep term: DISABLED (all zeros)")
    print(f"  prefactor: min={prefactor.min():.4f}  max={prefactor.max():.4f}")
    print(f"  Pairs attractive (pre < 0): {(prefactor < 0).sum()}")
    print(f"  Pairs repulsive  (pre > 0): {(prefactor > 0).sum()}")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"\n[Step 3] Saving to {out_path} ...")
    np.savez_compressed(
        out_path,
        prefactor=prefactor,
        eps_ij=eps_ij,
        sep_term=sep_term,
        phi=phi,
        phi_bar=np.array([phi_bar]),
        beta=np.array([args.beta]),
        eps0=np.array([args.eps0]),
        alpha=np.array([args.alpha]),
        n_beads=np.array([n_beads]),
        chrom=np.array([chrom]),
        region_start=np.array([r_start]),
        region_end=np.array([r_end]),
        binsize=np.array([binsize]),
    )

    print(f"\n✓  Energy matrix ({n_beads}×{n_beads}) saved  →  {out_path}\n")


if __name__ == "__main__":
    main()
