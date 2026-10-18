#!/usr/bin/env python3
"""
Stacked triangle contact map comparison.
Top:    MCF10A WT Hi-C (mcool, 20 kb, log scale, 95th pct)
Middle: Simulation — 20U supernatant rep1 (npy, 25 kb)
Bottom: Simulation — Genomic rep1 sensitivity (npy, 25 kb)
Shared x-axis for direct TAD-boundary comparison.

Outputs
-------
comparison_triangle_clean.png
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import rotate
from pathlib import Path
import cooler

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MCOOL_PATH      = Path("/home/gdbishoplab/jxb377/condense-seq/GSE246947_hg38_MCF10A_WT_merged.mcool")
SUP_REP1_NPY    = Path("/home/gdbishoplab/jxb377/condense-seq/results_supernatant_rep1_chr2_25kb_gpu/contact_matrix_avg.npy")
SENS_NPY        = Path("/home/gdbishoplab/jxb377/condense-seq/results_rep1_chr2_25kb_gpu/contact_matrix_avg.npy")
OUT_PATH        = Path("/home/gdbishoplab/jxb377/condense-seq/comparison_triangle_clean.png")

CHROM       = "chr2"
START_BP    = 10_000_000
END_BP      = 14_000_000
MCOOL_RES   = 20_000
SIM_BINSIZE = 25_000

PCT_HIC     = 95   # percentile vmax for Hi-C panel
PCT_SIM     = 90   # percentile vmax for simulation panels

# Fraction of triangle height to show from base up (0.35 = bottom 35%)
TRIM_FRAC   = 0.35

# ---------------------------------------------------------------------------


def prepare_triangle(mat, pct, log_scale=False):
    """
    Symmetrise, clean, optionally log-transform, rotate 45 deg,
    trim to TRIM_FRAC of height, and mask outside triangle.

    Returns
    -------
    img_masked : np.ma.MaskedArray
    vmax       : float
    w          : int   (rotated image width, used for tick placement)
    n          : int   (number of bins)
    trim_row   : int
    full_h     : int
    """
    n = mat.shape[0]
    c = mat.copy().astype(float)
    c = np.triu(c) + np.triu(c, 1).T
    c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(c, 0.0)   # zero diagonal so it doesn't dominate vmax

    # Log transform BEFORE computing vmax so scale is correct
    if log_scale:
        c = np.log1p(c)

    # vmax from off-diagonal positive values
    off_diag     = c[~np.eye(n, dtype=bool)]
    off_diag_pos = off_diag[off_diag > 0]
    vmax = np.percentile(off_diag_pos, pct) if len(off_diag_pos) > 0 else 1.0

    # Clip to vmax before rotation to avoid rotation artefacts
    c = np.clip(c, 0, vmax)

    # Upper triangle only for rotation
    upper = np.triu(c)

    # Rotate 45 degrees — diagonal becomes base of triangle
    rotated  = rotate(upper, angle=45, reshape=True, cval=0.0, order=1)
    h, w     = rotated.shape
    full_h   = h // 2
    triangle = rotated[:full_h, :]

    # Trim: keep only bottom TRIM_FRAC rows (short-range contacts)
    trim_row = int(full_h * (1 - TRIM_FRAC))
    trimmed  = triangle[trim_row:, :]

    # Build clean triangle mask using full geometry
    rows, cols = np.indices(trimmed.shape)
    center = w / 2
    slope  = center / full_h   # slope based on FULL height for correct edges

    # Absolute row position in full triangle
    abs_rows = rows + trim_row
    left  = center - abs_rows * slope
    right = center + abs_rows * slope

    # Small buffer to avoid clipping edge pixels
    buf     = 1.5
    outside = (cols < left - buf) | (cols > right + buf)

    img_masked = np.ma.array(trimmed, mask=outside)
    return img_masked, vmax, w, n, trim_row, full_h


def add_ticks(ax, w, n, bin_starts, binsize, fontsize=9):
    """Add genomic Mb tick labels to a triangle axes."""
    mb      = bin_starts / 1e6
    mb_end  = mb[-1] + binsize / 1e6
    tick_mb = np.linspace(mb[0], mb_end, 7)
    tidx    = [(np.abs(mb - t)).argmin() for t in tick_mb]
    tidx    = sorted(set(np.clip(tidx, 0, n - 1)))
    ax.set_xticks([int(i * w / n) for i in tidx])
    ax.set_xticklabels([f"{mb[i]:.1f}" for i in tidx],
                       fontsize=fontsize, rotation=45, ha="right")


def draw_panel(ax, img_masked, vmax, w, n, bin_starts, binsize,
               label, cbar_label, show_xticks=False):
    """Render one triangle panel."""
    im = ax.imshow(img_masked, cmap="Reds", vmin=0, vmax=vmax,
                   aspect="auto", interpolation="nearest", origin="upper")
    cbar = plt.colorbar(im, ax=ax, label=cbar_label,
                        shrink=0.7, pad=0.01, orientation="vertical")
    cbar.ax.tick_params(labelsize=8)

    if show_xticks:
        add_ticks(ax, w, n, bin_starts, binsize)
        ax.set_xlabel("Genomic position (Mb)", fontsize=11)
    else:
        ax.set_xticks([])

    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.text(0.01, 0.97, label, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top", ha="left")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading Hi-C from mcool...")
clr          = cooler.Cooler(f"{MCOOL_PATH}::resolutions/{MCOOL_RES}")
region       = f"{CHROM}:{START_BP}-{END_BP}"
hic_mat      = clr.matrix(balance=True).fetch(region)
hic_bins     = clr.bins().fetch(region)["start"].values

print("Loading 20U supernatant rep1 simulation...")
sup_rep1_mat  = np.load(SUP_REP1_NPY)
sup_rep1_bins = np.arange(sup_rep1_mat.shape[0]) * SIM_BINSIZE + START_BP

print("Loading genomic sensitivity simulation...")
sens_mat      = np.load(SENS_NPY)
sens_bins     = np.arange(sens_mat.shape[0]) * SIM_BINSIZE + START_BP

# ---------------------------------------------------------------------------
# Prepare triangles
# ---------------------------------------------------------------------------
hic_img,       hic_vmax,       hic_w,       hic_n,       *_ = prepare_triangle(
    hic_mat,      PCT_HIC, log_scale=True)
sup_rep1_img,  sup_rep1_vmax,  sup_rep1_w,  sup_rep1_n,  *_ = prepare_triangle(
    sup_rep1_mat, PCT_SIM, log_scale=False)
sens_img,      sens_vmax,      sens_w,      sens_n,      *_ = prepare_triangle(
    sens_mat,     PCT_SIM, log_scale=False)

# ---------------------------------------------------------------------------
# Plot — equal-height panels
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(14, 8))
gs  = GridSpec(3, 1, figure=fig, hspace=0.1,
               height_ratios=[1, 1, 1])

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])

draw_panel(ax1, hic_img,      hic_vmax,      hic_w,      hic_n,
           hic_bins,     MCOOL_RES,
           label=f"MCF10A WT Hi-C — log scale ({PCT_HIC}th pct)",
           cbar_label="log(1 + freq.)",
           show_xticks=False)

draw_panel(ax2, sup_rep1_img, sup_rep1_vmax, sup_rep1_w, sup_rep1_n,
           sup_rep1_bins, SIM_BINSIZE,
           label=f"Simulation — 20U supernatant rep1 ({PCT_SIM}th pct)",
           cbar_label="Contact prob.",
           show_xticks=False)

draw_panel(ax3, sens_img,     sens_vmax,     sens_w,     sens_n,
           sens_bins,    SIM_BINSIZE,
           label=f"Simulation — Genomic sensitivity ({PCT_SIM}th pct)",
           cbar_label="Contact prob.",
           show_xticks=True)

fig.suptitle(f"{CHROM}:{START_BP//1_000_000}–{END_BP//1_000_000} Mb",
             fontsize=13, fontweight="bold", y=0.99)

fig.savefig(OUT_PATH, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT_PATH}")
