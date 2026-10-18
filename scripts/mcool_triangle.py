#!/usr/bin/env python3
"""
Plot triangle contact map from .mcool file using cooler.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from pathlib import Path
import cooler

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MCOOL_PATH = Path("/home/gdbishoplab/jxb377/condense-seq/GSE246947_hg38_MCF10A_WT_merged.mcool")
RESOLUTION = 20_000       # 20 kb
CHROM      = "chr2"
START_BP   = 10_000_000
END_BP     = 14_000_000
OUT_PATH   = Path("/home/gdbishoplab/jxb377/condense-seq/contact_map_triangle_mcool.png")

# ---------------------------------------------------------------------------

def plot_triangle_contact_map(contact: np.ndarray, out_path: Path,
                               title: str = "Contact probability map",
                               bin_starts: np.ndarray = None,
                               binsize: int = None):
    n = contact.shape[0]
    c = contact.copy().astype(float)

    # Symmetrise
    c = np.triu(c) + np.triu(c, 1).T

    # Replace NaN/inf
    c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)

    # Clip to 99th percentile (excluding diagonal)
    off_diag = c[~np.eye(n, dtype=bool)]
    vmax = np.percentile(off_diag[off_diag > 0], 99) if (off_diag > 0).any() else 1.0

    # Zero out lower triangle
    mat = c.copy()
    mat[np.tril_indices(n, -1)] = 0.0

    # Rotate 45 degrees — diagonal becomes base of triangle
    rotated = rotate(mat, angle=45, reshape=True, cval=0.0, order=1)

    # Crop to top half (the triangle)
    h, w = rotated.shape
    triangle = rotated[:h // 2, :]

    # Mask outside-triangle regions
    rows, cols = np.indices(triangle.shape)
    tri_h, tri_w = triangle.shape
    center      = tri_w / 2
    slope       = center / tri_h
    left_bound  = center - rows * slope
    right_bound = center + rows * slope
    outside     = (cols < left_bound) | (cols > right_bound)
    img_masked  = np.ma.array(triangle, mask=outside)

    fig, ax = plt.subplots(figsize=(12, 5))
    im = ax.imshow(
        img_masked,
        cmap="Reds",
        vmin=0,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest",
        origin="upper",
    )
    cbar = plt.colorbar(im, ax=ax, label="Contact frequency",
                        shrink=0.6, pad=0.02, orientation="vertical")
    cbar.ax.tick_params(labelsize=9)

    if bin_starts is not None and binsize is not None:
        mb       = bin_starts / 1e6
        mb_start = mb[0]
        mb_end   = mb[-1] + binsize / 1e6
        tick_mb  = np.linspace(mb_start, mb_end, 7)
        tick_idx = [(np.abs(mb - t)).argmin() for t in tick_mb]
        tick_idx = sorted(set(np.clip(tick_idx, 0, n - 1)))
        pixel_ticks = [int(i * w / n) for i in tick_idx]
        ax.set_xticks(pixel_ticks)
        ax.set_xticklabels([f"{mb[i]:.1f}" for i in tick_idx],
                           fontsize=9, rotation=45, ha="right")
        ax.set_xlabel("Genomic position (Mb)", fontsize=11)
    else:
        ax.set_xlabel("Bin index", fontsize=10)

    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title(title, fontsize=11, pad=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Load from .mcool
# ---------------------------------------------------------------------------
# Open at 20kb resolution
clr = cooler.Cooler(f"{MCOOL_PATH}::resolutions/{RESOLUTION}")

# Fetch the matrix for the region
region = f"{CHROM}:{START_BP}-{END_BP}"
print(f"Fetching {region} at {RESOLUTION//1000}kb resolution...")

matrix = clr.matrix(balance=True).fetch(region)

# Get bin starts for tick labels
bins = clr.bins().fetch(region)
bin_starts = bins["start"].values

title = f"MCF10A WT Hi-C — {CHROM}:{START_BP//1_000_000}-{END_BP//1_000_000} Mb (20 kb)"

plot_triangle_contact_map(
    matrix, OUT_PATH,
    title=title,
    bin_starts=bin_starts,
    binsize=RESOLUTION,
)

print("Done.")
