#!/usr/bin/env python3
"""
Compare contrast enhancement strategies for triangle contact maps.
Generates 4 versions: original + 3 contrast variants.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from pathlib import Path
import cooler

# ---------------------------------------------------------------------------
MCOOL_PATH = Path("/home/gdbishoplab/jxb377/condense-seq/GSE246947_hg38_MCF10A_WT_merged.mcool")
RESOLUTION = 20_000
CHROM      = "chr2"
START_BP   = 10_000_000
END_BP     = 14_000_000
OUT_DIR    = Path("/home/gdbishoplab/jxb377/condense-seq")

# ---------------------------------------------------------------------------

def _rotate_and_mask(mat):
    """Shared rotation + masking logic."""
    rotated = rotate(mat, angle=45, reshape=True, cval=0.0, order=1)
    h, w    = rotated.shape
    triangle = rotated[:h // 2, :]
    rows, cols = np.indices(triangle.shape)
    tri_h, tri_w = triangle.shape
    center      = tri_w / 2
    slope       = center / tri_h
    left_bound  = center - rows * slope
    right_bound = center + rows * slope
    outside     = (cols < left_bound) | (cols > right_bound)
    return np.ma.array(triangle, mask=outside), w


def _add_ticks(ax, w, n, bin_starts, binsize):
    mb       = bin_starts / 1e6
    mb_end   = mb[-1] + binsize / 1e6
    tick_mb  = np.linspace(mb[0], mb_end, 7)
    tick_idx = [(np.abs(mb - t)).argmin() for t in tick_mb]
    tick_idx = sorted(set(np.clip(tick_idx, 0, n - 1)))
    ax.set_xticks([int(i * w / n) for i in tick_idx])
    ax.set_xticklabels([f"{mb[i]:.1f}" for i in tick_idx],
                       fontsize=9, rotation=45, ha="right")
    ax.set_xlabel("Genomic position (Mb)", fontsize=11)


def plot_triangle(img_masked, w, vmin, vmax, cmap, cbar_label,
                  title, out_path, bin_starts, binsize, n):
    fig, ax = plt.subplots(figsize=(12, 5))
    im = ax.imshow(img_masked, cmap=cmap, vmin=vmin, vmax=vmax,
                   aspect="auto", interpolation="nearest", origin="upper")
    cbar = plt.colorbar(im, ax=ax, label=cbar_label,
                        shrink=0.6, pad=0.02, orientation="vertical")
    cbar.ax.tick_params(labelsize=9)
    _add_ticks(ax, w, n, bin_starts, binsize)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=11, pad=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
clr    = cooler.Cooler(f"{MCOOL_PATH}::resolutions/{RESOLUTION}")
region = f"{CHROM}:{START_BP}-{END_BP}"
print(f"Fetching {region} at {RESOLUTION//1000}kb...")

raw    = clr.matrix(balance=True).fetch(region)
bins   = clr.bins().fetch(region)
bin_starts = bins["start"].values

n = raw.shape[0]

# Symmetrise + clean
c = np.triu(raw) + np.triu(raw, 1).T
c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)

# Zero lower triangle for rotation
mat = c.copy()
mat[np.tril_indices(n, -1)] = 0.0

off_diag = c[~np.eye(n, dtype=bool)]
off_diag_pos = off_diag[off_diag > 0]

base_title = f"MCF10A WT — {CHROM}:{START_BP//1_000_000}-{END_BP//1_000_000} Mb (20 kb)"

# ---------------------------------------------------------------------------
# 1. Original — 99th percentile vmax
# ---------------------------------------------------------------------------
vmax_99 = np.percentile(off_diag_pos, 99)
img1, w  = _rotate_and_mask(mat)
plot_triangle(img1, w, 0, vmax_99, "Reds", "Contact frequency",
              f"{base_title} | 99th pct vmax",
              OUT_DIR / "triangle_99pct.png",
              bin_starts, RESOLUTION, n)

# ---------------------------------------------------------------------------
# 2. Lower percentile vmax — 90th percentile
# ---------------------------------------------------------------------------
vmax_90 = np.percentile(off_diag_pos, 90)
img2, _  = _rotate_and_mask(mat)
plot_triangle(img2, w, 0, vmax_90, "Reds", "Contact frequency",
              f"{base_title} | 90th pct vmax",
              OUT_DIR / "triangle_90pct.png",
              bin_starts, RESOLUTION, n)

# ---------------------------------------------------------------------------
# 3. Fixed vmax — median of top 5% (robust high-end clip)
# ---------------------------------------------------------------------------
vmax_fixed = np.percentile(off_diag_pos, 95)
img3, _    = _rotate_and_mask(mat)
plot_triangle(img3, w, 0, vmax_fixed, "Reds", "Contact frequency",
              f"{base_title} | 95th pct vmax",
              OUT_DIR / "triangle_95pct.png",
              bin_starts, RESOLUTION, n)

# ---------------------------------------------------------------------------
# 4. Log scale
# ---------------------------------------------------------------------------
log_mat = mat.copy()
log_mat[log_mat > 0] = np.log1p(log_mat[log_mat > 0])
img4, _  = _rotate_and_mask(log_mat)

log_off  = log_mat[~np.eye(n, dtype=bool)]
vmax_log = np.percentile(log_off[log_off > 0], 99)
plot_triangle(img4, w, 0, vmax_log, "Reds", "log(1 + contact frequency)",
              f"{base_title} | log scale",
              OUT_DIR / "triangle_log.png",
              bin_starts, RESOLUTION, n)

print("Done. Generated: triangle_99pct.png, triangle_90pct.png, triangle_95pct.png, triangle_log.png")
