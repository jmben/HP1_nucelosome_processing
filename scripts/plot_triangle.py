#!/usr/bin/env python3
"""
Standalone triangle contact map plotter.
Loads contact_matrix_avg.npy from each results directory and saves triangle maps.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration — adjust these to match your setup
# ---------------------------------------------------------------------------
BASE_DIR   = Path(".")   # directory containing your results_* folders
RESULTS    = [
    "results_supernatant_rep1_chr2_25kb_gpu",
    "results_supernatant_rep2_chr2_25kb_gpu",
]
CHROM      = "chr2"
START_BP   = 10_000_000   # genomic start of the region in bp
BINSIZE    = 25_000       # 25 kb bins

# ---------------------------------------------------------------------------

def plot_triangle_contact_map(contact: np.ndarray, out_path: Path,
                               title: str = "Contact probability map",
                               bin_starts: np.ndarray = None,
                               binsize: int = None):
    """
    Triangle Hi-C contact map.
    Apex at top center, diagonal runs along the bottom.
    Uses scipy rotate to achieve the standard published triangle orientation.
    """
    n = contact.shape[0]
    c = contact.copy().astype(float)

    # Symmetrise (in case matrix is upper-triangle only)
    c = np.triu(c) + np.triu(c, 1).T

    # Clip to 99th percentile (excluding diagonal)
    off_diag = c[~np.eye(n, dtype=bool)]
    vmax = np.percentile(off_diag[off_diag > 0], 99) if (off_diag > 0).any() else 1.0

    # Zero out lower triangle, keep upper only
    mat = c.copy()
    mat[np.tril_indices(n, -1)] = 0.0

    # Rotate 45 degrees clockwise — diagonal becomes the base of the triangle
    # cval=0 fills background with 0 (will appear white on Reds colormap)
    rotated = rotate(mat, angle=45, reshape=True, cval=0.0, order=1)

    # After rotation the triangle sits in the top half of the rotated image
    h, w = rotated.shape
    triangle = rotated[:h // 2, :]

    # Build a simple triangle mask based on geometry
    rows, cols = np.indices(triangle.shape)
    tri_h, tri_w = triangle.shape
    # Triangle: apex at top center, base at bottom
    # Left edge:  col >= center - row * (center/tri_h)
    # Right edge: col <= center + row * (center/tri_h)
    center = tri_w / 2
    slope  = center / tri_h
    left_bound  = center - rows * slope
    right_bound = center + rows * slope
    outside_triangle = (cols < left_bound) | (cols > right_bound)

    img_masked = np.ma.array(triangle, mask=outside_triangle)

    # --- Plot ---
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

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label="Contact probability",
                        shrink=0.6, pad=0.02, orientation="vertical")
    cbar.ax.tick_params(labelsize=9)

    # X-axis ticks
    # After 45° rotation of an n x n matrix, the base of the triangle
    # spans the full width w of the rotated image.
    # Bin i on the diagonal maps linearly: pixel_x = i * (w / n)
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
# Main
# ---------------------------------------------------------------------------
for res_dir_name in RESULTS:
    res_dir  = BASE_DIR / res_dir_name
    npy_path = res_dir / "contact_matrix_avg.npy"

    if not npy_path.exists():
        print(f"  Skipping {res_dir_name}: contact_matrix_avg.npy not found")
        continue

    print(f"Processing {res_dir_name}...")
    matrix = np.load(npy_path)
    n = matrix.shape[0]

    bin_starts = np.arange(n) * BINSIZE + START_BP

    label = res_dir_name.replace("results_", "").replace("_gpu", "").replace("_", " ")
    title = f"Contact probability map — {CHROM} ({label})"

    out_path = res_dir / "contact_map_triangle.png"
    plot_triangle_contact_map(
        matrix, out_path,
        title=title,
        bin_starts=bin_starts,
        binsize=BINSIZE,
    )

print("Done.")
