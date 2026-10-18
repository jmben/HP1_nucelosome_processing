#!/usr/bin/env python3
"""
TAD boundary calling and comparison between simulation and Hi-C.
Clean publication figure version.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from scipy.ndimage import rotate, gaussian_filter1d
from scipy.signal import find_peaks
from pathlib import Path
import cooler

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MCOOL_PATH          = Path("/home/gdbishoplab/jxb377/condense-seq/GSE246947_hg38_MCF10A_WT_merged.mcool")
SUPERNATANT_REP1_NPY = Path("/home/gdbishoplab/jxb377/condense-seq/results_supernatant_rep1_chr2_25kb_gpu/contact_matrix_avg.npy")
SENS_NPY            = Path("/home/gdbishoplab/jxb377/condense-seq/results_rep1_chr2_25kb_gpu/contact_matrix_avg.npy")
OUT_COMPARE         = Path("/home/gdbishoplab/jxb377/condense-seq/tad_boundary_comparison_clean.png")
OUT_STATS           = Path("/home/gdbishoplab/jxb377/condense-seq/tad_boundary_stats.txt")

CHROM       = "chr2"
START_BP    = 10_000_000
END_BP      = 14_000_000
MCOOL_RES   = 20_000
SIM_BINSIZE = 25_000
TRIM_FRAC   = 0.35

WINDOW_HIC      = 10
WINDOW_SIM      = 8
SMOOTH_SIGMA    = 1.5
MIN_PROMINENCE  = 0.05
MIN_DISTANCE    = 3
BOUNDARY_TOL_KB = 200

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def symmetrise(mat):
    c = mat.copy().astype(float)
    c = np.triu(c) + np.triu(c, 1).T
    c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(c, 0.0)
    return c


def insulation_score(mat, window):
    n      = mat.shape[0]
    scores = np.full(n, np.nan)
    for i in range(window, n - window):
        block = mat[i - window:i, i + 1:i + window + 1]
        vals  = block[block > 0]
        if len(vals) > 0:
            scores[i] = vals.mean()
    valid = np.isfinite(scores) & (scores > 0)
    scores[valid] = np.log2(scores[valid])
    mean_score    = np.nanmean(scores[valid])
    scores[valid] -= mean_score
    scores = gaussian_filter1d(scores, sigma=SMOOTH_SIGMA, mode='nearest')
    return scores


def call_boundaries(scores, binsize,
                    min_prominence=MIN_PROMINENCE,
                    min_distance=MIN_DISTANCE):
    valid = np.isfinite(scores)
    s     = scores.copy()
    s[~valid] = np.nanmax(s)
    peaks, _ = find_peaks(-s, prominence=min_prominence,
                          distance=min_distance)
    return peaks, peaks * binsize


def boundary_overlap(bnd1_bp, bnd2_bp, tol_bp):
    matched1, matched2, used2 = [], [], set()
    for b1 in bnd1_bp:
        dists = np.abs(bnd2_bp - b1)
        best  = np.argmin(dists)
        if dists[best] <= tol_bp and best not in used2:
            matched1.append(b1)
            matched2.append(bnd2_bp[best])
            used2.add(best)
    n_m  = len(matched1)
    rec  = n_m / len(bnd1_bp) if len(bnd1_bp) > 0 else 0
    prec = n_m / len(bnd2_bp) if len(bnd2_bp) > 0 else 0
    f1   = (2 * prec * rec / (prec + rec)
            if (prec + rec) > 0 else 0)
    return {"n_bnd1": len(bnd1_bp), "n_bnd2": len(bnd2_bp),
            "n_matched": n_m, "recall": rec,
            "precision": prec, "f1": f1}


def prepare_triangle(mat, pct=95, log_scale=False, trim_frac=TRIM_FRAC):
    n = mat.shape[0]
    c = mat.copy().astype(float)
    c = np.triu(c) + np.triu(c, 1).T
    c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(c, 0.0)
    if log_scale:
        c = np.log1p(c)
    off_diag     = c[~np.eye(n, dtype=bool)]
    off_diag_pos = off_diag[off_diag > 0]
    vmax = np.percentile(off_diag_pos, pct) if len(off_diag_pos) > 0 else 1.0
    c    = np.clip(c, 0, vmax)
    upper    = np.triu(c)
    rotated  = rotate(upper, angle=45, reshape=True, cval=0.0, order=1)
    h, w     = rotated.shape
    full_h   = h // 2
    triangle = rotated[:full_h, :]
    trim_row = int(full_h * (1 - trim_frac))
    trimmed  = triangle[trim_row:, :]
    rows, cols = np.indices(trimmed.shape)
    center     = w / 2
    slope      = center / full_h
    abs_rows   = rows + trim_row
    outside    = (cols < center - abs_rows * slope - 1.5) | \
                 (cols > center + abs_rows * slope + 1.5)
    return np.ma.array(trimmed, mask=outside), vmax, w, n


def bin_to_pixel(bnd_mb, mb_start, mb_range, w):
    """Convert genomic Mb position to pixel x coordinate."""
    return [(b - mb_start) / mb_range * w for b in bnd_mb]


def set_mb_ticks(ax, mb_start, mb_end, n_ticks=7, fontsize=9):
    ticks = np.linspace(mb_start, mb_end, n_ticks)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.1f}" for t in ticks],
                       fontsize=fontsize, rotation=45, ha="right")
    ax.set_xlim(mb_start, mb_end)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading Hi-C...")
clr      = cooler.Cooler(f"{MCOOL_PATH}::resolutions/{MCOOL_RES}")
region   = f"{CHROM}:{START_BP}-{END_BP}"
hic_raw  = clr.matrix(balance=True).fetch(region)
hic_bins = clr.bins().fetch(region)["start"].values
hic_mat  = symmetrise(hic_raw)
n_hic    = hic_mat.shape[0]
hic_mb   = hic_bins / 1e6

print("Loading simulations...")
rep1_mat = symmetrise(np.load(SUPERNATANT_REP1_NPY))
sens_mat = symmetrise(np.load(SENS_NPY))
n_sim    = rep1_mat.shape[0]
sim_bins = np.arange(n_sim) * SIM_BINSIZE + START_BP
sim_mb   = sim_bins / 1e6

# ---------------------------------------------------------------------------
# Insulation scores + boundary calls
# ---------------------------------------------------------------------------
print("Computing insulation scores and boundaries...")
hic_ins  = insulation_score(hic_mat,  WINDOW_HIC)
rep1_ins = insulation_score(rep1_mat, WINDOW_SIM)
sens_ins = insulation_score(sens_mat, WINDOW_SIM)

hic_bnd_idx,  hic_bnd_bp  = call_boundaries(hic_ins,  MCOOL_RES)
rep1_bnd_idx, rep1_bnd_bp = call_boundaries(rep1_ins, SIM_BINSIZE)
sens_bnd_idx, sens_bnd_bp = call_boundaries(sens_ins, SIM_BINSIZE)

hic_bnd_mb  = hic_bnd_bp  / 1e6 + START_BP / 1e6
rep1_bnd_mb = rep1_bnd_bp / 1e6 + START_BP / 1e6
sens_bnd_mb = sens_bnd_bp / 1e6 + START_BP / 1e6

tol_bp       = BOUNDARY_TOL_KB * 1000
rep1_overlap = boundary_overlap(hic_bnd_bp, rep1_bnd_bp, tol_bp)
sens_overlap = boundary_overlap(hic_bnd_bp, sens_bnd_bp, tol_bp)

# ---------------------------------------------------------------------------
# Prepare triangle images
# ---------------------------------------------------------------------------
hic_img,  hic_vmax,  hic_w,  hic_n  = prepare_triangle(
    hic_mat,  pct=95, log_scale=True)
rep1_img, rep1_vmax, rep1_w, rep1_n = prepare_triangle(
    rep1_mat, pct=90)
sens_img, sens_vmax, sens_w, sens_n = prepare_triangle(
    sens_mat, pct=90)

# Genomic extents
hic_mb_start  = hic_mb[0]
hic_mb_end    = hic_mb[-1] + MCOOL_RES / 1e6
hic_mb_range  = hic_mb_end - hic_mb_start

sim_mb_start  = sim_mb[0]
sim_mb_end    = sim_mb[-1] + SIM_BINSIZE / 1e6
sim_mb_range  = sim_mb_end - sim_mb_start

# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
print("Building figure...")
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         10,
    "axes.linewidth":    1.1,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
})

colors = {"hic": "#222222", "rep1": "#e41a1c", "sens": "#377eb8"}

# Layout: 3 groups of (map + insulation), plus stats row at bottom
fig = plt.figure(figsize=(14, 13))
gs  = GridSpec(
    7, 1, figure=fig,
    height_ratios=[2.2, 0.8, 2.2, 0.8, 2.2, 0.8, 0.5],
    hspace=0.10,
    left=0.05, right=0.87,
    top=0.95,  bottom=0.08,
)

# Separate colorbar axes for Hi-C vs simulation
cbar_hic_ax = fig.add_axes([0.895, 0.72, 0.016, 0.21])
cbar_sim_ax = fig.add_axes([0.895, 0.18, 0.016, 0.42])

panels = [
    # (map_row, ins_row, img, vmax, w, n,
    #  mb_start, mb_end, mb_range, bins, binsize,
    #  ins, bnd_mb, color, label, show_x, cbar_key)
    (0, 1,
     hic_img,  hic_vmax,  hic_w,  hic_n,
     hic_mb_start, hic_mb_end, hic_mb_range,
     hic_bins,  MCOOL_RES,
     hic_ins,   hic_bnd_mb,
     colors["hic"], "MCF10A WT Hi-C (20 kb, log scale)",
     False, "hic"),

    (2, 3,
     rep1_img, rep1_vmax, rep1_w, rep1_n,
     sim_mb_start, sim_mb_end, sim_mb_range,
     sim_bins,  SIM_BINSIZE,
     rep1_ins,  rep1_bnd_mb,
     colors["rep1"], "Simulation — 20U supernatant rep1 (25 kb)",
     False, "sim"),

    (4, 5,
     sens_img, sens_vmax, sens_w, sens_n,
     sim_mb_start, sim_mb_end, sim_mb_range,
     sim_bins,  SIM_BINSIZE,
     sens_ins,  sens_bnd_mb,
     colors["sens"], "Simulation — Sensitivity rep2 (25 kb)",
     True, "sim"),
]

hic_im_ref = None
sim_im_ref = None

for (map_row, ins_row,
     img, vmax, w, n,
     mb_start, mb_end, mb_range,
     bins, binsize,
     ins, bnd_mb,
     color, label, show_x, cbar_key) in panels:

    # --- Contact map ---
    ax_map = fig.add_subplot(gs[map_row])

    im = ax_map.imshow(img, cmap="Reds", vmin=0, vmax=vmax,
                       aspect="auto", interpolation="nearest",
                       origin="upper")
    if cbar_key == "hic" and hic_im_ref is None:
        hic_im_ref = im
    if cbar_key == "sim" and sim_im_ref is None:
        sim_im_ref = im

    # Boundary vertical lines spanning full map height
    px_list = bin_to_pixel(bnd_mb, mb_start, mb_range, w)
    for px in px_list:
        ax_map.axvline(px, color=color, lw=1.2,
                       linestyle="--", alpha=0.8, zorder=4)

    ax_map.set_xticks([])
    ax_map.set_yticks([])
    for spine in ax_map.spines.values():
        spine.set_visible(False)
    ax_map.text(0.005, 0.97, label,
                transform=ax_map.transAxes,
                fontsize=10, fontweight="bold",
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.2",
                          fc="white", ec="none", alpha=0.7))

    # --- Insulation score ---
    ax_ins = fig.add_subplot(gs[ins_row])
    mb_arr = bins / 1e6

    ax_ins.plot(mb_arr, ins, color=color, lw=1.5, zorder=2)
    ax_ins.axhline(0, color="grey", lw=0.7, linestyle="--", zorder=1)
    ax_ins.fill_between(mb_arr, ins, 0,
                        where=(ins < 0),
                        alpha=0.2, color=color, zorder=1)

    # Mark called boundaries with downward triangles on the track
    for b in bnd_mb:
        y_val = np.interp(b, mb_arr, ins)
        ax_ins.axvline(b, color=color, lw=1.0,
                       linestyle=":", alpha=0.6, zorder=3)
        ax_ins.plot(b, y_val, "v", color=color,
                    ms=6, zorder=5, markeredgewidth=0.5,
                    markeredgecolor="white")

    ax_ins.set_xlim(mb_start, mb_end)
    ax_ins.set_ylabel("Ins.\nscore", fontsize=8)
    ax_ins.spines["top"].set_visible(False)
    ax_ins.spines["right"].set_visible(False)
    ax_ins.tick_params(axis="y", labelsize=8)

    if show_x:
        set_mb_ticks(ax_ins, mb_start, mb_end, fontsize=9)
        ax_ins.set_xlabel("Genomic position (Mb)", fontsize=11)
    else:
        ax_ins.set_xticks([])

# --- Stats text in dedicated bottom row ---
ax_stats = fig.add_subplot(gs[6])
ax_stats.axis("off")

stats_lines = [
    (f"TAD boundary overlap (tolerance ±{BOUNDARY_TOL_KB} kb)  |  "
     f"Hi-C boundaries: {rep1_overlap['n_bnd1']}"),
    (f"20U supernatant rep1:  "
     f"F1 = {rep1_overlap['f1']:.2f}    "
     f"Recall = {rep1_overlap['recall']:.2f}    "
     f"Precision = {rep1_overlap['precision']:.2f}    "
     f"({rep1_overlap['n_matched']} / {rep1_overlap['n_bnd1']} Hi-C boundaries matched)"),
    (f"Sensitivity:           "
     f"F1 = {sens_overlap['f1']:.2f}    "
     f"Recall = {sens_overlap['recall']:.2f}    "
     f"Precision = {sens_overlap['precision']:.2f}    "
     f"({sens_overlap['n_matched']} / {sens_overlap['n_bnd1']} Hi-C boundaries matched)"),
]

for i, line in enumerate(stats_lines):
    ax_stats.text(
        0.0, 1.0 - i * 0.35, line,
        transform=ax_stats.transAxes,
        fontsize=9,
        fontweight="bold" if i == 0 else "normal",
        color="#222222" if i == 0 else "#444444",
        va="top", ha="left",
        family="monospace",
    )

# --- Separate colorbars ---
cbar_hic = fig.colorbar(hic_im_ref, cax=cbar_hic_ax)
cbar_hic.set_label("log(1 + freq.)\nHi-C", fontsize=8, labelpad=6)
cbar_hic.ax.tick_params(labelsize=7)

cbar_sim = fig.colorbar(sim_im_ref, cax=cbar_sim_ax)
cbar_sim.set_label("Contact prob.\nSimulation", fontsize=8, labelpad=6)
cbar_sim.ax.tick_params(labelsize=7)

# --- Legend for boundary markers ---
legend_elements = [
    Line2D([0], [0], color=colors["hic"],  lw=1.5,
           linestyle="--", label="Hi-C boundaries"),
    Line2D([0], [0], color=colors["rep1"], lw=1.5,
           linestyle="--", label="20U sup. rep1 boundaries"),
    Line2D([0], [0], color=colors["sens"], lw=1.5,
           linestyle="--", label="Sensitivity boundaries"),
]
fig.legend(handles=legend_elements,
           loc="upper right",
           bbox_to_anchor=(0.87, 0.97),
           fontsize=8.5, frameon=True,
           framealpha=0.9, edgecolor="#cccccc")

fig.suptitle(
    f"{CHROM}:{START_BP//1_000_000}–{END_BP//1_000_000} Mb  |  "
    f"TAD boundary comparison",
    fontsize=13, fontweight="bold", x=0.46,
)

fig.savefig(OUT_COMPARE, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUT_COMPARE}")
print("Done.")
