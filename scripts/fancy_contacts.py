#!/usr/bin/env python3
"""
Quantitative comparison of simulation vs Hi-C contact maps.
- Normalised P(s) curves
- O/E correlation statistics
- Rank-normalised difference maps (publication quality)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
from scipy.ndimage import rotate
from scipy.stats import pearsonr, spearmanr, rankdata
from pathlib import Path
import cooler

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MCOOL_PATH          = Path("/home/gdbishoplab/jxb377/condense-seq/GSE246947_hg38_MCF10A_WT_merged.mcool")
SUPERNATANT_REP1_NPY = Path("/home/gdbishoplab/jxb377/condense-seq/results_supernatant_rep1_chr2_25kb_gpu/contact_matrix_avg.npy")
SENS_NPY            = Path("/home/gdbishoplab/jxb377/condense-seq/results_rep2_chr2_25kb_gpu/contact_matrix_avg.npy")
OUT_PS              = Path("/home/gdbishoplab/jxb377/condense-seq/ps_curves_norm.png")
OUT_DIFF            = Path("/home/gdbishoplab/jxb377/condense-seq/difference_maps_pub.png")
OUT_STATS           = Path("/home/gdbishoplab/jxb377/condense-seq/correlation_stats_norm.txt")

CHROM       = "chr2"
START_BP    = 10_000_000
END_BP      = 14_000_000
MCOOL_RES   = 20_000
SIM_BINSIZE = 25_000
TRIM_FRAC   = 0.35

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def symmetrise(mat):
    c = mat.copy().astype(float)
    c = np.triu(c) + np.triu(c, 1).T
    c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(c, 0.0)
    return c


def oe_normalise(mat):
    """Divide each diagonal by its mean — removes distance decay."""
    n  = mat.shape[0]
    oe = np.zeros_like(mat, dtype=float)
    for d in range(1, n):
        diag = np.diag(mat, k=d)
        pos  = diag[diag > 0]
        mu   = pos.mean() if len(pos) > 0 else 1.0
        idx  = np.arange(n - d)
        val  = diag / mu
        oe[idx,     idx + d] = val
        oe[idx + d, idx    ] = val
    return oe


def rank_normalise(mat):
    """
    Convert upper-triangle values to percentile ranks in [0, 1].
    Zeros (no contact) are kept as zero — only positive values are ranked.
    This makes two matrices directly comparable regardless of scale.
    """
    n    = mat.shape[0]
    out  = np.zeros_like(mat, dtype=float)
    idx  = np.triu_indices(n, k=1)
    vals = mat[idx]

    pos_mask         = vals > 0
    ranked           = np.zeros_like(vals)
    ranked[pos_mask] = rankdata(vals[pos_mask]) / pos_mask.sum()

    out[idx]            = ranked
    out[idx[1], idx[0]] = ranked   # symmetrise
    return out


def ps_curve_normalised(mat, binsize):
    """Mean contact at each separation, normalised so first point = 1."""
    n   = mat.shape[0]
    sep, avg = [], []
    for d in range(1, n):
        diag = np.diag(mat, k=d)
        pos  = diag[diag > 0]
        if len(pos) > 0:
            sep.append(d * binsize)
            avg.append(pos.mean())
    sep = np.array(sep, dtype=float)
    avg = np.array(avg, dtype=float)
    if avg[0] > 0:
        avg = avg / avg[0]
    return sep / 1e6, avg


def interpolate_to_common_bins(mat_sim, sim_binsize, n_hic, hic_binsize):
    """Nearest-neighbour resample simulation matrix to Hi-C bin size."""
    n_sim = mat_sim.shape[0]
    ratio = sim_binsize / hic_binsize
    out   = np.zeros((n_hic, n_hic), dtype=float)
    for i in range(n_hic):
        for j in range(n_hic):
            si = min(int(i / ratio), n_sim - 1)
            sj = min(int(j / ratio), n_sim - 1)
            out[i, j] = mat_sim[si, sj]
    return out


def scc_score(mat1, mat2, max_diag=None):
    """Stratum-adjusted correlation coefficient."""
    n = min(mat1.shape[0], mat2.shape[0])
    if max_diag is None:
        max_diag = n - 1
    rhos, weights = [], []
    for d in range(1, max_diag + 1):
        d1    = np.diag(mat1, k=d)[:n - d]
        d2    = np.diag(mat2, k=d)[:n - d]
        valid = np.isfinite(d1) & np.isfinite(d2) & ((d1 > 0) | (d2 > 0))
        if valid.sum() < 3:
            continue
        r, _ = pearsonr(d1[valid], d2[valid])
        if np.isnan(r):
            continue
        rhos.append(r)
        weights.append(valid.sum())
    if not rhos:
        return np.nan, [], []
    return np.average(rhos, weights=weights), rhos, weights


def prepare_diff_triangle(diff_mat, trim_frac=TRIM_FRAC):
    """Rotate, trim and mask a difference matrix for triangle display."""
    n      = diff_mat.shape[0]
    upper  = np.triu(diff_mat, 1)
    vals   = upper[np.abs(upper) > 0]
    vmax   = np.percentile(np.abs(vals), 95) if len(vals) > 0 else 1.0

    rotated  = rotate(upper, angle=45, reshape=True, cval=np.nan, order=1)
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
    img        = np.where(outside, np.nan, trimmed)
    return np.ma.array(img, mask=np.isnan(img)), vmax, w, n


def add_ticks(ax, w, n, bin_starts, binsize, fontsize=10):
    mb      = bin_starts / 1e6
    mb_end  = mb[-1] + binsize / 1e6
    tick_mb = np.linspace(mb[0], mb_end, 7)
    tidx    = [(np.abs(mb - t)).argmin() for t in tick_mb]
    tidx    = sorted(set(np.clip(tidx, 0, n - 1)))
    ax.set_xticks([int(i * w / n) for i in tidx])
    ax.set_xticklabels([f"{mb[i]:.1f}" for i in tidx],
                       fontsize=fontsize, rotation=45, ha="right")


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

print("Loading simulations...")
rep1_mat = symmetrise(np.load(SUPERNATANT_REP1_NPY))
sens_mat = symmetrise(np.load(SENS_NPY))

print("Resampling simulations to Hi-C resolution...")
rep1_res = interpolate_to_common_bins(rep1_mat, SIM_BINSIZE, n_hic, MCOOL_RES)
sens_res = interpolate_to_common_bins(sens_mat, SIM_BINSIZE, n_hic, MCOOL_RES)

# ---------------------------------------------------------------------------
# 1. Normalised P(s) curves
# ---------------------------------------------------------------------------
print("Computing P(s) curves...")
hic_sep,  hic_ps  = ps_curve_normalised(hic_mat,  MCOOL_RES)
rep1_sep, rep1_ps = ps_curve_normalised(rep1_mat, SIM_BINSIZE)
sens_sep, sens_ps = ps_curve_normalised(sens_mat, SIM_BINSIZE)

plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.size":          11,
    "axes.linewidth":     1.2,
    "xtick.major.width":  1.2,
    "ytick.major.width":  1.2,
})

colors = {"hic": "#222222", "rep1": "#e41a1c", "sens": "#377eb8"}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, logscale in zip(axes, [False, True]):
    ax.plot(hic_sep,  hic_ps,  color=colors["hic"],
            lw=2,   label="MCF10A WT Hi-C")
    ax.plot(rep1_sep, rep1_ps, color=colors["rep1"],
            lw=2,   linestyle="--", label="Simulation — 20U supernatant rep1")
    ax.plot(sens_sep, sens_ps, color=colors["sens"],
            lw=2,   linestyle=":",  label="Simulation — Sensitivity (rep2)")
    ax.set_xlabel("Genomic separation (Mb)", fontsize=12)
    ax.set_ylabel("Normalised contact probability",  fontsize=12)
    ax.legend(fontsize=9, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if logscale:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title("P(s) — log-log scale", fontsize=12)
    else:
        ax.axvspan(0, 2, alpha=0.05, color="grey", zorder=0)
        ax.set_xlim(left=0)
        ax.set_title("P(s) — linear scale", fontsize=12)

fig.suptitle(
    f"Normalised contact probability vs separation — "
    f"{CHROM}:{START_BP//1_000_000}–{END_BP//1_000_000} Mb",
    fontsize=13, fontweight="bold"
)
fig.tight_layout()
fig.savefig(OUT_PS, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUT_PS}")

# ---------------------------------------------------------------------------
# 2. Correlation statistics — O/E matrices
# ---------------------------------------------------------------------------
print("Computing O/E correlations...")
hic_oe  = oe_normalise(hic_mat)
rep1_oe = oe_normalise(rep1_res)
sens_oe = oe_normalise(sens_res)

def upper_vals(mat):
    idx = np.triu_indices(mat.shape[0], k=1)
    v   = mat[idx]
    return v[np.isfinite(v) & (v != 0)]

hic_v  = upper_vals(hic_oe)
rep1_v = upper_vals(rep1_oe)
sen_v  = upper_vals(sens_oe)
ml     = min(len(hic_v), len(rep1_v), len(sen_v))
hic_v, rep1_v, sen_v = hic_v[:ml], rep1_v[:ml], sen_v[:ml]

pear_rep1, _ = pearsonr(hic_v, rep1_v)
pear_sen,  _ = pearsonr(hic_v, sen_v)
spear_rep1   = spearmanr(hic_v, rep1_v).statistic
spear_sen    = spearmanr(hic_v, sen_v).statistic

max_diag       = int(2e6 / MCOOL_RES)
scc_rep1, _, _ = scc_score(hic_oe, rep1_oe, max_diag=max_diag)
scc_sen,  _, _ = scc_score(hic_oe, sens_oe, max_diag=max_diag)

# Also compute rank-normalised Pearson/Spearman for reporting
hic_rn  = rank_normalise(hic_mat)
rep1_rn = rank_normalise(rep1_res)
sen_rn  = rank_normalise(sens_res)

hic_rv  = upper_vals(hic_rn)
rep1_rv = upper_vals(rep1_rn)
sen_rv  = upper_vals(sen_rn)
ml2     = min(len(hic_rv), len(rep1_rv), len(sen_rv))
hic_rv, rep1_rv, sen_rv = hic_rv[:ml2], rep1_rv[:ml2], sen_rv[:ml2]

pear_rep1_rn,  _ = pearsonr(hic_rv, rep1_rv)
pear_sen_rn,   _ = pearsonr(hic_rv, sen_rv)
spear_rep1_rn    = spearmanr(hic_rv, rep1_rv).statistic
spear_sen_rn     = spearmanr(hic_rv, sen_rv).statistic

with open(OUT_STATS, "w") as f:
    f.write(f"Correlation statistics — {CHROM}:{START_BP}-{END_BP}\n")
    f.write(f"Hi-C: {MCOOL_RES} bp  |  Simulation: {SIM_BINSIZE} bp\n\n")
    f.write(f"{'Metric':<35} {'20U sup. rep1':>14} {'Sensitivity rep2':>18}\n")
    f.write("-" * 69 + "\n")
    f.write(f"{'SCC — O/E (up to 2 Mb)':<35} {scc_rep1:>14.4f} {scc_sen:>18.4f}\n")
    f.write(f"{'Pearson r — O/E':<35} {pear_rep1:>14.4f} {pear_sen:>18.4f}\n")
    f.write(f"{'Spearman rho — O/E':<35} {spear_rep1:>14.4f} {spear_sen:>18.4f}\n")
    f.write(f"{'Pearson r — rank normalised':<35} {pear_rep1_rn:>14.4f} {pear_sen_rn:>18.4f}\n")
    f.write(f"{'Spearman rho — rank normalised':<35} {spear_rep1_rn:>14.4f} {spear_sen_rn:>18.4f}\n")

print(open(OUT_STATS).read())

# ---------------------------------------------------------------------------
# 3. Publication-quality difference maps — rank normalised
# ---------------------------------------------------------------------------
print("Computing rank-normalised difference maps...")

diff_rep1 = np.triu(rep1_rn - hic_rn, 1)
diff_sen  = np.triu(sen_rn  - hic_rn, 1)

rep1_tri, rep1_vmax, rep1_w, rep1_n = prepare_diff_triangle(diff_rep1)
sens_tri, sens_vmax, sens_w, sens_n = prepare_diff_triangle(diff_sen)
vmax_shared = max(rep1_vmax, sens_vmax)

fig = plt.figure(figsize=(14, 9))

# Layout: two rows, generous top margin for titles + stats
gs = GridSpec(2, 1, figure=fig,
              hspace=0.55,
              left=0.03, right=0.87,
              top=0.91,  bottom=0.09)

# Single shared colorbar axis
cbar_ax = fig.add_axes([0.895, 0.15, 0.018, 0.65])

panels = [
    (gs[0], rep1_tri, rep1_w, rep1_n,
     "20U supernatant rep1 − Hi-C",
     scc_rep1, pear_rep1_rn, spear_rep1_rn,
     False),
    (gs[1], sens_tri, sens_w, sens_n,
     "Sensitivity (rep2) − Hi-C",
     scc_sen, pear_sen_rn, spear_sen_rn,
     True),
]

im_ref = None
for spec, tri, w, n, ptitle, scc, pear, spear, show_x in panels:
    ax = fig.add_subplot(spec)

    im = ax.imshow(tri, cmap="RdBu_r",
                   vmin=-vmax_shared, vmax=vmax_shared,
                   aspect="auto", interpolation="nearest", origin="upper")
    im_ref = im

    if show_x:
        add_ticks(ax, w, n, hic_bins, MCOOL_RES, fontsize=10)
        ax.set_xlabel("Genomic position (Mb)", fontsize=11)
    else:
        ax.set_xticks([])

    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Bold panel title — sits inside axes top-left
    ax.text(0.0, 1.13, ptitle,
            transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="bottom", ha="left")

    # Stats line — sits just below panel title
    stats_str = (f"SCC (O/E) = {scc:.3f}     "
                 f"Pearson r = {pear:.3f}     "
                 f"Spearman ρ = {spear:.3f}   "
                 f"[rank-normalised]")
    ax.text(0.0, 1.05, stats_str,
            transform=ax.transAxes,
            fontsize=8.5, color="#444444", va="bottom", ha="left")

cbar = fig.colorbar(im_ref, cax=cbar_ax)
cbar.set_label("Rank-normalised difference\n(simulation − Hi-C)",
               fontsize=10, labelpad=10)
cbar.ax.tick_params(labelsize=9)

fig.suptitle(
    f"{CHROM}:{START_BP//1_000_000}–{END_BP//1_000_000} Mb  |  "
    f"Simulation vs MCF10A WT Hi-C",
    fontsize=13, fontweight="bold", x=0.45, y=0.97
)

fig.savefig(OUT_DIFF, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUT_DIFF}")
print("Done.")
