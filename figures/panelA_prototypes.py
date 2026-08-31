"""
Prototypes to replace Figure 1 Panel A (the metrics table) with a graph.

Values are taken verbatim from the curated 1A table so the graph matches what
has already been vetted (titer / Shannon H / Top-1%, ss & ds, R0-R4).

Outputs (panels/):
  panelA_v1_diversity_twinaxis.png   - Shannon H (left) + Top-1% (right), ss/ds
  panelA_v2_combined_titer.png       - titer (log, top) + diversity twin-axis (bottom)
  panelA_v3_smallmultiples.png       - titer | Shannon H | Top-1% (1x3)
"""
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec

CODE_DIR = Path(__file__).resolve().parent
PANELS = CODE_DIR.parent / "panels"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "lines.linewidth": 2.0, "lines.markersize": 6,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

R = ["R0", "R1", "R2", "R3", "R4"]
x = np.arange(5)
titer = np.array([np.nan, 2.6e-4, 2.7e-7, 4.5e-6, 4.6e-5])
shannon_ss = np.array([6.12, 6.12, 6.04, 3.65, 3.96])
shannon_ds = np.array([6.17, 6.15, 6.03, 4.88, 5.35])
# Top-1% = cumulative read fraction of the top 1% of clonotypes (the ~11 most-abundant
# per library) within the analysis set; same denominator basis as Shannon H.
top1_ss = np.array([19.4, 19.6, 23.7, 68.3, 62.2])
top1_ds = np.array([19.1, 19.3, 24.4, 44.7, 37.5])
# "all" = combined ssDNA+dsDNA library, shared CDRH3 merged (reads summed),
# computed by the same method as ss/ds (raw reads, natural-log Shannon, top-1%-of-clones frac).
shannon_all = np.array([6.16, 6.15, 6.05, 4.38, 4.69])
top1_all = np.array([19.9, 20.1, 24.5, 57.1, 51.4])

C_SH = "#3D5A80"     # Shannon (slate blue)
C_T1 = "#E07A1F"     # Top-1% (orange)
C_TI = "#222222"     # titer (near-black)
C_ALL = "#000000"    # combined "all" library (matches Fig 2 "All" = black)


def _diversity_twin(ax_left):
    axR = ax_left.twinx()
    axR.spines["top"].set_visible(False)
    # Shannon H on the left axis
    l1, = ax_left.plot(x, shannon_ss, "o-", color=C_SH, label="Shannon H (ss)")
    l2, = ax_left.plot(x, shannon_ds, "s--", color=C_SH, mfc="white",
                       label="Shannon H (ds)")
    # Top-1% on the right axis
    l3, = axR.plot(x, top1_ss, "o-", color=C_T1, label="Top-1% (ss)")
    l4, = axR.plot(x, top1_ds, "s--", color=C_T1, mfc="white", label="Top-1% (ds)")
    ax_left.set_ylabel("Shannon H", color=C_SH)
    axR.set_ylabel("Top-1% of reads (%)", color=C_T1)
    ax_left.tick_params(axis="y", colors=C_SH)
    axR.tick_params(axis="y", colors=C_T1)
    ax_left.set_ylim(3, 7)
    axR.set_ylim(0, 40)
    ax_left.set_xticks(x)
    ax_left.set_xticklabels(R)
    return [l1, l2, l3, l4], axR


# ---- v1: diversity twin-axis only --------------------------------------
fig, ax = plt.subplots(figsize=(5.2, 4.6))
lines, axR = _diversity_twin(ax)
ax.set_xlabel("selection round")
ax.legend(lines, [l.get_label() for l in lines], frameon=False,
          loc="center left", fontsize=9.5)
ax.set_title("Diversity collapses at R3", fontweight="bold")
fig.tight_layout()
fig.savefig(PANELS / "panelA_v1_diversity_twinaxis.png")
plt.close(fig)

# ---- v2: titer (top) + diversity twin-axis (bottom) --------------------
fig = plt.figure(figsize=(5.2, 6.0))
gs = gridspec.GridSpec(2, 1, height_ratios=[1.0, 2.4], hspace=0.12)
axT = fig.add_subplot(gs[0])
axT.plot(x, titer, "D-", color=C_TI)
axT.set_yscale("log")
axT.set_ylabel("output/input\n(titer)")
axT.set_xticks(x)
axT.set_xticklabels([])
axT.set_xlim(-0.4, 4.4)
axB = fig.add_subplot(gs[1])
lines, axR = _diversity_twin(axB)
axB.set_xlim(-0.4, 4.4)
axB.set_xlabel("selection round")
axB.legend(lines, [l.get_label() for l in lines], frameon=False,
           loc="center left", fontsize=9)
fig.savefig(PANELS / "panelA_v2_combined_titer.png")
plt.close(fig)

# ---- v3: small multiples (titer | Shannon | Top-1%) --------------------
fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.6))
a0, a1, a2 = axes
a0.plot(x, titer, "D-", color=C_TI)
a0.set_yscale("log")
a0.set_title("output/input (titer)")
a0.set_ylabel("ratio")
a1.plot(x, shannon_ss, "o-", color=C_SH, label="ssDNA")
a1.plot(x, shannon_ds, "s--", color=C_SH, mfc="white", label="dsDNA")
a1.plot(x, shannon_all, "^-", color=C_ALL, label="all")
a1.set_title("clonal diversity")
a1.set_ylabel("Shannon H")
a1.set_ylim(3, 7)
a1.legend(frameon=False)
a2.plot(x, top1_ss, "o-", color=C_T1, label="ssDNA")
a2.plot(x, top1_ds, "s--", color=C_T1, mfc="white", label="dsDNA")
a2.plot(x, top1_all, "^-", color=C_ALL, label="all")
a2.set_title("Top-1% of reads (%)")
a2.set_ylabel("% of reads")
a2.set_ylim(0, 75)
a2.legend(frameon=False)
for a in axes:
    a.set_xticks(x)
    a.set_xticklabels(R)
    a.set_xlabel("selection round")
fig.tight_layout()
fig.savefig(PANELS / "panelA_v3_smallmultiples.png")
plt.close(fig)

print("wrote panelA_v1_diversity_twinaxis.png, panelA_v2_combined_titer.png, panelA_v3_smallmultiples.png")
