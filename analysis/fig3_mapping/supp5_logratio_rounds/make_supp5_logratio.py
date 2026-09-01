"""
New Supplementary figure: log10(ssDNA/dsDNA) of antigen-reactive vs non-reactive HCDR3
clonotypes at EVERY round (R0-R4). Fig 3B shows R3 only; this puts all five rounds on the
record so the round choice is not a selection.

Same 146 assayed clonotypes (35 ELISA+, 111 ELISA-) at every round -- no per-round dropout.
House style: Arial, 13 / title 14 / label 14 / tick 12 / legend 10, 3.6-in panels.
Numbered as Supplementary Figure 5 (appended; no renumbering of Supp 1-4).

Outputs (this folder):
  suppfig4_logratio_R0..R4.png   - five per-round histograms (Supp 2 layout convention)
  suppfig4_logratio_summary.png  - single panel: median (IQR) by round, both groups
  preview_*.png
Run with /opt/anaconda3/bin/python.
"""
import csv, collections
from pathlib import Path
import numpy as np
import openpyxl
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.image import imread
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "analysis/fig3_mapping/supp5_logratio_rounds"
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "savefig.dpi": 300, "savefig.bbox": "tight",
})
C_POS, C_NEG = "#6A1B9A", "#9aa0a6"
BINS = np.linspace(-2, 2, 21)

wb = openpyxl.load_workbook(ROOT / "data/phage_elisa_wells.xlsx",
                            read_only=True, data_only=True)
raw = collections.defaultdict(list)
for r in list(wb["ELISA"].iter_rows(values_only=True))[1:]:
    if r[2] is None:
        continue
    try:
        p = int(r[5])
    except (TypeError, ValueError):
        continue
    raw[str(r[2]).strip()].append(p)
y = {s: (1 if sum(v) / len(v) >= 0.5 else 0) for s, v in raw.items() if v}
tj = collections.defaultdict(dict)
with open(ROOT / "data/report_q40_ppm100_list.csv") as f:
    for row in csv.DictReader(f):
        tj[row["lib"]][row["seq"]] = [float(row[f"PPM{i}"]) for i in range(5)]
ss, ds = tj["ss"], tj["ds"]

S = {}
for R in range(5):
    mapped = [s for s in y if s in ss and s in ds and ds[s][R] > 0]
    lp = np.array([np.log10(ss[s][R] / ds[s][R]) for s in mapped if y[s] == 1])
    ln = np.array([np.log10(ss[s][R] / ds[s][R]) for s in mapped if y[s] == 0])
    S[R] = (lp, ln, mannwhitneyu(lp, ln).pvalue)

print(f"{'R':>2} {'n+':>4} {'n-':>4} {'med+':>7} {'med-':>7} {'diff':>7} {'fold+':>6} {'P':>10}")
for R in range(5):
    lp, ln, p = S[R]
    print(f"{R:>2} {len(lp):>4} {len(ln):>4} {np.median(lp):>+7.3f} {np.median(ln):>+7.3f} "
          f"{np.median(lp)-np.median(ln):>+7.3f} {10**np.median(lp):>6.2f} {p:>10.2e}")

# ---------- per-round histograms ----------
peak = 0
for R in range(5):
    lp, ln, _ = S[R]
    for v in (lp, ln):
        h, _ = np.histogram(v, bins=BINS, density=True)
        peak = max(peak, h.max())
print(f"\nmax density across all rounds = {peak:.1f} (shared y would crush R3/R4)")

def hist_panel(R, ylabel=True, ymax=None):
    lp, ln, pval = S[R]
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    ax.hist(ln, bins=BINS, density=True, color=C_NEG, alpha=0.6, label="ELISA −")
    ax.hist(lp, bins=BINS, density=True, color=C_POS, alpha=0.7, label="ELISA +")
    ax.axvline(0, color="k", ls=":", lw=1)
    ax.axvline(np.median(ln), color=C_NEG, ls=(0, (1, 1.6)), lw=2.4)
    ax.axvline(np.median(lp), color=C_POS, ls=(0, (1, 1.6)), lw=2.4)
    ax.set_ylim(0, ymax if ymax else ax.get_ylim()[1] * 1.38)
    ax.set_xlabel(f"log$_{{10}}$(ssDNA / dsDNA) at R{R}")
    if ylabel:
        ax.set_ylabel("density")
    ax.set_title(f"R{R}")
    ax.text(0.98, 0.97,
            f"median + : {np.median(lp):+.2f}\nmedian − : {np.median(ln):+.2f}\nP = {pval:.0e}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    fn = OUT / f"suppfig4_logratio_R{R}.png"
    fig.savefig(fn); plt.close(fig)
    return fn

SHARED = peak * 1.38   # peak density is ~4.5 in every round (early rounds are one
                       # narrow spike), so a shared y-axis is affordable and honest
files = [hist_panel(R, ylabel=True, ymax=SHARED) for R in range(5)]
print(f"shared y-limit = {SHARED:.2f}")

# ---------- summary panel: median (IQR) by round ----------
fig, ax = plt.subplots(figsize=(3.6, 3.6))
x = np.arange(5)
for lab, idx, col, off in (("ELISA +", 0, C_POS, -0.07), ("ELISA −", 1, C_NEG, +0.07)):
    med = np.array([np.median(S[R][idx]) for R in range(5)])
    q1 = np.array([np.percentile(S[R][idx], 25) for R in range(5)])
    q3 = np.array([np.percentile(S[R][idx], 75) for R in range(5)])
    ax.errorbar(x + off, med, yerr=[med - q1, q3 - med], color=col, marker="o",
                ms=6, lw=2.0, capsize=3, label=lab)
ax.axhline(0, color="k", ls=":", lw=1)
ax.set_xticks(x); ax.set_xticklabels([f"R{i}" for i in range(5)])
ax.set_xlabel("panning round")
ax.set_ylabel("log$_{10}$(ssDNA / dsDNA)")
ax.set_ylim(-0.85, 1.35)
for R in range(5):
    p = S[R][2]
    ax.text(R, 1.22, f"{p:.0e}".replace("e-0", "e−"), ha="center", fontsize=8, color="#444")
ax.text(-0.62, 1.22, "P", ha="center", fontsize=8, color="#444")
ax.legend(frameon=False, loc="lower right", fontsize=9)
fig.tight_layout()
fn_sum = OUT / "suppfig4_logratio_summary.png"
fig.savefig(fn_sum); plt.close(fig)

# ---------- previews ----------
fig, axes = plt.subplots(1, 5, figsize=(21, 4.4))
for ax, fn in zip(axes, files):
    ax.imshow(imread(fn)); ax.axis("off")
fig.tight_layout(); fig.savefig(OUT / "preview_supp4_histograms.png", dpi=110); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.2, 4.4))
ax.imshow(imread(fn_sum)); ax.axis("off")
fig.tight_layout(); fig.savefig(OUT / "preview_supp4_summary.png", dpi=160); plt.close(fig)
print("\nwrote:", *[f.name for f in files], fn_sum.name, sep="\n  ")
