"""
Supplementary figure for Fig 3: HCDR3 physicochemical properties barely distinguish
ELISA+ from ELISA- clones (descriptive, NO machine learning). All 568 ELISA HCDR3.
6 panels: length, cysteine, hydrophobicity (GRAVY), net charge, aromaticity, pI.
House style. Run with /opt/anaconda3/bin/python.
"""
import collections
from pathlib import Path
import numpy as np
import openpyxl
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "out/figures"
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 13.5,
    "axes.labelsize": 13.5, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "savefig.dpi": 300, "savefig.bbox": "tight",
})
C_POS, C_NEG = "#6A1B9A", "#9aa0a6"
KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
      'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
      'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

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
P = [s for s in y if y[s] == 1]; N = [s for s in y if y[s] == 0]


def gravy(s):
    return float(np.mean([KD.get(a, 0) for a in s]))
def charge(s):
    return sum(s.count(a) for a in "KR") - sum(s.count(a) for a in "DE")
def arom(s):
    return 100 * sum(s.count(a) for a in "FWY") / len(s)


def violin(ax, fn, ylab):
    vn = [fn(s) for s in N]; vp = [fn(s) for s in P]
    parts = ax.violinplot([vn, vp], positions=[0, 1], showmedians=True, widths=0.75)
    for pc, c in zip(parts["bodies"], [C_NEG, C_POS]):
        pc.set_facecolor(c); pc.set_alpha(0.6)
    for k in ("cbars", "cmins", "cmaxes", "cmedians"):
        parts[k].set_color("#444444"); parts[k].set_linewidth(1.0)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["ELISA −", "ELISA +"])
    ax.set_ylabel(ylab)
    ax.set_title(f"p = {mannwhitneyu(vp, vn).pvalue:.2f}")


fig, axes = plt.subplots(2, 3, figsize=(10.8, 6.6))

# A length (kept)
violin(axes[0, 0], len, "HCDR3 length (aa)")
# B cysteine (kept as grouped bar)
ax = axes[0, 1]
def cfrac(grp):
    cc = [min(s.count("C"), 2) for s in grp]
    return [100 * cc.count(k) / len(cc) for k in (0, 1, 2)]
fp, fn_ = cfrac(P), cfrac(N)
x = np.arange(3); w = 0.38
ax.bar(x - w / 2, fn_, w, color=C_NEG, label="ELISA −")
ax.bar(x + w / 2, fp, w, color=C_POS, label="ELISA +")
ax.set_xticks(x); ax.set_xticklabels(["0", "1", "≥2"]); ax.set_xlabel("cysteines in HCDR3")
ax.set_ylabel("% of clones")
pc = mannwhitneyu([s.count("C") for s in P], [s.count("C") for s in N]).pvalue
ax.set_title(f"p = {pc:.2f}"); ax.legend(frameon=False, fontsize=9)
# C-F physicochemical
violin(axes[0, 2], gravy, "hydrophobicity (GRAVY)")
violin(axes[1, 0], charge, "net charge (pH 7)")

# E, F: low-dim embeddings (UMAP) colored by ELISA -- physicochemical (Atchley) + protein LM (ESM-2)
def umap_panel(ax, npz, title):
    d = np.load(OUT / npz)
    U, ylm, auc = d["U"], d["y"], float(d["auc"])
    ax.scatter(U[ylm == 0, 0], U[ylm == 0, 1], s=10, c=C_NEG, alpha=0.6, label="ELISA −")
    ax.scatter(U[ylm == 1, 0], U[ylm == 1, 1], s=22, c=C_POS, edgecolors="k", linewidths=0.3, label="ELISA +")
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(f"{title} (AUC {auc:.2f})")

umap_panel(axes[1, 1], "atchley_umap_coords.npz", "Atchley-factor embedding")
umap_panel(axes[1, 2], "esm_umap_coords.npz", "ESM-2 embedding")

fig.tight_layout()
fn_out = OUT / "supp_sequence_sized.png"
fig.savefig(fn_out); plt.close(fig)
print("wrote", fn_out.name)
for name, f in [("length", len), ("cys", lambda s: s.count("C")), ("gravy", gravy),
                ("charge", charge), ("arom", arom)]:
    print(f"  {name:8} p={mannwhitneyu([f(s) for s in P], [f(s) for s in N]).pvalue:.3f}")
