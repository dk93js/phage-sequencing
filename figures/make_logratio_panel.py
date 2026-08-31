"""
Quantitative enrichment panel: distribution of log10(ssDNA/dsDNA) at R3, split by
ELISA status. House style (Arial; 13 / title 14 / label 14 / tick 12 / legend 10;
3.6-in panel height). ELISA colors match Fig 3 scatter (purple +, grey -).
Writes into version 3/panels/. Run with /opt/anaconda3/bin/python.
"""
import csv, collections
from pathlib import Path
import numpy as np
import openpyxl
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
PANELS = ROOT / "out/figures"
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "savefig.dpi": 300, "savefig.bbox": "tight",
})
C_POS, C_NEG = "#6A1B9A", "#9aa0a6"   # ELISA + / - (match Fig 3 scatter)

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
R = 3
mapped = [s for s in y if s in ss and s in ds and ds[s][R] > 0]
lp = [np.log10(ss[s][R] / ds[s][R]) for s in mapped if y[s] == 1]
ln = [np.log10(ss[s][R] / ds[s][R]) for s in mapped if y[s] == 0]
pval = mannwhitneyu(lp, ln).pvalue

fig, ax = plt.subplots(figsize=(3.6, 3.6))
bins = np.linspace(-2, 2, 21)
ax.hist(ln, bins=bins, density=True, color=C_NEG, alpha=0.6, label="ELISA −")
ax.hist(lp, bins=bins, density=True, color=C_POS, alpha=0.7, label="ELISA +")
ax.axvline(0, color="k", ls=":", lw=1)
# medians drawn dotted (JC): solid vertical lines read as histogram bars
ax.axvline(np.median(ln), color=C_NEG, ls=(0, (1, 1.6)), lw=2.4)
ax.axvline(np.median(lp), color=C_POS, ls=(0, (1, 1.6)), lw=2.4)
ax.set_ylim(0, ax.get_ylim()[1] * 1.35)   # headroom so legend clears the bars
ax.set_xlabel("log$_{10}$(ssDNA / dsDNA) at R3")
ax.set_ylabel("density")
ax.set_title("binders shift toward the ssDNA pool")
ax.text(0.98, 0.97,
        f"median + : {np.median(lp):+.2f}\nmedian − : {np.median(ln):+.2f}\nP = {pval:.0e}",
        transform=ax.transAxes, ha="right", va="top", fontsize=9)
ax.text(0.02, 0.45, "ds ←", transform=ax.transAxes, fontsize=9, color="#666")
ax.text(0.90, 0.45, "→ ss", transform=ax.transAxes, fontsize=9, color="#666")
ax.legend(frameon=False, loc="upper left", fontsize=9)
fig.tight_layout()
out = PANELS / "logratio_by_elisa_R3_sized.png"
fig.savefig(out)
plt.close(fig)
print("wrote", out)
print(f"median + {np.median(lp):+.2f} / - {np.median(ln):+.2f}  MWU p={pval:.2e}  n={len(lp)}+/{len(ln)}-")
