"""
Fig 3 (mapping) in the locked house style: Arial, font 13 / title 14 non-bold /
label 14 / tick 12 / legend 10, 3.6-in square panel scale, lowercase rule, red y=x.
Panels: 3A scatter ss-vs-ds with ELISA +/- overlay (R3 and R4); 3C worked-clone table.
ELISA ground truth = Yoo 2020/180528_TR_MET...xlsx (sheet ELISA). Run with /opt/anaconda3/bin/python.
"""
import csv, collections
from pathlib import Path
import numpy as np
import openpyxl
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "out/figures"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "lines.linewidth": 2.0, "lines.markersize": 6,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})
PANEL_H = 3.6
C_POS = "#6A1B9A"   # ELISA+ : purple (distinct from cluster red/blue/green & library orange/slate)
C_NEG = "#7f7f7f"   # ELISA- : grey open
C_BG = "#dcdcdc"    # all shared clones

# ---- load ELISA labels (majority per HCDR3) ----
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

# ---- load ss/ds trajectories ----
tj = collections.defaultdict(dict)
with open(ROOT / "data/report_q40_ppm100_list.csv") as f:
    for row in csv.DictReader(f):
        tj[row["lib"]][row["seq"]] = [float(row[f"PPM{i}"]) for i in range(5)]
ss, ds = tj["ss"], tj["ds"]
shared = [s for s in ss if s in ds]
mapped = [s for s in y if s in ss and s in ds]


def render_scatter(r):
    fig, ax = plt.subplots(figsize=(PANEL_H, PANEL_H))
    gx = np.array([ss[s][r] for s in shared]); gy = np.array([ds[s][r] for s in shared])
    m = (gx > 0) & (gy > 0)
    ax.scatter(gx[m], gy[m], s=7, color=C_BG, edgecolor="none", label="all shared clones", zorder=1)
    pos = [s for s in mapped if y[s] == 1 and ss[s][r] > 0 and ds[s][r] > 0]
    neg = [s for s in mapped if y[s] == 0 and ss[s][r] > 0 and ds[s][r] > 0]
    ax.scatter([ss[s][r] for s in neg], [ds[s][r] for s in neg], s=24, facecolors="none",
               edgecolors=C_NEG, linewidths=1.0, label=f"ELISA − (n={len(neg)})", zorder=3)
    ax.scatter([ss[s][r] for s in pos], [ds[s][r] for s in pos], s=40, color=C_POS,
               edgecolors="black", linewidths=0.4, label=f"ELISA + (n={len(pos)})", zorder=5)
    lo = min(gx[m].min(), gy[m].min()) * 0.7
    hi = max(gx[m].max(), gy[m].max()) * 1.4
    ax.plot([lo, hi], [lo, hi], color="red", linewidth=1.1, label="y = x", zorder=2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ssDNA clonal frequency"); ax.set_ylabel("dsDNA clonal frequency")
    ax.set_title(f"R{r}")
    biased = [s for s in mapped if ds[s][r] > 0 and ss[s][r] / ds[s][r] > 2]
    npos = sum(y[s] for s in biased)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fn = OUT / f"fig3A_scatter_elisa_R{r}_sized.png"
    fig.savefig(fn); plt.close(fig)
    # report
    lr = [np.log10(ss[s][r] / ds[s][r]) for s in mapped if ds[s][r] > 0]
    lp = [np.log10(ss[s][r] / ds[s][r]) for s in mapped if y[s] == 1 and ds[s][r] > 0]
    ln = [np.log10(ss[s][r] / ds[s][r]) for s in mapped if y[s] == 0 and ds[s][r] > 0]
    print(f"R{r}: ss-biased ELISA+ {npos}/{len(biased)} | median logratio + {np.median(lp):+.2f} - {np.median(ln):+.2f} | MWU p={mannwhitneyu(lp, ln).pvalue:.1e}")
    return fn


def render_table():
    def ratio(s, r):
        return ss[s][r] / ds[s][r] if ds[s][r] else 99

    def rowdata(s, tag):
        return [s, len(s), f"{ratio(s,3):.1f}", f"{ratio(s,4):.1f}",
                "+" if y[s] == 1 else "−", tag]
    ssb = sorted([s for s in mapped if y[s] == 1 and ratio(s, 3) > 2],
                 key=lambda s: -ratio(s, 3))[:5]
    val = [s for s in ["SADSCATCATYPSEIDA"] if s in mapped]
    cn = sorted([s for s in mapped if y[s] == 0 and 0.5 <= ratio(s, 3) <= 2],
                key=lambda s: -(ss[s][3] + ds[s][3]))[:3]
    rows = ([rowdata(s, "ssDNA-biased binder") for s in ssb]
            + [rowdata(s, "concordant binder") for s in val]
            + [rowdata(s, "concordant non-binder") for s in cn])
    cols = ["CDRH3", "len", "ss/ds R3", "ss/ds R4", "ELISA", "class"]
    fig, ax = plt.subplots(figsize=(7.6, 0.42 * (len(rows) + 1) + 0.4))
    ax.axis("off")
    t = ax.table(cellText=rows, colLabels=cols, loc="center", cellLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(11); t.scale(1, 1.55)
    fill = {"ssDNA-biased binder": "#ece1f3", "concordant binder": "#e9eef5",
            "concordant non-binder": "#f2f2f2"}
    for (ri, ci), cell in t.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if ri == 0:
            cell.set_facecolor("#1d3557"); cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor(fill[rows[ri - 1][-1]])
            if ci == 0:
                cell.get_text().set_fontfamily("monospace"); cell.get_text().set_fontsize(9.5)
    t.auto_set_column_width(list(range(len(cols))))
    fn = OUT / "fig3C_clone_table_sized.png"
    fig.savefig(fn, bbox_inches="tight"); plt.close(fig)
    print("wrote", fn.name, "| clones:", [r[0] for r in rows])
    return fn


for rr in (0, 1, 2, 3, 4):
    render_scatter(rr)
render_table()
print("done")
