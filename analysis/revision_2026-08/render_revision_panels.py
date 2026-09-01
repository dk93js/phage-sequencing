"""Panels for the two new supplementary figures of the Biomolecules revision.

Supp Fig 1C-E (R2 point 2)  - choice of k and cluster stability
Supp Fig 3E-F (R2 point 3)  - the 215 cluster-discordant clonotypes

Style follows the published panels: Arial 13 / title 14 / label 14 / tick 12 /
legend 10, 3.6-in panel height, 300 dpi, no top-right spines, cluster palette
C1 #E63946, C2 #17375E, C3 #56C1B0.
"""
from pathlib import Path
import os
os.makedirs("out/figures", exist_ok=True)
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "lines.linewidth": 2.0, "lines.markersize": 6,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})
H = 3.6
OUT = Path("out/figures")
R = Path("out")
C1, C2, C3 = "#E63946", "#17375E", "#56C1B0"
GREY = "#8A8A8A"

def save(fig, name):
    fig.savefig(OUT / f"{name}.png"); plt.close(fig)
    print("wrote", OUT / f"{name}.png")

# ---------------- Supp Fig 6 -------------------------------------------------
idx  = pd.read_csv(R / "r2_02_k_indices.csv")
stab = pd.read_csv(R / "r2_02_bootstrap_stability.csv")
conc = pd.read_csv(R / "r2_02_conclusion_vs_k.csv")

fig, ax = plt.subplots(figsize=(H * 1.15, H))
ax.plot(idx.k, idx.silhouette, "o-", color=C1, label="silhouette")
ax.set_xlabel("number of clusters (k)"); ax.set_ylabel("silhouette", color=C1)
ax.tick_params(axis="y", colors=C1); ax.set_ylim(0, 0.65)
ax2 = ax.twinx(); ax2.spines["top"].set_visible(False)
ax2.plot(idx.k, idx.gap, "s--", color=C2, label="gap statistic")
ax2.fill_between(idx.k, idx.gap - idx.gap_se, idx.gap + idx.gap_se, color=C2, alpha=0.15, lw=0)
ax2.set_ylabel("gap statistic", color=C2); ax2.tick_params(axis="y", colors=C2)
ax.axvline(3, color=GREY, ls=":", lw=1.5)
ax.set_xticks(idx.k)
ax.text(3.1, 0.03, "k = 3\n(used)", fontsize=10, color="#555555")
ax.set_title("Internal validity indices vs. k")
save(fig, "suppfig1C_k_indices")

fig, ax = plt.subplots(figsize=(H * 1.05, H))
w = 0.38
ax.bar(stab.k - w/2, stab.mean_jaccard, w, color=C2, label="mean over clusters")
ax.bar(stab.k + w/2, stab.min_cluster,  w, color=C3, label="least stable cluster")
ax.axhline(0.75, color=C1, ls="--", lw=1.5)
ax.text(6.35, 0.765, "stable", color=C1, fontsize=11, ha="right")
ax.set_xlabel("number of clusters (k)"); ax.set_ylabel("bootstrap Jaccard")
ax.set_ylim(0, 1.28); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_xticks(stab.k)
ax.legend(frameon=False, loc="upper right", ncol=1)
ax.set_title("Cluster stability\n(200 bootstrap resamples)")
save(fig, "suppfig1D_bootstrap_jaccard")

fig, ax = plt.subplots(figsize=(H * 1.05, H))
for rnd, col, mk in [("R2", C3, "^"), ("R3", C1, "o"), ("R4", C2, "s")]:
    s = conc[conc["round"] == rnd]
    ax.plot(s.k, s.slope, mk + "-", color=col, label=rnd)
ax.axhline(1.0, color=GREY, ls=":", lw=1.5)
for k, g in conc[conc["round"] == "R3"].groupby("k"):
    ax.text(k, 1.10, f"n={int(g.n.iloc[0])}", ha="center", fontsize=10, color="#555555")
ax.set_xlabel("number of clusters (k)"); ax.set_ylabel("log–log slope (dsDNA on ssDNA)")
ax.set_ylim(0.3, 1.18); ax.set_xlim(1.6, 6.4); ax.set_xticks(sorted(conc.k.unique()))
ax.legend(frameon=False, title=None, loc="center left", bbox_to_anchor=(0.02, 0.45))
ax.set_title("Enriching-cluster slope across k")
save(fig, "suppfig1E_slope_vs_k")

# ---------------- Supp Fig 7 -------------------------------------------------
d = pd.read_csv(R / "r2_03_shared_clonotype_table.csv", index_col=0)
order = [(1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3)]
keep  = [t for t in order if ((d.ds_lab==t[0]) & (d.ss_lab==t[1])).sum() >= 5]

fig, ax = plt.subplots(figsize=(H * 1.7, H))
data, labs, cols = [], [], []
for ds_c, ss_c in keep:
    g = d[(d.ds_lab==ds_c) & (d.ss_lab==ss_c)]
    data.append(g.lr3.dropna().values)
    labs.append(f"C{ds_c}→C{ss_c}\n(n={len(g)})")
    cols.append({1:C1,2:C2,3:C3}[ss_c] if ds_c != ss_c else "#BFBFBF")
bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6)
for p, c in zip(bp["boxes"], cols):
    p.set_facecolor(c); p.set_alpha(0.65); p.set_edgecolor("#333333")
for e in ("whiskers","caps","medians"):
    for l in bp[e]: l.set_color("#333333")
ax.axhline(0, color=GREY, ls=":", lw=1.5)
ax.set_xticklabels(labs, fontsize=10)
ax.set_xlabel("dsDNA-copy cluster → ssDNA-copy cluster")
ax.set_ylabel("log$_{10}$(ssDNA / dsDNA) at R3")
ax.set_title("Discordant clonotypes by transition class")
save(fig, "suppfig3E_transition_logratio")

a = d[d.elisa.notna()]
ds1 = a[a.ds_lab == 1]
grp = [("ssDNA agrees\n(C1 → C1)", ds1[ds1.ss_lab == 1]),
       ("ssDNA disagrees\n(C1 → C2/C3)", ds1[ds1.ss_lab != 1])]
fig, ax = plt.subplots(figsize=(H * 0.95, H))
vals = [(g.elisa == 1).mean() * 100 for _, g in grp]
ns   = [f"{int((g.elisa==1).sum())}/{len(g)}" for _, g in grp]
ax.bar([0, 1], vals, 0.55, color=[C1, GREY])
for i, (v, t) in enumerate(zip(vals, ns)):
    ax.text(i, v + 2.5, t, ha="center", fontsize=12)
tab = [[int((g.elisa==1).sum()), int((g.elisa==0).sum())] for _, g in grp]
p = stats.fisher_exact(tab)[1]
ax.plot([0, 0, 1, 1], [88, 91, 91, 88], color="#333333", lw=1.2)
exp = int(np.floor(np.log10(p))); mant = p / 10 ** exp
ax.text(0.5, 92.5, f"P = {mant:.1f} × 10$^{{{exp}}}$", ha="center", fontsize=12)
ax.set_xticks([0, 1]); ax.set_xticklabels([g[0] for g in grp], fontsize=11)
ax.set_ylabel("antigen-reactive (%)"); ax.set_ylim(0, 105)
ax.set_title("dsDNA-called enriching clonotypes")
save(fig, "suppfig3F_elisa_concordance")
