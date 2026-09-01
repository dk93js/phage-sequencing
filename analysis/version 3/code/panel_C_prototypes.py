"""
Prototype candidates for Figure 1 Panel C.

Builds two options from data/clusters_pooled.csv (PPM per round, pooled ss+ds):
  (1) clonal_takeover_stack_pooled.png  - 100% stacked clonal-composition bar,
      top clones shown individually and colored by cluster, rest = grey "other".
  (2) rank_abundance_by_round_pooled.png - rank-abundance (Zipf) curves per round.

Standalone / non-destructive: does not touch make_panels.py outputs.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parent
PANELS = ROOT / "panels"
DATA = ROOT / "data"

C_CLUSTER = {1: "#E63946", 2: "#457B9D", 3: "#2A9D8F"}
CLUSTER_LABEL = {1: "C1 enriching", 2: "C2 intermediate", 3: "C3 depleted"}
ROUNDS = ["R0", "R1", "R2", "R3", "R4"]

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 13,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def lighten(color, f):
    """Blend color toward white by fraction f in [0,1]."""
    r, g, b = to_rgb(color)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


df = pd.read_csv(DATA / "clusters_pooled.csv")
ppm = df[[f"PPM{i}" for i in range(5)]].to_numpy(dtype=float)   # (n_clones, 5)
# per-round fraction of total reads (equal ss/ds weighting; pooled like Panel B)
frac = ppm / ppm.sum(axis=0, keepdims=True)                      # columns sum to 1

# ----------------------------------------------------------------------
# (1) clonal-takeover stacked bar
# ----------------------------------------------------------------------
TOP = 20
peak = frac.max(axis=1)
top_idx = np.argsort(-peak)[:TOP]
# order tracked clones: grouped by cluster (1,2,3), then by R4 fraction desc
clu = df["cluster"].to_numpy()
top_sorted = sorted(top_idx, key=lambda i: (clu[i], -frac[i, 4]))

# within-cluster lightness ramp so adjacent slices are distinguishable
shade = {}
for c in (1, 2, 3):
    members = [i for i in top_sorted if clu[i] == c]
    for rank, i in enumerate(members):
        f = 0.15 + 0.45 * (rank / max(1, len(members) - 1)) if len(members) > 1 else 0.2
        shade[i] = lighten(C_CLUSTER[c], f)

fig, ax = plt.subplots(figsize=(4.6, 4.4))
x = np.arange(5)
bottom = np.zeros(5)
for i in top_sorted:
    ax.bar(x, frac[i], bottom=bottom, width=0.78, color=shade[i],
           edgecolor="white", linewidth=0.4)
    bottom += frac[i]
other = 1 - bottom
ax.bar(x, other, bottom=bottom, width=0.78, color="#dcdcdc",
       edgecolor="white", linewidth=0.4, label="other")

ax.set_xticks(x)
ax.set_xticklabels(ROUNDS)
ax.set_xlabel("Selection round")
ax.set_ylabel("Fraction of reads")
ax.set_ylim(0, 1)
ax.set_title("Clonal composition (top 20 clones)")
# legend by cluster family
from matplotlib.patches import Patch
handles = [Patch(facecolor=C_CLUSTER[c], label=CLUSTER_LABEL[c]) for c in (1, 2, 3)]
handles.append(Patch(facecolor="#dcdcdc", label="other"))
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.0),
          frameon=False, fontsize=10)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(PANELS / "clonal_takeover_stack_pooled.png")
plt.close(fig)

# report the top-clone share per round
share = bottom.copy()
print("Top-20 clone cumulative share of reads by round:")
for r, s in zip(ROUNDS, share):
    print(f"  {r}: {s*100:5.1f}%   (other {100-s*100:5.1f}%)")

# ----------------------------------------------------------------------
# (2) rank-abundance (Zipf) curves
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4.6, 4.4))
cmap = plt.get_cmap("viridis")
for j in range(5):
    vals = np.sort(frac[:, j])[::-1]
    vals = vals[vals > 0]
    rank = np.arange(1, len(vals) + 1)
    ax.plot(rank, vals * 100, color=cmap(j / 4), lw=1.8, label=ROUNDS[j])
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Clone rank")
ax.set_ylabel("Frequency (% of reads)")
ax.set_title("Rank-abundance by round")
ax.legend(title="round", frameon=False, fontsize=10)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(PANELS / "rank_abundance_by_round_pooled.png")
plt.close(fig)

print("\nSaved:")
print("  panels/clonal_takeover_stack_pooled.png")
print("  panels/rank_abundance_by_round_pooled.png")
