"""Main Fig 1B — pooled / dsDNA / ssDNA cluster trajectories (option B palette).

Supp Fig 1C is promoted to main Fig 1B, so this renders all three trajectory
panels in one consistent style:

  * descriptive cluster labels on every panel (the dsDNA panel used to carry
    terse "C1"/"C3" labels while ssDNA carried descriptive ones)
  * NO "(n=...)" in the panel legend — the counts belong in the figure caption
  * legend outside upper-right, as in the current panels

Run:  python "figures/panels/new color/optionB/code/render_recolored.py"
"""
import shutil
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_panels as mp
import make_panels_sized as mps

OUT = HERE.parent
mpl.rcParams.update({"axes.titlesize": 14, "legend.fontsize": 10})

LABEL = {1: "C1 enriching", 2: "C2 intermediate: peaks at R2", 3: "C3 depleted"}
OUTNAME = {"pooled": "fig1B_trajectory_pooled.png",
           "dsDNA": "fig1B_trajectory_dsDNA.png",
           "ssDNA": "fig1B_trajectory_ssDNA.png"}


def render(ds_name, info):
    z, labels = info["z"], info["labels"]
    fig, ax = plt.subplots(figsize=mps.B_FIGSIZE)
    counts = {}
    for k in (1, 2, 3):
        sel = labels == k
        counts[k] = int(sel.sum())
        if not sel.any():
            continue
        m, s = z[sel].mean(axis=0), z[sel].std(axis=0)
        ax.plot(mp.ROUNDS, m, "o-", color=mps.C_CLUSTER[k], label=LABEL[k],
                linewidth=2.4, markersize=8)
        ax.fill_between(mp.ROUNDS, m - s, m + s, color=mps.C_CLUSTER[k],
                        alpha=0.18, linewidth=0)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle=":")
    ax.set_xticks(mp.ROUNDS)
    ax.set_xticklabels([f"R{r}" for r in mp.ROUNDS])
    ax.set_xlabel("selection round")
    ax.set_ylabel("depth-normalized (z)")
    ax.set_title(mp.DATASET_LABEL[ds_name])
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.savefig(mp.PANELS / OUTNAME[ds_name])
    plt.close(fig)
    shutil.copy2(mp.PANELS / OUTNAME[ds_name], OUT / OUTNAME[ds_name])
    return counts


df, info = mp.load_and_cluster()

print("\ncluster sizes for the figure caption (n is no longer on the panels):")
for ds in ("pooled", "dsDNA", "ssDNA"):
    c = render(ds, info[ds])
    print(f"  {mp.DATASET_LABEL[ds]:22} C1 {c[1]:5}   C2 {c[2]:5}   C3 {c[3]:5}   total {sum(c.values()):5}")
print(f"\nwrote 3 panels to {OUT}")
