"""
Resolve the 788-vs-1003 shared-clone gap (figure_checklist_2026-06-01 blocking item #1).

Hypothesis (from make_panels.compute_per_cluster_stats): df["cluster"] is the pooled
Ward label assigned PER ROW. A shared clonotype (present in both ss and ds) has two
rows and therefore two independent cluster labels. The "all" line counts a shared
clone once; a per-cluster line only retains it if BOTH copies carry that cluster.
So shared clones whose ss-copy and ds-copy disagree on cluster are in "all" but in
no colored line. This script confirms that and quantifies the discordance.

Reads the already-computed data/clusters_pooled.csv (no re-clustering).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

CODE_DIR = Path(__file__).resolve().parent
DATA = CODE_DIR.parent / "data"
PANELS = CODE_DIR.parent / "panels"
ROUNDS = [0, 1, 2, 3, 4]
CL_LABEL = {1: "C1 enriching", 2: "C2 intermediate", 3: "C3 depleted"}
C_CLUSTER = {1: "#E63946", 2: "#457B9D", 3: "#2A9D8F"}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 16,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

pooled = pd.read_csv(DATA / "clusters_pooled.csv")   # seq, lib, cluster, PPM0..4, z0..4

# pivot cluster label per seq per lib
cl = pooled.pivot_table(index="seq", columns="lib", values="cluster", aggfunc="first")
ppm = {lib: pooled[pooled.lib == lib].set_index("seq") for lib in ["ss", "ds"]}

shared_seqs = cl.dropna(subset=["ss", "ds"]).index          # in both libraries at all
cl_shared = cl.loc[shared_seqs].astype(int)
n_shared_total = len(cl_shared)
concordant = (cl_shared["ss"] == cl_shared["ds"])
print(f"Sequences in BOTH libraries (trajectory-level): {n_shared_total}")
print(f"  concordant (same pooled cluster ss & ds): {concordant.sum()} "
      f"({100*concordant.mean():.1f}%)")
print(f"  DISCORDANT (ss-copy and ds-copy differ):   {(~concordant).sum()} "
      f"({100*(~concordant).mean():.1f}%)")

print("\nConfusion matrix  rows = ssDNA-copy cluster, cols = dsDNA-copy cluster")
conf = pd.crosstab(cl_shared["ss"], cl_shared["ds"])
conf.index = [CL_LABEL[i] for i in conf.index]
conf.columns = [CL_LABEL[i] for i in conf.columns]
print(conf.to_string())
print(f"  off-diagonal total = {conf.values.sum() - np.trace(conf.values)}")

# Now reproduce per-round n_shared exactly as make_panels does, and the by-cluster sum
print("\nPer-round reconciliation (matches per_cluster_stats.csv):")
print(f"{'round':>5} {'all(PPM>0 both)':>16} {'C1':>5} {'C2':>5} {'C3':>5} {'C1+C2+C3':>9} {'gap':>5}")
recon = []
for r in ROUNDS:
    # "all": seq present (PPM_r>0) in BOTH libs
    ss_pos = set(ppm["ss"].index[ppm["ss"][f"PPM{r}"] > 0])
    ds_pos = set(ppm["ds"].index[ppm["ds"][f"PPM{r}"] > 0])
    both = ss_pos & ds_pos
    n_all = len(both)
    # by cluster: of `both`, keep where BOTH copies are in cluster k
    percl = {}
    for k in (1, 2, 3):
        keep = [s for s in both
                if cl.loc[s, "ss"] == k and cl.loc[s, "ds"] == k]
        percl[k] = len(keep)
    s = sum(percl.values())
    recon.append((r, n_all, percl[1], percl[2], percl[3], s, n_all - s))
    print(f"R{r:<4} {n_all:>16} {percl[1]:>5} {percl[2]:>5} {percl[3]:>5} {s:>9} {n_all-s:>5}")

# Break the gap down by which kind of discordance, at R3 (the headline round)
r = 3
ss_pos = set(ppm["ss"].index[ppm["ss"][f"PPM{r}"] > 0])
ds_pos = set(ppm["ds"].index[ppm["ds"][f"PPM{r}"] > 0])
both_r = sorted(ss_pos & ds_pos)
disc = [(s, int(cl.loc[s, "ss"]), int(cl.loc[s, "ds"]))
        for s in both_r if cl.loc[s, "ss"] != cl.loc[s, "ds"]]
print(f"\nAt R3: {len(disc)} discordant shared clones. Breakdown (ss→ds):")
from collections import Counter
for (a, b), n in sorted(Counter((a, b) for _, a, b in disc).items()):
    print(f"  {CL_LABEL[a]} (ss)  ->  {CL_LABEL[b]} (ds):  {n}")

# Of the discordant clones, how many involve C1 (the enriching cluster) on one side?
involves_c1 = [d for d in disc if 1 in (d[1], d[2])]
print(f"\n  discordant clones involving C1 on either side: {len(involves_c1)} / {len(disc)}")

out = pd.DataFrame(disc, columns=["seq", "ss_cluster", "ds_cluster"])
out.to_csv(DATA / "shared_clone_cluster_discordance_R3.csv", index=False)
print(f"\nwrote {DATA / 'shared_clone_cluster_discordance_R3.csv'}")


# =========================================================
# Panel: ss-vs-ds cluster confusion matrix for shared clonotypes
# =========================================================
def render_confusion_panel(conf_counts):
    """conf_counts: 3x3 np array, rows=ss cluster (1..3), cols=ds cluster (1..3)."""
    M = conf_counts.astype(int)
    total = M.sum()
    conc = np.trace(M)
    disc = total - conc

    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    # color by row-fraction so the dominant C3 diagonal does not wash out the rest
    row_frac = M / M.sum(axis=1, keepdims=True)
    im = ax.imshow(row_frac, cmap="Blues", vmin=0, vmax=1)

    for i in range(3):
        for j in range(3):
            on_diag = i == j
            ax.text(j, i, f"{M[i, j]}",
                    ha="center", va="center",
                    fontsize=15 if on_diag else 13,
                    fontweight="bold" if on_diag else "normal",
                    color="white" if row_frac[i, j] > 0.55 else "#222")
            if on_diag:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor="#E63946", linewidth=2.4))

    ax.set_xticks(range(3)); ax.set_yticks(range(3))
    ax.set_xticklabels([CL_LABEL[k] for k in (1, 2, 3)], rotation=20, ha="right")
    ax.set_yticklabels([CL_LABEL[k] for k in (1, 2, 3)])
    # color the tick labels by cluster
    for t, k in zip(ax.get_xticklabels(), (1, 2, 3)):
        t.set_color(C_CLUSTER[k])
    for t, k in zip(ax.get_yticklabels(), (1, 2, 3)):
        t.set_color(C_CLUSTER[k])
    ax.set_xlabel("dsDNA-copy cluster", fontweight="bold")
    ax.set_ylabel("ssDNA-copy cluster", fontweight="bold")
    ax.set_title(f"Shared clonotypes: ss vs ds cluster\n"
                 f"concordant {conc} ({100*conc/total:.0f}%), "
                 f"discordant {disc} ({100*disc/total:.0f}%)",
                 fontsize=13)
    ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    fig.savefig(PANELS / "shared_clone_cluster_confusion.png")
    plt.close(fig)
    print(f"wrote {PANELS / 'shared_clone_cluster_confusion.png'}")


# build 3x3 from trajectory-level concordance (matches confusion matrix printed above)
conf_counts = np.zeros((3, 3))
for _, a, b in zip(cl_shared.index, cl_shared["ss"], cl_shared["ds"]):
    conf_counts[a - 1, b - 1] += 1
render_confusion_panel(conf_counts)
