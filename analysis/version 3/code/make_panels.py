"""
Atomic panels for the phage-seq paper (version 3).

All panels are written to ../panels/ with descriptive filenames so they can be
assembled freely into figures or supplements. Each panel renders as a
standalone, self-titled PNG at 300 DPI.

Panels produced
---------------
Methodology / data overview:
    read_depth_per_round_table.png
    enrichment_heatmap_pooled.png
    enrichment_heatmap_ssDNA.png
    enrichment_heatmap_dsDNA.png
    cluster_trajectory_pooled.png
    cluster_trajectory_ssDNA.png
    cluster_trajectory_dsDNA.png

ssDNA-vs-dsDNA bias finding:
    spearman_rho_by_cluster.png
    loglog_slope_by_cluster.png
    scatter_C1_enriching_R3.png       (3 clusters × 2 rounds = 6 atomic scatters)
    scatter_C1_enriching_R4.png
    scatter_C2_intermediate_R3.png
    scatter_C2_intermediate_R4.png
    scatter_C3_depleted_R3.png
    scatter_C3_depleted_R4.png
    clone_trajectory_GASGYSSGNIDA_vs_SADSCATCATYPSEIDA.png

Derived data
------------
data/clusters_pooled.csv, clusters_ssDNA.csv, clusters_dsDNA.csv
data/read_depth_per_round.csv
data/per_cluster_stats.csv         (Spearman ρ and log-log slope per cluster per round)
data/cluster_size_summary.csv      (n per cluster per dataset)

Methods
-------
- Clustering: Ward hierarchical on row-z-scored PPM trajectories, k = 3.
- Canonical labels: cluster 1 = highest R4–R0 mean z (enriching),
  cluster 3 = lowest (depleted), cluster 2 = remaining.
- Per-cluster ρ and slope (Panels: spearman/slope/scatter): use the pooled
  Ward clustering; for each cluster and round, restrict to clonotypes
  detected (PPM > 0) in both libraries.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import spearmanr, linregress

CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parent
PANELS = ROOT / "panels"
DATA = ROOT / "data"
PANELS.mkdir(exist_ok=True)
DATA.mkdir(exist_ok=True)

SRC = Path(__file__).resolve().parents[3] / "data" / "report_q40_ppm100_list.csv"
assert SRC.exists(), f"source CSV not found: {SRC}"

ROUNDS = [0, 1, 2, 3, 4]
PPM_COLS = [f"PPM{r}" for r in ROUNDS]
P_COLS = [f"P{r}" for r in ROUNDS]
LIBS = ["ss", "ds"]
LIB_LABEL = {"ss": "ssDNA", "ds": "dsDNA"}

C_CLUSTER = {1: "#E63946", 2: "#457B9D", 3: "#2A9D8F"}
CLUSTER_LABEL = {1: "C1 enriching", 2: "C2 intermediate", 3: "C3 depleted"}
CLUSTER_SHORT_NAME = {1: "C1_enriching", 2: "C2_intermediate", 3: "C3_depleted"}
DATASETS = ["pooled", "ssDNA", "dsDNA"]
DATASET_LABEL = {"pooled": "all (ssDNA + dsDNA)", "ssDNA": "ssDNA", "dsDNA": "dsDNA"}

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 13,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.1,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

CMAP_DIV = LinearSegmentedColormap.from_list("rwb", ["#2166AC", "#F7F7F7", "#B2182B"])


# =========================================================
# data + clustering
# =========================================================
def zscore_rows(mat: np.ndarray) -> np.ndarray:
    mu = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, keepdims=True)
    sd[sd == 0] = 1
    return (mat - mu) / sd


def _lighten(color, f):
    """Blend `color` toward white by fraction f in [0, 1]."""
    r, g, b = to_rgb(color)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def canonicalize(labels: np.ndarray, z: np.ndarray) -> np.ndarray:
    uniq = [int(u) for u in np.unique(labels)]
    slopes = {u: float(z[labels == u, 4].mean() - z[labels == u, 0].mean()) for u in uniq}
    order = sorted(uniq, key=lambda u: slopes[u], reverse=True)
    mapping = {old: new for new, old in zip([1, 2, 3], order)}
    out = np.full_like(labels, -1, dtype=int)
    for old, new in mapping.items():
        out[labels == old] = new
    return out


def cluster_ward(z: np.ndarray, k: int = 3) -> np.ndarray:
    Z = linkage(z, method="ward")
    labels = fcluster(Z, t=k, criterion="maxclust").astype(int)
    return canonicalize(labels, z)


def load_and_cluster() -> tuple:
    """Return df (with pooled-Ward cluster column) + per-dataset cluster info."""
    df = pd.read_csv(SRC)
    info = {}
    for ds_name in DATASETS:
        if ds_name == "pooled":
            sub = df.reset_index(drop=True)
        elif ds_name == "ssDNA":
            sub = df[df["lib"] == "ss"].reset_index(drop=True)
        else:
            sub = df[df["lib"] == "ds"].reset_index(drop=True)
        z = zscore_rows(sub[PPM_COLS].to_numpy(dtype=float))
        labels = cluster_ward(z)
        slope = z[:, 4] - z[:, 0]
        order = np.lexsort((-slope, labels))
        # save per-dataset cluster table
        out = sub[["seq", "lib"]].copy()
        out["cluster"] = labels
        for r in ROUNDS:
            out[f"PPM{r}"] = sub[f"PPM{r}"]
            out[f"z{r}"] = z[:, r]
        out.to_csv(DATA / f"clusters_{ds_name}.csv", index=False)
        info[ds_name] = {"df": sub, "z": z, "labels": labels, "order": order}

    # attach pooled cluster onto the master df for downstream bias panels
    pooled = pd.read_csv(DATA / "clusters_pooled.csv")[["seq", "lib", "cluster"]]
    df = df.merge(pooled, on=["seq", "lib"], how="left", suffixes=("_orig", ""))
    if "cluster_orig" in df.columns:
        df = df.drop(columns=["cluster_orig"])
    df["cluster"] = df["cluster"].astype(int)

    # cluster-size summary
    rows = []
    for ds_name in DATASETS:
        labels = info[ds_name]["labels"]
        for k in [1, 2, 3]:
            rows.append({"dataset": ds_name, "cluster": k,
                         "label": CLUSTER_LABEL[k],
                         "n": int((labels == k).sum())})
    pd.DataFrame(rows).to_csv(DATA / "cluster_size_summary.csv", index=False)
    return df, info


# =========================================================
# Panel: read-depth table
# =========================================================
def make_read_depth_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lib in LIBS:
        sub = df[df["lib"] == lib]
        for r in ROUNDS:
            p_sum = int(sub[f"P{r}"].sum())
            ppm_sum = float(sub[f"PPM{r}"].sum())
            nz = sub[(sub[f"P{r}"] > 0) & (sub[f"PPM{r}"] > 0)]
            implied_total = int((nz[f"P{r}"] / nz[f"PPM{r}"]).median() * 1e6) if len(nz) else 0
            rows.append({
                "library": LIB_LABEL[lib],
                "round": f"R{r}",
                "raw_reads_in_filtered_set": p_sum,
                "depth_normalized_filtered_set": ppm_sum,
                "Q40_filtered_total_reads": implied_total,
            })
    tbl = pd.DataFrame(rows)
    tbl.to_csv(DATA / "read_depth_per_round.csv", index=False)
    return tbl


def render_read_depth_table(df: pd.DataFrame, tbl: pd.DataFrame):
    n_clones = {LIB_LABEL[lib]: int((df["lib"] == lib).sum()) for lib in LIBS}

    fig, ax = plt.subplots(figsize=(10.5, 3.0))
    ax.axis("off")

    header = ["Library", "Metric"] + [f"R{r}" for r in ROUNDS] + ["Clonotypes"]
    rows_disp = []
    for lib_name in ["ssDNA", "dsDNA"]:
        sub = tbl[tbl["library"] == lib_name].set_index("round")
        rows_disp.append(
            [lib_name, "Raw reads"]
            + [f"{int(sub.loc[f'R{r}', 'raw_reads_in_filtered_set']):,}" for r in ROUNDS]
            + [f"{n_clones[lib_name]:,}"]
        )
        rows_disp.append(
            ["", "Depth-normalized"]
            + [f"{int(sub.loc[f'R{r}', 'depth_normalized_filtered_set']):,}" for r in ROUNDS]
            + ["—"]
        )
    table = ax.table(
        cellText=rows_disp, colLabels=header, loc="center", cellLoc="center", colLoc="center",
        colWidths=[0.10, 0.20] + [0.11] * len(ROUNDS) + [0.12],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.0, 1.7)
    for (r, c), cell in table.get_celld().items():
        cell.set_linewidth(0.6)
        if r == 0:
            cell.set_text_props(weight="bold", color="white", fontsize=10)
            cell.set_facecolor("#2a3a4a")
        else:
            cell.set_facecolor("#fdebec" if r in (1, 2) else "#e8eef5")
            if c == 0:
                cell.set_text_props(weight="bold")
            if c == 1:
                cell.set_text_props(style="italic")
            if r in (2, 4):
                cell.set_text_props(style="italic")
    ax.text(
        0.0, -0.05,
        r"Clonotype = unique CDRH3 amino-acid sequence per library, passing "
        r"q$\geq$40 (Phred over CDRH3) and depth-normalized value $\geq 100$ in at least one round."
        "\n"
        r"Depth-normalized = (clonotype reads at round $r$) / (Q40-filtered total reads at round $r$) $\times\,10^6$, "
        r"summed over the analysis set.",
        transform=ax.transAxes, va="top", ha="left", fontsize=10, color="#444",
    )
    fig.savefig(PANELS / "read_depth_per_round_table.png")
    plt.close(fig)
    print("wrote panels/read_depth_per_round_table.png")


# =========================================================
# Panel: enrichment heatmap (one per dataset)
# =========================================================
def render_enrichment_heatmap(ds_name: str, info: dict):
    z_ord = info["z"][info["order"]]
    cl_ord = info["labels"][info["order"]]
    n = len(cl_ord)

    fig = plt.figure(figsize=(4.2, 5.4))
    gs = gridspec.GridSpec(1, 3, width_ratios=[0.10, 1.0, 0.06], wspace=0.07)
    ax_anno = fig.add_subplot(gs[0])
    ax_hm = fig.add_subplot(gs[1], sharey=ax_anno)
    ax_cb = fig.add_subplot(gs[2])

    cmap_anno = mpl.colors.ListedColormap([C_CLUSTER[1], C_CLUSTER[2], C_CLUSTER[3]])
    ax_anno.imshow(np.array([cl_ord]).T, aspect="auto", cmap=cmap_anno,
                   interpolation="nearest", extent=[0, 1, n, 0])
    ax_anno.set_xticks([])
    ax_anno.set_yticks([])
    ax_anno.set_xlim(0, 1)
    ax_anno.set_ylabel(f"n = {n}")
    # label each contiguous cluster block with its ID (C1/C2/C3)
    for c in (1, 2, 3):
        idx = np.where(cl_ord == c)[0]
        if len(idx):
            ax_anno.text(0.5, idx.mean() + 0.5, f"C{c}", ha="center", va="center",
                         color="white", fontweight="bold", fontsize=11, rotation=90)

    vmax = 2.5
    im = ax_hm.imshow(z_ord, aspect="auto", cmap=CMAP_DIV, vmin=-vmax, vmax=vmax,
                      interpolation="nearest", extent=[-0.5, 4.5, n, 0])
    ax_hm.set_xticks(ROUNDS)
    ax_hm.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax_hm.set_yticks([])
    ax_hm.set_xlabel("selection round")
    ax_hm.set_title(DATASET_LABEL[ds_name], loc="center", fontsize=14)

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("depth-normalized (z)")
    cb.outline.set_linewidth(0)
    fig.savefig(PANELS / f"enrichment_heatmap_{ds_name}.png")
    plt.close(fig)
    print(f"wrote panels/enrichment_heatmap_{ds_name}.png")


# =========================================================
# Panel: cluster trajectory (one per dataset)
# =========================================================
def render_cluster_trajectory(ds_name: str, info: dict, short_labels: bool = False):
    z = info["z"]
    labels = info["labels"]
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    for k in [1, 2, 3]:
        sel = labels == k
        if not sel.any():
            continue
        m = z[sel].mean(axis=0)
        s = z[sel].std(axis=0)
        label = f"C{k} (n={int(sel.sum())})" if short_labels else f"{CLUSTER_LABEL[k]} (n={int(sel.sum())})"
        ax.plot(ROUNDS, m, "o-", color=C_CLUSTER[k], label=label,
                linewidth=2.4, markersize=8)
        ax.fill_between(ROUNDS, m - s, m + s, color=C_CLUSTER[k], alpha=0.18, linewidth=0)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle=":")
    ax.set_xticks(ROUNDS)
    ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("Selection round")
    ax.set_ylabel("Normalized PPM")
    ax.set_title(DATASET_LABEL[ds_name], loc="center", fontweight="bold", fontsize=14)
    ax.legend(frameon=False, loc="upper left", fontsize=10)
    fig.savefig(PANELS / f"cluster_trajectory_{ds_name}.png")
    plt.close(fig)
    print(f"wrote panels/cluster_trajectory_{ds_name}.png")


# =========================================================
# Panel: clonal-takeover stacked composition (Figure 1C)
# =========================================================
def render_clonal_takeover_stack(df: pd.DataFrame, top_n: int = 20,
                                 annotate_top: bool = True):
    """100% stacked clonal composition per round (pooled ss+ds).

    The `top_n` clonotypes with the highest peak frequency are drawn as
    individual slices, colored by their pooled Ward cluster with a
    within-cluster lightness ramp so adjacent clones stay distinguishable;
    all remaining clonotypes are collapsed into a grey 'other' slice.
    Each round column is normalized to the round's total reads, so the
    panel reads as the fraction of the sequenced pool captured by the
    enriching (C1), intermediate (C2) and depleted (C3) clusters.
    """
    ppm = df[PPM_COLS].to_numpy(dtype=float)
    frac = ppm / ppm.sum(axis=0, keepdims=True)          # each round column sums to 1
    clu = df["cluster"].to_numpy()
    seqs = df["seq"].to_numpy()

    peak = frac.max(axis=1)
    top_idx = np.argsort(-peak)[:top_n]
    # stacking order: grouped by cluster (1,2,3), largest R4 slice at the base
    top_sorted = sorted(top_idx, key=lambda i: (clu[i], -frac[i, 4]))

    shade = {}
    for c in (1, 2, 3):
        members = [i for i in top_sorted if clu[i] == c]
        for rank, i in enumerate(members):
            f = 0.15 + 0.45 * (rank / max(1, len(members) - 1)) if len(members) > 1 else 0.2
            shade[i] = _lighten(C_CLUSTER[c], f)

    fig, ax = plt.subplots(figsize=(4.8, 5.0))
    x = np.arange(5)
    bottom = np.zeros(5)
    for i in top_sorted:
        ax.bar(x, frac[i], bottom=bottom, width=0.80, color=shade[i],
               edgecolor="white", linewidth=0.4)
        bottom += frac[i]
    ax.bar(x, 1 - bottom, bottom=bottom, width=0.80, color="#dcdcdc",
           edgecolor="white", linewidth=0.4)

    ax.set_xticks(x)
    ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("Selection round")
    ax.set_ylabel("Fraction of reads")
    ax.set_ylim(0, 1)
    ax.set_title("Clonal composition", loc="center", fontweight="bold", fontsize=14)

    handles = [Patch(facecolor=C_CLUSTER[c], label=CLUSTER_LABEL[c]) for c in (1, 2, 3)]
    handles.append(Patch(facecolor="#dcdcdc", label=f"other (rank >{top_n})"))
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=9.5)

    if annotate_top:
        dom = max(top_sorted, key=lambda i: frac[i, 4])
        below = sum(frac[i, 4] for i in top_sorted[:top_sorted.index(dom)])
        yc = below + frac[dom, 4] / 2
        ax.annotate(f"{seqs[dom]}\n{frac[dom, 4] * 100:.1f}% at R4",
                    xy=(4.40, yc), xytext=(4.62, yc),
                    va="center", ha="left", fontsize=8, color="#333",
                    arrowprops=dict(arrowstyle="-", color="#888", lw=0.8),
                    annotation_clip=False)

    fig.savefig(PANELS / "clonal_takeover_stack_pooled.png")
    plt.close(fig)
    print("wrote panels/clonal_takeover_stack_pooled.png")


# =========================================================
# Bias panels (use pooled Ward clusters on the master df)
# =========================================================
def compute_per_cluster_stats(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cl in [None, 1, 2, 3]:
        for r in ROUNDS:
            if cl is None:
                sub = df.copy(); label = "all"
            else:
                sub = df[df["cluster"] == cl].copy()
                label = f"C{cl}"
            piv = sub.pivot_table(index="seq", columns="lib", values=f"PPM{r}",
                                  aggfunc="first").fillna(0)
            if "ss" not in piv or "ds" not in piv:
                continue
            mask = (piv["ss"] > 0) & (piv["ds"] > 0)
            ss = piv.loc[mask, "ss"].to_numpy()
            ds = piv.loc[mask, "ds"].to_numpy()
            if len(ss) < 5:
                continue
            rho, _ = spearmanr(ss, ds)
            lr = linregress(np.log10(ss), np.log10(ds))
            rows.append({"cluster": label, "round": r, "n_shared": int(len(ss)),
                         "spearman": rho, "slope_loglog": lr.slope,
                         "intercept_loglog": lr.intercept})
    out = pd.DataFrame(rows)
    out.to_csv(DATA / "per_cluster_stats.csv", index=False)
    return out


def render_spearman_by_cluster(stats: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    for cl_key, color, label in [
        ("all", "black", "All"),
        ("C1", C_CLUSTER[1], "C1 enriching"),
        ("C2", C_CLUSTER[2], "C2 intermediate"),
        ("C3", C_CLUSTER[3], "C3 depleted"),
    ]:
        sub = stats[stats["cluster"] == cl_key].sort_values("round")
        ax.plot(sub["round"], sub["spearman"], "o-", color=color, label=label,
                linewidth=2.4, markersize=8)
    ax.set_xticks(ROUNDS)
    ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("Selection round")
    ax.set_ylabel("Spearman ρ (ds vs ss)")
    ax.set_ylim(0.6, 1.02)
    ax.axhline(1.0, color="gray", linewidth=0.8, linestyle=":")
    ax.set_title("Rank correlation by cluster", loc="center", fontweight="bold")
    ax.legend(frameon=False, loc="lower left", fontsize=10)
    fig.savefig(PANELS / "spearman_rho_by_cluster.png")
    plt.close(fig)
    print("wrote panels/spearman_rho_by_cluster.png")


def render_loglog_slope_by_cluster(stats: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    for cl_key, color, label in [
        ("all", "black", "All"),
        ("C1", C_CLUSTER[1], "C1 enriching"),
        ("C2", C_CLUSTER[2], "C2 intermediate"),
        ("C3", C_CLUSTER[3], "C3 depleted"),
    ]:
        sub = stats[stats["cluster"] == cl_key].sort_values("round")
        ax.plot(sub["round"], sub["slope_loglog"], "o-", color=color, label=label,
                linewidth=2.4, markersize=8)
    ax.set_xticks(ROUNDS)
    ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("Selection round")
    ax.set_ylabel("Log-log slope (ds vs ss)")
    ax.set_ylim(0.45, 1.1)
    ax.axhline(1.0, color="gray", linewidth=0.8, linestyle=":")
    ax.set_title("Quantitative bias by cluster", loc="center", fontweight="bold")
    ax.legend(frameon=False, loc="lower left", fontsize=10)
    fig.savefig(PANELS / "loglog_slope_by_cluster.png")
    plt.close(fig)
    print("wrote panels/loglog_slope_by_cluster.png")


def render_scatter_cluster_round(df: pd.DataFrame, cl: int, r: int):
    sub = df[df["cluster"] == cl]
    piv = sub.pivot_table(index="seq", columns="lib", values=f"PPM{r}",
                          aggfunc="first").fillna(0)
    mask = (piv["ss"] > 0) & (piv["ds"] > 0)
    ss = piv.loc[mask, "ss"].to_numpy()
    ds = piv.loc[mask, "ds"].to_numpy()

    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    ax.scatter(ss, ds, s=14, color=C_CLUSTER[cl], alpha=0.6, edgecolor="none")
    lo = min(ss.min(), ds.min()) * 0.7
    hi = max(ss.max(), ds.max()) * 1.4
    ax.plot([lo, hi], [lo, hi], color="red", linewidth=1.2, label="y = x")
    lr = linregress(np.log10(ss), np.log10(ds))
    xs = np.array([lo, hi])
    ys = 10 ** (lr.intercept + lr.slope * np.log10(xs))
    ax.plot(xs, ys, "--", color="black", linewidth=1.2, label="log-log fit")
    rho, _ = spearmanr(ss, ds)
    ax.text(0.04, 0.96,
            f"n = {mask.sum()}\nρ = {rho:.3f}\nslope = {lr.slope:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=11,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", linewidth=0.6))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ssDNA PPM")
    ax.set_ylabel("dsDNA PPM")
    ax.set_title(f"{CLUSTER_LABEL[cl]} — R{r}", loc="center",
                 fontweight="bold", color=C_CLUSTER[cl])
    ax.legend(frameon=False, loc="lower right", fontsize=10)
    fname = f"scatter_{CLUSTER_SHORT_NAME[cl]}_R{r}.png"
    fig.savefig(PANELS / fname)
    plt.close(fig)
    print(f"wrote panels/{fname}")


def render_all_scatters(df: pd.DataFrame):
    for cl in [1, 2, 3]:
        for r in [3, 4]:
            render_scatter_cluster_round(df, cl, r)


# =========================================================
# Panel: clone trajectory contrast
# =========================================================
CLONES_OF_INTEREST = [
    ("GASGYSSGNIDA", "12 aa, C1 — ss-biased", "#9D0208"),
    ("SADSCATCATYPSEIDA", "17 aa, C1 — unbiased", "#003566"),
]


def render_clone_trajectory_contrast(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    linestyles = {"ss": "-", "ds": "--"}
    markers = {"ss": "o", "ds": "s"}
    for seq, anno, color in CLONES_OF_INTEREST:
        for lib in ["ss", "ds"]:
            row = df[(df["seq"] == seq) & (df["lib"] == lib)]
            if len(row) == 0:
                continue
            ppm = [row[f"PPM{r}"].iloc[0] for r in ROUNDS]
            ax.plot(ROUNDS, ppm, linestyles[lib], marker=markers[lib],
                    color=color, markersize=9, linewidth=2.4,
                    label=f"{seq} ({anno}) — {LIB_LABEL[lib]}")
    ax.set_yscale("log")
    ax.set_xticks(ROUNDS)
    ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("Selection round")
    ax.set_ylabel("PPM (log scale)")
    ax.set_title("Two C1 binders: one is strongly ss-biased, the other is not",
                 loc="center", fontweight="bold")
    ax.legend(frameon=False, loc="lower right", fontsize=10)
    fig.savefig(PANELS / "clone_trajectory_GASGYSSGNIDA_vs_SADSCATCATYPSEIDA.png")
    plt.close(fig)
    print("wrote panels/clone_trajectory_GASGYSSGNIDA_vs_SADSCATCATYPSEIDA.png")


if __name__ == "__main__":
    df, info = load_and_cluster()
    # methodology / overview panels
    tbl = make_read_depth_table(df)
    render_read_depth_table(df, tbl)
    for ds_name in DATASETS:
        render_enrichment_heatmap(ds_name, info[ds_name])
        render_cluster_trajectory(ds_name, info[ds_name])
    render_clonal_takeover_stack(df)
    # bias panels
    stats = compute_per_cluster_stats(df)
    render_spearman_by_cluster(stats)
    render_loglog_slope_by_cluster(stats)
    render_all_scatters(df)
    render_clone_trajectory_contrast(df)
    print("\nDone. Outputs:")
    print(f"  panels : {PANELS}")
    print(f"  data   : {DATA}")
