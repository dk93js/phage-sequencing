"""
Unified-style panels for Figure 2, Supplementary Figure 2, Supplementary Figure 3,
following the Figure-1 rules: A's font (13 / title 14 non-bold / label 14 / tick 12 /
legend 10), 3.6 in panel-height scale, lowercase case rule. Reuses make_panels for
data + clustering + per-cluster stats.

Panels (-> panels/*_sized.png):
  Fig 2 : scatter_{C1,C2,C3}_R{3,4}_sized  + cdrh3_length_by_round_sized (CORRECTED from v3 data)
  Supp 2: scatter_pooled_R{0..4}_sized  (pooled all-clones ss-vs-ds per round)
  Supp 3: morisita_horn_overlap_sized, shared_clone_cluster_confusion_sized,
          spearman_rho_by_cluster_sized, loglog_slope_by_cluster_sized

Then composes three 1:1 layout previews (not deliverables):
  fig2_preview.png, supp_fig2_preview.png, supp_fig3_preview.png
"""
from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import spearmanr, linregress

import make_panels as mp

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "lines.linewidth": 2.0, "lines.markersize": 6,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

PANELS = mp.PANELS
ROUNDS = mp.ROUNDS
PPM_COLS = mp.PPM_COLS
C_CLUSTER = mp.C_CLUSTER
CLUSTER_LABEL = mp.CLUSTER_LABEL
CLUSTER_SHORT = mp.CLUSTER_SHORT_NAME
PANEL_H = 3.6
C_SS = "#E07A1F"     # ssDNA (length panel) — Fig 1A orange (panelA_prototypes C_T1)
C_DS = "#3D5A80"     # dsDNA (length panel) — Fig 1A slate-blue (panelA_prototypes C_SH)
DPI = 300


# ============ scatter (Fig 2 A/B by cluster; Supp 2 pooled) ============
def render_scatter(df, r, cl=None):
    if cl is None:
        sub, color = df, "#555555"
        title, fname = f"all — R{r}", f"scatter_pooled_R{r}_sized.png"
    else:
        sub, color = df[df["cluster"] == cl], C_CLUSTER[cl]
        title, fname = f"{CLUSTER_LABEL[cl]} — R{r}", f"scatter_{CLUSTER_SHORT[cl]}_R{r}_sized.png"
    piv = sub.pivot_table(index="seq", columns="lib", values=f"PPM{r}", aggfunc="first").fillna(0)
    if "ss" not in piv or "ds" not in piv:
        return None
    mask = (piv["ss"] > 0) & (piv["ds"] > 0)
    ss, ds = piv.loc[mask, "ss"].to_numpy(), piv.loc[mask, "ds"].to_numpy()
    if len(ss) < 5:
        return None
    fig, ax = plt.subplots(figsize=(PANEL_H, PANEL_H))
    ax.scatter(ss, ds, s=11, color=color, alpha=0.55, edgecolor="none")
    lo, hi = min(ss.min(), ds.min()) * 0.7, max(ss.max(), ds.max()) * 1.4
    ax.plot([lo, hi], [lo, hi], color="red", linewidth=1.1, label="y = x")
    lr = linregress(np.log10(ss), np.log10(ds))
    xs = np.array([lo, hi])
    ax.plot(xs, 10 ** (lr.intercept + lr.slope * np.log10(xs)), "--", color="black",
            linewidth=1.1, label="log-log fit")
    rho, _ = spearmanr(ss, ds)
    ax.text(0.04, 0.96, f"n = {int(mask.sum())}\nρ = {rho:.3f}\nslope = {lr.slope:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", linewidth=0.6))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ssDNA clonal frequency"); ax.set_ylabel("dsDNA clonal frequency")
    ax.set_title(title, color=color)
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.savefig(PANELS / fname); plt.close(fig)
    return fname


# ============ CDRH3 length, corrected from v3 data (Fig 2C) ============
def render_length_panel(df):
    d = df.copy()
    d["L"] = d["seq"].str.len()
    lengths = list(range(7, 21))
    fig, axes = plt.subplots(1, 5, figsize=(12.5, 3.2), sharey=True)
    for r, ax in zip(ROUNDS, axes):
        for lib, color, lab in [("ss", C_SS, "ssDNA"), ("ds", C_DS, "dsDNA")]:
            sub = d[d.lib == lib]
            tot = sub[f"PPM{r}"].sum()
            frac = sub.groupby("L")[f"PPM{r}"].sum() / tot
            ax.plot(lengths, [frac.get(L, 0) for L in lengths], "o-", color=color,
                    markersize=4, linewidth=1.6, label=lab)
        ax.set_title(f"R{r}")
        ax.set_xlabel("HCDR3 length (aa)")
        if r == 0:
            ax.set_ylabel("fraction of reads")
            ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(PANELS / "cdrh3_length_by_round_sized.png"); plt.close(fig)
    print("wrote cdrh3_length_by_round_sized.png")


# ============ Morisita-Horn overlap (Supp 3A) ============
def _mh(x, y):
    X, Y = x.sum(), y.sum()
    if X == 0 or Y == 0:
        return np.nan
    return float(2 * np.sum(x * y) / ((np.sum(x ** 2) / X ** 2 + np.sum(y ** 2) / Y ** 2) * X * Y))


def render_morisita(df):
    piv = df.pivot_table(index="seq", columns="lib", values=PPM_COLS, aggfunc="first").fillna(0)
    mh = [_mh(piv[(f"PPM{r}", "ss")].to_numpy(), piv[(f"PPM{r}", "ds")].to_numpy()) for r in ROUNDS]
    fig, ax = plt.subplots(figsize=(4.2, PANEL_H))
    ax.plot(ROUNDS, mh, "o-", color="black", linewidth=2.2, markersize=7)
    ax.axhline(1.0, color="gray", lw=0.8, ls=":")
    ax.set_ylim(0, 1.05)
    ax.set_xticks(ROUNDS); ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("selection round")
    ax.set_ylabel("Morisita–Horn (abundance overlap)")
    ax.set_title("ss/ds overlap")
    fig.savefig(PANELS / "morisita_horn_overlap_sized.png"); plt.close(fig)
    print("wrote morisita_horn_overlap_sized.png  MH =", [round(v, 3) for v in mh])


# ============ by-cluster line (Supp 3 C spearman, D slope) ============
def render_line_by_cluster(stats, col, ylabel, title, ylim, fname):
    fig, ax = plt.subplots(figsize=(4.3, PANEL_H))
    for key, color, label in [("all", "black", "all"), ("C1", C_CLUSTER[1], "C1 enriching"),
                              ("C2", C_CLUSTER[2], "C2 intermediate"), ("C3", C_CLUSTER[3], "C3 depleted")]:
        s = stats[stats["cluster"] == key].sort_values("round")
        ax.plot(s["round"], s[col], "o-", color=color, label=label, linewidth=2.2, markersize=6)
    ax.axhline(1.0, color="gray", lw=0.8, ls=":")
    ax.set_ylim(*ylim)
    ax.set_xticks(ROUNDS); ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("selection round"); ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, loc="lower left", fontsize=9)
    fig.savefig(PANELS / fname); plt.close(fig)
    print(f"wrote {fname}")


# ============ confusion matrix (Supp 3B) ============
def render_confusion(df):
    cl = df.pivot_table(index="seq", columns="lib", values="cluster", aggfunc="first")
    sh = cl.dropna(subset=["ss", "ds"]).astype(int)
    M = np.zeros((3, 3), dtype=int)
    for a, b in zip(sh["ss"], sh["ds"]):
        M[a - 1, b - 1] += 1
    total, conc = M.sum(), int(np.trace(M))
    disc = total - conc
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    rf = M / M.sum(axis=1, keepdims=True)
    ax.imshow(rf, cmap="Blues", vmin=0, vmax=1)
    for i in range(3):
        for j in range(3):
            on = i == j
            ax.text(j, i, f"{M[i, j]}", ha="center", va="center",
                    fontsize=14 if on else 12, fontweight="bold" if on else "normal",
                    color="white" if rf[i, j] > 0.55 else "#222")
            if on:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor="#E63946", linewidth=2.2))
    ax.set_xticks(range(3)); ax.set_yticks(range(3))
    ax.set_xticklabels([CLUSTER_LABEL[k] for k in (1, 2, 3)], rotation=20, ha="right")
    ax.set_yticklabels([CLUSTER_LABEL[k] for k in (1, 2, 3)])
    for t, k in zip(ax.get_xticklabels(), (1, 2, 3)):
        t.set_color(C_CLUSTER[k])
    for t, k in zip(ax.get_yticklabels(), (1, 2, 3)):
        t.set_color(C_CLUSTER[k])
    ax.set_xlabel("dsDNA"); ax.set_ylabel("ssDNA")
    ax.set_title(f"shared clonotypes: ss vs ds cluster\n"
                 f"concordant {conc} ({100*conc/total:.0f}%), discordant {disc} ({100*disc/total:.0f}%)")
    ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", length=0)
    fig.savefig(PANELS / "shared_clone_cluster_confusion_sized.png"); plt.close(fig)
    print("wrote shared_clone_cluster_confusion_sized.png")


# ============ compose previews (native 1:1) ============
def compose(rows, out, gap_h=0.30, gap_v=0.55):
    imgs = [[mpimg.imread(PANELS / p) for p in row] for row in rows]
    win = lambda im: (im.shape[1] / DPI, im.shape[0] / DPI)
    row_w = [sum(win(im)[0] for im in row) + gap_h * (len(row) - 1) for row in imgs]
    row_h = [max(win(im)[1] for im in row) for row in imgs]
    FIG_W, FIG_H = max(row_w), sum(row_h) + gap_v * (len(rows) - 1)
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    y = FIG_H
    for ri, row in enumerate(imgs):
        y -= row_h[ri]
        x = (FIG_W - row_w[ri]) / 2
        for im in row:
            w, h = win(im)
            a = fig.add_axes([x / FIG_W, (y + (row_h[ri] - h)) / FIG_H, w / FIG_W, h / FIG_H])
            a.imshow(im); a.axis("off")
            x += w + gap_h
        y -= gap_v
    fig.savefig(PANELS / out, dpi=180, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    df, info = mp.load_and_cluster()
    stats = mp.compute_per_cluster_stats(df)

    # Fig 2 scatters (by cluster, R3 + R4)
    for r in (3, 4):
        for cl in (1, 2, 3):
            render_scatter(df, r, cl)
    render_length_panel(df)
    # Supp 2 pooled scatters (all clones, every round)
    for r in ROUNDS:
        render_scatter(df, r, cl=None)
    # Supp 3
    render_morisita(df)
    render_confusion(df)
    render_line_by_cluster(stats, "spearman", "Spearman ρ (ds vs ss)",
                           "rank correlation by cluster", (0.6, 1.02),
                           "spearman_rho_by_cluster_sized.png")
    render_line_by_cluster(stats, "slope_loglog", "log-log slope (ds vs ss)",
                           "frequency bias by cluster", (0.45, 1.1),
                           "loglog_slope_by_cluster_sized.png")

    # previews
    compose([["scatter_C1_enriching_R3_sized.png", "scatter_C2_intermediate_R3_sized.png", "scatter_C3_depleted_R3_sized.png"],
             ["scatter_C1_enriching_R4_sized.png", "scatter_C2_intermediate_R4_sized.png", "scatter_C3_depleted_R4_sized.png"],
             ["cdrh3_length_by_round_sized.png"]], "fig2_preview.png")
    compose([[f"scatter_pooled_R{r}_sized.png" for r in ROUNDS]], "supp_fig2_preview.png")
    compose([["morisita_horn_overlap_sized.png", "shared_clone_cluster_confusion_sized.png"],
             ["spearman_rho_by_cluster_sized.png", "loglog_slope_by_cluster_sized.png"]], "supp_fig3_preview.png")
    print("done")
