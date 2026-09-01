"""
Font-matched re-render of Figure 1 panels B and C so their text matches Panel A.

Panel 1A (panelA_v3_smallmultiples) was approved as-is; it uses:
    font 13, axes.titlesize 14 (NON-bold), axes.labelsize 14,
    xtick 12, ytick 12, legend 10.
make_panels.py renders B/C with titlesize 16 + BOLD titles + legend 11, which
reads larger when the panels are placed next to A.

This script re-renders ONLY those two panels with A's settings, writing:
    panels/cluster_trajectory_pooled_sized.png      (Fig 1B)
    panels/clonal_takeover_stack_pooled_sized.png   (Fig 1C)

Content, data and normalization are identical to make_panels.py — the only
changes are font sizes and title weight. Labels are reproduced to match the
assembled figure ("ssDNA + dsDNA"; C2 legend notes "peaks at R2").
"""
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import make_panels as mp  # reuse data loading, clustering, constants, helpers

# ---- Panel A's approved font settings (verbatim from panelA_prototypes.py) ----
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "lines.linewidth": 2.0, "lines.markersize": 6,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

ROUNDS = mp.ROUNDS
PPM_COLS = mp.PPM_COLS
C_CLUSTER = mp.C_CLUSTER
CLUSTER_LABEL = mp.CLUSTER_LABEL
PANELS = mp.PANELS

# --- Match 1A's SCALE, not just its font spec ---
# panelA_v3 is figsize (11, 3.6) for three panels -> each panel is 3.6 in tall.
# Render B and C at that same 3.6 in height so, placed at equal zoom next to 1A,
# their text is the same physical size. Width follows each panel's own aspect
# ratio (so the plot shapes are preserved), per DK's "match the ratio".
A_PANEL_H = 3.6
B_FIGSIZE = (A_PANEL_H * 5.6 / 5.0, A_PANEL_H)   # B keeps its 5.6:5.0 aspect -> (4.03, 3.6)
C_FIGSIZE = (A_PANEL_H * 4.8 / 5.0, A_PANEL_H)   # C keeps its 4.8:5.0 aspect -> (3.46, 3.6)


def render_B_sized(info: dict):
    """Fig 1B — pooled (ss+ds) cluster trajectory, A-matched fonts."""
    z = info["z"]
    labels = info["labels"]
    legend_label = {
        1: f"C1 enriching (n={int((labels == 1).sum())})",
        2: f"C2 intermediate: peaks at R2 (n={int((labels == 2).sum())})",
        3: f"C3 depleted (n={int((labels == 3).sum())})",
    }
    fig, ax = plt.subplots(figsize=B_FIGSIZE)
    for k in [1, 2, 3]:
        sel = labels == k
        if not sel.any():
            continue
        m = z[sel].mean(axis=0)
        s = z[sel].std(axis=0)
        ax.plot(ROUNDS, m, "o-", color=C_CLUSTER[k], label=legend_label[k],
                linewidth=2.4, markersize=8)
        ax.fill_between(ROUNDS, m - s, m + s, color=C_CLUSTER[k], alpha=0.18, linewidth=0)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle=":")
    ax.set_xticks(ROUNDS)
    ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("selection round")
    ax.set_ylabel("depth-normalized (z)")   # values are per-clone z-scores
    ax.set_title("all (ssDNA + dsDNA)")    # non-bold, size 14 -> matches 1A
    # legend outside (upper-right) so it no longer overlaps the C2/C3 lines
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.savefig(PANELS / "cluster_trajectory_pooled_sized.png")
    plt.close(fig)
    print("wrote panels/cluster_trajectory_pooled_sized.png")


def render_C_sized(df, top_n: int = 20):
    """Fig 1C — clonal composition stack, A-matched fonts.

    Normalization is IDENTICAL to make_panels.render_clonal_takeover_stack
    (per-round renormalization over all rows, both libraries). Only fonts
    and title weight change here.
    """
    ppm = df[PPM_COLS].to_numpy(dtype=float)
    frac = ppm / ppm.sum(axis=0, keepdims=True)
    clu = df["cluster"].to_numpy()

    peak = frac.max(axis=1)
    top_idx = np.argsort(-peak)[:top_n]
    top_sorted = sorted(top_idx, key=lambda i: (clu[i], -frac[i, 4]))

    shade = {}
    for c in (1, 2, 3):
        members = [i for i in top_sorted if clu[i] == c]
        for rank, i in enumerate(members):
            # RECOLOUR (option B): narrower lighten span — a light C3 base would
            # otherwise fade into the grey "other" band at the pale end.
            f = 0.10 + 0.30 * (rank / max(1, len(members) - 1)) if len(members) > 1 else 0.15
            shade[i] = mp._lighten(C_CLUSTER[c], f)

    fig, ax = plt.subplots(figsize=C_FIGSIZE)
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
    ax.set_xlabel("selection round")
    ax.set_ylabel("fraction of reads")
    ax.set_ylim(0, 1)
    ax.set_title("clonal composition")     # non-bold, size 14 -> matches 1A

    # Legend and worked-clone (GASGYSSGNIDA) annotation removed per DK;
    # cluster colors are defined by the shared legend in Panel B.
    fig.savefig(PANELS / "clonal_takeover_stack_pooled_sized.png")
    plt.close(fig)
    print("wrote panels/clonal_takeover_stack_pooled_sized.png")


def render_supp_trajectory_sized(ds_name: str, info: dict, short_labels: bool):
    """Supp Fig 1C — per-library (ss/ds) cluster trajectory at 1A's font + scale.

    Same font rule (A's rcParams, non-bold title) and scale rule (3.6 in height)
    as Fig 1B, plus the lowercase case rule. ssDNA keeps descriptive cluster
    labels; dsDNA uses terse C1/C2/C3 (the locked-in legend style). y-axis matches
    the Fig 1B trajectory ("depth-normalized (z)"), not the heatmap colorbar.
    """
    z = info["z"]
    labels = info["labels"]
    fig, ax = plt.subplots(figsize=B_FIGSIZE)
    for k in [1, 2, 3]:
        sel = labels == k
        if not sel.any():
            continue
        m = z[sel].mean(axis=0)
        s = z[sel].std(axis=0)
        if k == 2:
            label = f"C2 intermediate: peaks at R2 (n={int(sel.sum())})"
        else:
            label = (f"C{k} (n={int(sel.sum())})" if short_labels
                     else f"{CLUSTER_LABEL[k]} (n={int(sel.sum())})")
        ax.plot(ROUNDS, m, "o-", color=C_CLUSTER[k], label=label,
                linewidth=2.4, markersize=8)
        ax.fill_between(ROUNDS, m - s, m + s, color=C_CLUSTER[k], alpha=0.18, linewidth=0)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle=":")
    ax.set_xticks(ROUNDS)
    ax.set_xticklabels([f"R{r}" for r in ROUNDS])
    ax.set_xlabel("selection round")
    ax.set_ylabel("depth-normalized (z)")
    ax.set_title(ds_name)                  # "ssDNA"/"dsDNA", non-bold 14 -> matches 1A
    # legend outside (upper-right) — same fix as Fig 1B; keeps ss/ds symmetric
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.savefig(PANELS / f"cluster_trajectory_{ds_name}_sized.png")
    plt.close(fig)
    print(f"wrote panels/cluster_trajectory_{ds_name}_sized.png")


if __name__ == "__main__":
    df, info = mp.load_and_cluster()
    render_B_sized(info["pooled"])
    render_C_sized(df)
    # Supp Fig 1C — per-library trajectories (ssDNA descriptive, dsDNA terse)
    render_supp_trajectory_sized("ssDNA", info["ssDNA"], short_labels=False)
    render_supp_trajectory_sized("dsDNA", info["dsDNA"], short_labels=True)
    print("done — _sized panels written to panels/")
