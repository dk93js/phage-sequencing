"""Fig 1C split into dsDNA / ssDNA panels (mirrors Fig 1B's per-library split).

The deployed Fig 1C (render_fig1C_caps.py) pools ssDNA+dsDNA into one stack,
which can't show the ssDNA-vs-dsDNA divergence the paper is actually about.
This renders two per-library stacks instead: each library's own top-20 most
abundant HCDR3 clonotypes (by that library's own PPM), colored by that
library's own Ward cluster (same clustering already used for the dsDNA/ssDNA
panels in Fig 1B, so cluster sizes here match the Fig 1B caption: C1/C2/C3 =
262/191/640 for dsDNA, 160/75/817 for ssDNA). Fraction of reads is normalized
within each library separately, so the top of the colored stack in each panel
reproduces Panel A's per-library Top-1% curve.

Run:  python "figures/panels/new color/optionB/code/render_recolored.py"
"""
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PROJECT = HERE.parents[3]
OPTIONB_CODE = PROJECT / "figures" / "panels" / "new color" / "optionB" / "code"

import sys
sys.path.insert(0, str(OPTIONB_CODE))
import make_panels as mp            # noqa: E402
import make_panels_sized as mps     # noqa: E402

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "savefig.dpi": 300, "savefig.bbox": None,
})

TITLE = {"dsDNA": "Clonal composition\ndsDNA", "ssDNA": "Clonal composition\nssDNA"}
OUTNAME = {"dsDNA": "fig1C_composition_dsDNA.png", "ssDNA": "fig1C_composition_ssDNA.png"}
ZOOM_ROUND = 2   # R2, same round zoomed in the pooled version
PAD = 0.006


def render(ds_name, sub, labels, top_n=20):
    ppm = sub[mp.PPM_COLS].to_numpy(dtype=float)
    frac = ppm / ppm.sum(axis=0, keepdims=True)
    clu = labels

    top_idx = np.argsort(-frac.max(axis=1))[:top_n]
    order = sorted(top_idx, key=lambda i: (clu[i], -frac[i, 4]))

    shade = {i: mp.C_CLUSTER[clu[i]] for i in order}

    x = np.arange(5)

    # cumulative bottoms, and the y-span occupied by C2 at the zoom round
    bottoms, cum = {}, np.zeros(5)
    for i in order:
        bottoms[i] = cum.copy()
        cum = cum + frac[i]
    c2 = [i for i in order if clu[i] == 2]
    lo = bottoms[c2[0]][ZOOM_ROUND]
    hi = bottoms[c2[-1]][ZOOM_ROUND] + frac[c2[-1], ZOOM_ROUND]
    y0, y1 = lo - PAD, hi + PAD

    fig = plt.figure(figsize=(5.30, 3.6))
    ax = fig.add_axes([0.12, 0.14, 0.53, 0.72])
    axz = fig.add_axes([0.755, 0.28, 0.095, 0.44])

    for a in (ax, axz):
        bottom = np.zeros(5)
        for i in order:
            a.bar(x, frac[i], bottom=bottom, width=0.80, color=shade[i],
                  edgecolor="white", linewidth=0.4)
            bottom += frac[i]
        a.bar(x, 1 - bottom, bottom=bottom, width=0.80, color="#dcdcdc",
              edgecolor="white", linewidth=0.4)

    ax.set_xticks(x)
    ax.set_xticklabels([f"R{r}" for r in mp.ROUNDS])
    ax.set_xlabel("selection round")
    ax.set_ylabel("fraction of reads")
    ax.set_ylim(0, 1)
    ax.set_title(TITLE[ds_name])

    # zoom window marker on the main axes
    ax.add_patch(Rectangle((ZOOM_ROUND - 0.42, y0), 0.84, y1 - y0,
                           fill=False, edgecolor="#333333", linewidth=1.0, zorder=5))

    # inset: same bars, y restricted to the C2 span at ZOOM_ROUND
    axz.set_xlim(ZOOM_ROUND - 0.42, ZOOM_ROUND + 0.42)
    axz.set_ylim(y0, y1)
    axz.set_xticks([])
    axz.set_yticks([])
    for s in ("top", "right", "bottom", "left"):
        axz.spines[s].set_visible(False)
    axz.add_patch(Rectangle((ZOOM_ROUND - 0.42, y0), 0.84, y1 - y0, fill=False,
                            edgecolor="#333333", linewidth=1.0, zorder=5))

    # legend, same wording/order as Fig 1B, with this library's C2 zoom span appended
    c2_pct = sum(frac[i, ZOOM_ROUND] for i in c2) * 100
    LEGEND_LABEL = {
        1: "C1 enriching",
        2: f"C2 intermediate: peaks at R2 ({c2_pct:.1f}% at R2)",
        3: "C3 depleted",
    }
    handles = [Patch(facecolor=mp.C_CLUSTER[c], edgecolor="none", label=LEGEND_LABEL[c])
               for c in (1, 2, 3)]
    fig.legend(handles=handles, frameon=False, loc="lower left",
               bbox_to_anchor=(0.755, 0.02), fontsize=10)

    for yy, ya in ((y0, y0), (y1, y1)):
        con = ConnectionPatch(xyA=(ZOOM_ROUND + 0.42, ya), coordsA=ax.transData,
                              xyB=(ZOOM_ROUND - 0.42, yy), coordsB=axz.transData,
                              color="#333333", linewidth=0.7, zorder=6)
        fig.add_artist(con)

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / OUTNAME[ds_name]
    fig.savefig(out)
    plt.close(fig)

    counts = {c: int((clu == c).sum()) for c in (1, 2, 3)}
    top_counts = {c: sum(1 for i in order if clu[i] == c) for c in (1, 2, 3)}
    return out, counts, top_counts


def main():
    df, info = mp.load_and_cluster()
    for ds_name in ("dsDNA", "ssDNA"):
        sub = info[ds_name]["df"]
        labels = info[ds_name]["labels"]
        out, counts, top_counts = render(ds_name, sub, labels)
        print(f"wrote {out}")
        print(f"  cluster sizes (library): C1 {counts[1]}  C2 {counts[2]}  C3 {counts[3]}  total {sum(counts.values())}")
        print(f"  of the top-20 shown:     C1 {top_counts[1]}  C2 {top_counts[2]}  C3 {top_counts[3]}")


if __name__ == "__main__":
    main()
