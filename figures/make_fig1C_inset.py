"""
Fig 1C with a magnified inset that resolves the C2 (intermediate) bands.

C2 is genuinely thin — it never holds more than ~1.7% of reads in any round
(R0 1.39, R1 1.71, R2 1.69, R3 0.63, R4 0.87 %) — so at full scale its three
clonotypes are only a few pixels tall. The inset magnifies the C1/C2/C3
boundary at R2 (the round named in the C2 label, "peaks at R2") so the three
C2 clonotypes are visible without distorting the stack itself.

Colors are the deck palette (C1 #E63946, C2 #17375E, C3 #56C1B0), drawn solid
(no lightening) so Fig 1C matches Fig 1B exactly.

Run with /opt/anaconda3/bin/python from this directory.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_panels as mp  # noqa: E402

PALETTE = {1: "#E63946", 2: "#17375E", 3: "#56C1B0"}
GREY = "#dcdcdc"
ZOOM_ROUND = 2          # R2
PAD = 0.006             # y-padding around the C2 span, in fraction-of-reads units

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "savefig.dpi": 300, "savefig.bbox": None,
})


def build(top_n: int = 20):
    df, _ = mp.load_and_cluster()
    ppm = df[mp.PPM_COLS].to_numpy(dtype=float)
    frac = ppm / ppm.sum(axis=0, keepdims=True)
    clu = df["cluster"].to_numpy()
    top = np.argsort(-frac.max(axis=1))[:top_n]
    order = sorted(top, key=lambda i: (clu[i], -frac[i, 4]))
    return frac, clu, order


def main():
    frac, clu, order = build()
    x = np.arange(5)

    # cumulative bottoms, and the y-span occupied by C2 in the zoom round
    bottoms, cum = {}, np.zeros(5)
    for i in order:
        bottoms[i] = cum.copy()
        cum = cum + frac[i]
    c2 = [i for i in order if clu[i] == 2]
    lo = bottoms[c2[0]][ZOOM_ROUND]
    hi = bottoms[c2[-1]][ZOOM_ROUND] + frac[c2[-1], ZOOM_ROUND]
    y0, y1 = lo - PAD, hi + PAD

    fig = plt.figure(figsize=(5.30, 3.6))
    ax = fig.add_axes([0.12, 0.14, 0.53, 0.78])
    axz = fig.add_axes([0.755, 0.30, 0.095, 0.44])

    for a in (ax, axz):
        bot = np.zeros(5)
        for i in order:
            a.bar(x, frac[i], bottom=bot, width=0.80, color=PALETTE[clu[i]],
                  edgecolor="white", linewidth=0.4)
            bot += frac[i]
        a.bar(x, 1 - bot, bottom=bot, width=0.80, color=GREY,
              edgecolor="white", linewidth=0.4)

    ax.set_xticks(x)
    ax.set_xticklabels([f"R{r}" for r in range(5)])
    ax.set_xlabel("selection round")
    ax.set_ylabel("fraction of reads")
    ax.set_ylim(0, 1)
    ax.set_title("clonal composition")

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

    labels = {1: "C1", 2: f"C2\n({(hi-lo)*100:.1f}%)", 3: "C3"}
    for c, yc in ((1, (y0 + lo) / 2), (2, (lo + hi) / 2), (3, (hi + y1) / 2)):
        axz.annotate(labels[c], xy=(ZOOM_ROUND + 0.44, yc), xytext=(5, 0),
                     textcoords="offset points", va="center", ha="left",
                     fontsize=10, color=PALETTE[c], annotation_clip=False,
                     linespacing=1.35)

    for yy, ya in ((y0, y0), (y1, y1)):
        con = ConnectionPatch(xyA=(ZOOM_ROUND + 0.42, ya), coordsA=ax.transData,
                              xyB=(ZOOM_ROUND - 0.42, yy), coordsB=axz.transData,
                              color="#333333", linewidth=0.7, zorder=6)
        fig.add_artist(con)

    out = mp.PANELS / "clonal_takeover_stack_pooled_sized_inset.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"C2 span at R{ZOOM_ROUND}: {lo*100:.2f}% -> {hi*100:.2f}% "
          f"(thickness {(hi-lo)*100:.2f}%), magnification x{1/(y1-y0):.0f}")


if __name__ == "__main__":
    main()
