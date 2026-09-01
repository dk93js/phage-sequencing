"""Fig 1C (clonal composition + C2 inset) with capitalised label text.

Copy of ``analysis/version 3/code/make_fig1C_inset.py`` with the label strings
and the output directory as the ONLY changes, so the stack, the inset geometry
and the palette stay identical to the deployed panel.

  title   "clonal composition" -> "Clonal composition"
  xlabel  "selection round"    -> "Selection round"
  ylabel  "fraction of reads"  -> "Fraction of reads"

Pass --verify to re-render the ORIGINAL lower-case labels into ``_verify/``
and md5-compare against the deployed panel.

Run:  python "figures/panels/new color/optionB/code/render_recolored.py" [--verify]
"""
import hashlib
import shutil
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PROJECT = HERE.parents[3]
V3_CODE = PROJECT / "analysis" / "version 3" / "code"
DEPLOYED = V3_CODE.parent / "panels" / "clonal_takeover_stack_pooled_sized_inset.png"

sys.path.insert(0, str(V3_CODE))
import make_panels as mp  # noqa: E402

PALETTE = {1: "#E63946", 2: "#17375E", 3: "#56C1B0"}
GREY = "#dcdcdc"
ZOOM_ROUND = 2          # R2
PAD = 0.006             # y-padding around the C2 span, in fraction-of-reads units
OUTNAME = "clonal_takeover_stack_pooled_sized_inset.png"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "savefig.dpi": 300, "savefig.bbox": None,
})

# --- the only difference from make_fig1C_inset.py ---------------------------
# Panel TITLE only. Axis labels stay exactly as published (lower case).
CAPS = {"title": "Clonal composition", "x": "selection round", "y": "fraction of reads"}
ORIG = {"title": "clonal composition", "x": "selection round", "y": "fraction of reads"}
# ----------------------------------------------------------------------------


def build(top_n: int = 20):
    df, _ = mp.load_and_cluster()
    ppm = df[mp.PPM_COLS].to_numpy(dtype=float)
    frac = ppm / ppm.sum(axis=0, keepdims=True)
    clu = df["cluster"].to_numpy()
    top = np.argsort(-frac.max(axis=1))[:top_n]
    order = sorted(top, key=lambda i: (clu[i], -frac[i, 4]))
    return frac, clu, order


def render(frac, clu, order, outdir, text):
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
    ax.set_xlabel(text["x"])
    ax.set_ylabel(text["y"])
    ax.set_ylim(0, 1)
    ax.set_title(text["title"])

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

    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / OUTNAME
    fig.savefig(out)
    plt.close(fig)
    return out, lo, hi, y1 - y0


def md5(p):
    return hashlib.md5(p.read_bytes()).hexdigest()


def main():
    frac, clu, order = build()

    if "--verify" in sys.argv:
        vdir = OUT / "_verify"
        out, *_ = render(frac, clu, order, vdir, ORIG)
        ok = DEPLOYED.exists() and md5(out) == md5(DEPLOYED)
        print("reproducing the DEPLOYED (lower-case) panel for md5 comparison:")
        print(f"  {OUTNAME:44} {'MATCH' if ok else 'DIFFERS'}  {md5(out)[:10]}\n")
        shutil.rmtree(vdir, ignore_errors=True)   # ignore_errors: exFAT "._" forks

    out, lo, hi, span = render(frac, clu, order, OUT, CAPS)
    print(f"wrote {out}")
    print(f"C2 span at R{ZOOM_ROUND}: {lo*100:.2f}% -> {hi*100:.2f}% "
          f"(thickness {(hi-lo)*100:.2f}%), magnification x{1/span:.0f}")


if __name__ == "__main__":
    main()
