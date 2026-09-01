"""Fig 1C split into dsDNA / ssDNA panels, colored set = each library's TOP 1%.

Replaces the fixed top-20 selection of render_fig1C_split.py. Rationale: Fig 1A
already reports a "Top-1%" concentration metric, and using a different head size
(top-20) in Fig 1C made the two panels look contradictory to a reader comparing
them (Prof. Chung, 2026-08-11).

Selection rule here is IDENTICAL to Fig 1A's Top-1%, verified by reverse-
engineering panelA_prototypes.py's hardcoded series:
  - k = ceil(0.01 * n_clonotypes_in_that_library)  -> 11 for both libraries
  - the top k are re-selected WITHIN EACH ROUND (a concentration metric, not a
    fixed cohort)
  - denominator = that library's analysis-set read total for that round
    (i.e. frac = PPM / PPM.sum(axis=0)), the same denominator Fig 1C already used

Consequence: the top of the colored stack in each panel now reproduces Panel A's
per-library Top-1% curve exactly, not approximately:
  ssDNA  19.4  19.6  23.7  68.3  62.2
  dsDNA  19.1  19.3  24.4  44.7  37.5
The script asserts this against Panel A's published series before writing.

Because the cohort is re-selected per round, a given stripe is not the same
clonotype across rounds; the figure reports the composition of each round's head,
which is what the accompanying Results sentence claims.

Everything else (figure size, axes boxes, palette, titles, C2 inset, legend
wording) is unchanged from render_fig1C_split.py.

Run:  python "figures/panels/new color/optionB/code/render_recolored.py"
"""
import math
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
OPTIONB_CODE = PROJECT / "figures" / "panels" / "new color" / "optionB" / "code"

sys.path.insert(0, str(OPTIONB_CODE))
import make_panels as mp            # noqa: E402

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "savefig.dpi": 300, "savefig.bbox": None,
})

TITLE = {"dsDNA": "Clonal composition\ndsDNA", "ssDNA": "Clonal composition\nssDNA"}
OUTNAME = {"dsDNA": "fig1C_composition_dsDNA_top1pct.png",
           "ssDNA": "fig1C_composition_ssDNA_top1pct.png"}
ZOOM_ROUND = 2
PAD = 0.006

# Panel A's published Top-1% series (panelA_prototypes.py), used as a cross-check.
PANEL_A = {"ssDNA": [19.4, 19.6, 23.7, 68.3, 62.2],
           "dsDNA": [19.1, 19.3, 24.4, 44.7, 37.5]}


def round_segments(frac, clu, k):
    """Per round, the top-k clonotypes of that round, ordered C1 -> C2 -> C3.

    Returns segs[r] = list of (cluster, height), bottom-to-top.
    """
    segs = []
    for r in range(frac.shape[1]):
        idx = np.argsort(-frac[:, r])[:k]
        idx = sorted(idx, key=lambda i: (clu[i], -frac[i, r]))
        segs.append([(int(clu[i]), float(frac[i, r])) for i in idx])
    return segs


def draw(ax, segs):
    """Stack one library's per-round segments plus the grey remainder."""
    for r, seg in enumerate(segs):
        bottom = 0.0
        for c, h in seg:
            ax.bar(r, h, bottom=bottom, width=0.80, color=mp.C_CLUSTER[c],
                   edgecolor="white", linewidth=0.4)
            bottom += h
        ax.bar(r, 1 - bottom, bottom=bottom, width=0.80, color="#dcdcdc",
               edgecolor="white", linewidth=0.4)


def c2_span(segs, r):
    """y-range occupied by C2 in round r (lo, hi); (nan, nan) if C2 absent."""
    bottom, lo, hi = 0.0, None, None
    for c, h in segs[r]:
        if c == 2:
            if lo is None:
                lo = bottom
            hi = bottom + h
        bottom += h
    return (lo, hi) if lo is not None else (float("nan"), float("nan"))


def render(ds_name, sub, labels):
    clu = np.asarray(labels)
    ppm = sub[mp.PPM_COLS].to_numpy(dtype=float)
    frac = ppm / ppm.sum(axis=0, keepdims=True)
    k = math.ceil(0.01 * len(clu))

    segs = round_segments(frac, clu, k)

    totals = [100 * sum(h for _, h in seg) for seg in segs]
    ref = PANEL_A[ds_name]
    assert np.allclose(totals, ref, atol=0.06), \
        f"{ds_name}: colored total {totals} != Panel A Top-1% {ref}"

    lo, hi = c2_span(segs, ZOOM_ROUND)
    y0, y1 = lo - PAD, hi + PAD

    fig = plt.figure(figsize=(5.30, 3.6))
    ax = fig.add_axes([0.12, 0.14, 0.53, 0.72])
    axz = fig.add_axes([0.755, 0.28, 0.095, 0.44])

    for a in (ax, axz):
        draw(a, segs)

    ax.set_xticks(np.arange(5))
    ax.set_xticklabels([f"R{r}" for r in mp.ROUNDS])
    ax.set_xlabel("selection round")
    ax.set_ylabel("fraction of reads")
    ax.set_ylim(0, 1)
    ax.set_title(TITLE[ds_name])

    ax.add_patch(Rectangle((ZOOM_ROUND - 0.42, y0), 0.84, y1 - y0,
                           fill=False, edgecolor="#333333", linewidth=1.0, zorder=5))

    axz.set_xlim(ZOOM_ROUND - 0.42, ZOOM_ROUND + 0.42)
    axz.set_ylim(y0, y1)
    axz.set_xticks([])
    axz.set_yticks([])
    for s in ("top", "right", "bottom", "left"):
        axz.spines[s].set_visible(False)
    axz.add_patch(Rectangle((ZOOM_ROUND - 0.42, y0), 0.84, y1 - y0, fill=False,
                            edgecolor="#333333", linewidth=1.0, zorder=5))

    # inline cluster labels beside the inset + connector lines, matching the
    # deployed pooled Fig 1C (make_fig1C_inset.py) rather than a bottom legend
    labels = {1: "C1", 2: f"C2\n({(hi - lo) * 100:.1f}%)", 3: "C3"}
    for c, yc in ((1, (y0 + lo) / 2), (2, (lo + hi) / 2), (3, (hi + y1) / 2)):
        axz.annotate(labels[c], xy=(ZOOM_ROUND + 0.44, yc), xytext=(5, 0),
                     textcoords="offset points", va="center", ha="left",
                     fontsize=10, color=mp.C_CLUSTER[c], annotation_clip=False,
                     linespacing=1.35)

    for yy in (y0, y1):
        fig.add_artist(ConnectionPatch(
            xyA=(ZOOM_ROUND + 0.42, yy), coordsA=ax.transData,
            xyB=(ZOOM_ROUND - 0.42, yy), coordsB=axz.transData,
            color="#333333", linewidth=0.7, zorder=6))

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / OUTNAME[ds_name]
    fig.savefig(out)
    plt.close(fig)

    per_cluster = {}
    for c in (1, 2, 3):
        per_cluster[c] = ([100 * sum(h for cc, h in seg if cc == c) for seg in segs],
                          [sum(1 for cc, _ in seg if cc == c) for seg in segs])
    return out, k, totals, per_cluster


def main():
    df, info = mp.load_and_cluster()
    for ds_name in ("dsDNA", "ssDNA"):
        out, k, totals, per_cluster = render(
            ds_name, info[ds_name]["df"], info[ds_name]["labels"])
        print(f"wrote {out}")
        print(f"  top 1% = {k} clonotypes, re-selected per round")
        print("  colored total (== Panel A Top-1%): "
              + "  ".join(f"{v:.1f}" for v in totals))
        for c in (1, 2, 3):
            pct, n = per_cluster[c]
            print(f"    C{c}: " + "  ".join(f"{v:5.1f}%" for v in pct)
                  + "   n/round " + str(n))


if __name__ == "__main__":
    main()
