"""Fig 1A small multiples (titer | Shannon H | top-1%) with capitalised titles.

Copy of the v3 block of ``analysis/version 3/code/panelA_prototypes.py``.
Changes vs that block, and nothing else:

  1. panel titles take an initial capital
     "output/input (titer)"  -> "Output/input (titer)"
     "clonal diversity"      -> "Clonal diversity"
     "Top-1% of reads (%)"   -> "Read fraction of top 1% clones"
        (DK retitled this panel on 2026-08-02; the capital is the new part)
  2. output filename/dir

Axis labels, legend, colours, markers, limits, figsize and dpi are untouched.

  !! The numbers here are the ones hardcoded in panelA_prototypes.py, and they
  !! reproduce EXACTLY from report_q40_ppm100_list.csv (top 1% of clonotypes =
  !! 11 per library, cumulative read fraction):
  !!     ssDNA 19.4 19.6 23.7 68.3 62.2   dsDNA 19.1 19.3 24.4 44.7 37.5
  !! They match the Results text. The panel currently in the deck
  !! (fig1A_smallmultiples_retitled.png) plots a DIFFERENT, stale series
  !! (~7-8 early, ss ~34.5/36 at R3/R4) that no script on disk produces.
  !! See README.md.

Pass --verify to render with the pre-retitle titles and md5-compare against
``analysis/version 3/panels/panelA_v3_smallmultiples.png``.

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

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PROJECT = HERE.parents[3]
REFERENCE = PROJECT / "analysis" / "version 3" / "panels" / "panelA_v3_smallmultiples.png"
OUTNAME = "fig1A_smallmultiples_retitled.png"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 14,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.1, "lines.linewidth": 2.0, "lines.markersize": 6,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

R = ["R0", "R1", "R2", "R3", "R4"]
x = np.arange(5)
titer = np.array([np.nan, 2.6e-4, 2.7e-7, 4.5e-6, 4.6e-5])
shannon_ss = np.array([6.12, 6.12, 6.04, 3.65, 3.96])
shannon_ds = np.array([6.17, 6.15, 6.03, 4.88, 5.35])
top1_ss = np.array([19.4, 19.6, 23.7, 68.3, 62.2])
top1_ds = np.array([19.1, 19.3, 24.4, 44.7, 37.5])
shannon_all = np.array([6.16, 6.15, 6.05, 4.38, 4.69])
top1_all = np.array([19.9, 20.1, 24.5, 57.1, 51.4])

C_SH = "#3D5A80"
C_T1 = "#E07A1F"
C_TI = "#222222"
C_ALL = "#000000"

# --- the only difference from panelA_prototypes.py v3 ------------------------
CAPS = ("Output/input (titer)", "Clonal diversity", "Read fraction of top 1% clones")
ORIG = ("output/input (titer)", "clonal diversity", "Top-1% of reads (%)")
# -----------------------------------------------------------------------------


def render(titles, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.6))
    a0, a1, a2 = axes
    a0.plot(x, titer, "D-", color=C_TI)
    a0.set_yscale("log")
    a0.set_title(titles[0])
    a0.set_ylabel("ratio")
    a1.plot(x, shannon_ss, "o-", color=C_SH, label="ssDNA")
    a1.plot(x, shannon_ds, "s--", color=C_SH, mfc="white", label="dsDNA")
    a1.plot(x, shannon_all, "^-", color=C_ALL, label="all")
    a1.set_title(titles[1])
    a1.set_ylabel("Shannon H")
    a1.set_ylim(3, 7)
    a1.legend(frameon=False)
    a2.plot(x, top1_ss, "o-", color=C_T1, label="ssDNA")
    a2.plot(x, top1_ds, "s--", color=C_T1, mfc="white", label="dsDNA")
    a2.plot(x, top1_all, "^-", color=C_ALL, label="all")
    a2.set_title(titles[2])
    a2.set_ylabel("% of reads")
    a2.set_ylim(0, 75)
    a2.legend(frameon=False)
    for a in axes:
        a.set_xticks(x)
        a.set_xticklabels(R)
        a.set_xlabel("selection round")
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / OUTNAME
    fig.savefig(out)
    plt.close(fig)
    return out


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


if "--verify" in sys.argv:
    vdir = OUT / "_verify"
    out = render(ORIG, vdir)
    ok = REFERENCE.exists() and md5(out) == md5(REFERENCE)
    print("reproducing panelA_v3_smallmultiples.png (pre-retitle titles):")
    print(f"  {'panelA_v3_smallmultiples.png':44} {'MATCH' if ok else 'DIFFERS'}  {md5(out)[:10]}\n")
    shutil.rmtree(vdir, ignore_errors=True)   # ignore_errors: exFAT "._" forks

out = render(CAPS, OUT)
print(f"wrote {out}")
print("top-1% series plotted (matches Results text):")
print(f"  ssDNA {top1_ss}\n  dsDNA {top1_ds}\n  all   {top1_all}")
