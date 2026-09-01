"""Fig 1B (dsDNA / ssDNA / pooled) with capitalised label text.

Copy of ``optionB/code/render_fig1B_promoted.py`` with the label strings as the
ONLY change, so geometry, palette, fonts and clustering stay byte-identical to
the deployed panels.

  * panel titles stay "ssDNA" / "dsDNA" -- these are nomenclature, not prose;
    "SsDNA" would be wrong. Only the pooled panel's prose title is capitalised.
  * axis labels take an initial capital: "Selection round",
    "Depth-normalized (z)".
  * legend entries already start with C1/C2/C3.

Pass --verify to re-render the ORIGINAL lower-case labels into
``_verify/`` and md5-compare against the deployed panels; that proves this
script reproduces the deck panels exactly before the capitalisation is applied.

Run:  python "figures/panels/new color/optionB/code/render_recolored.py" [--verify]
"""
import hashlib
import shutil
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PROJECT = HERE.parents[3]
OPTIONB_CODE = PROJECT / "figures" / "panels" / "new color" / "optionB" / "code"
DEPLOYED = OPTIONB_CODE.parent

sys.path.insert(0, str(OPTIONB_CODE))
import make_panels as mp            # noqa: E402
import make_panels_sized as mps     # noqa: E402

mpl.rcParams.update({"axes.titlesize": 14, "legend.fontsize": 10})

LABEL = {1: "C1 enriching", 2: "C2 intermediate: peaks at R2", 3: "C3 depleted"}
OUTNAME = {"pooled": "fig1B_trajectory_pooled.png",
           "dsDNA": "fig1B_trajectory_dsDNA.png",
           "ssDNA": "fig1B_trajectory_ssDNA.png"}

# --- the only difference from render_fig1B_promoted.py -----------------------
# Panel TITLE only. Axis labels stay exactly as published (lower case).
CAPS_TITLE = {"pooled": "All (ssDNA + dsDNA)", "ssDNA": "ssDNA", "dsDNA": "dsDNA"}
CAPS_XLABEL = "selection round"
CAPS_YLABEL = "depth-normalized (z)"

ORIG_TITLE = mp.DATASET_LABEL
ORIG_XLABEL = "selection round"
ORIG_YLABEL = "depth-normalized (z)"
# ----------------------------------------------------------------------------


def render(ds_name, info, outdir, title, xlabel, ylabel):
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
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / OUTNAME[ds_name])
    plt.close(fig)
    return counts


def md5(p):
    return hashlib.md5(p.read_bytes()).hexdigest()


df, info = mp.load_and_cluster()

if "--verify" in sys.argv:
    vdir = OUT / "_verify"
    print("reproducing the DEPLOYED (lower-case) panels for md5 comparison:")
    for ds in ("pooled", "dsDNA", "ssDNA"):
        render(ds, info[ds], vdir, ORIG_TITLE[ds], ORIG_XLABEL, ORIG_YLABEL)
        new, old = vdir / OUTNAME[ds], DEPLOYED / OUTNAME[ds]
        ok = old.exists() and md5(new) == md5(old)
        print(f"  {OUTNAME[ds]:32} {'MATCH' if ok else 'DIFFERS'}  {md5(new)[:10]}")
    shutil.rmtree(vdir, ignore_errors=True)   # ignore_errors: exFAT "._" forks
    print()

print("cluster sizes (n stays in the caption, not on the panel):")
for ds in ("pooled", "dsDNA", "ssDNA"):
    c = render(ds, info[ds], OUT, CAPS_TITLE[ds], CAPS_XLABEL, CAPS_YLABEL)
    print(f"  {CAPS_TITLE[ds]:22} C1 {c[1]:5}   C2 {c[2]:5}   C3 {c[3]:5}   total {sum(c.values()):5}")
print(f"\nwrote 3 panels to {OUT}")
