"""Supp Fig 1B enrichment heatmaps (ssDNA / dsDNA / pooled) with capitalised titles.

Reproduces the Supp 1B block of ``optionB/code/render_recolored.py`` exactly --
same module, same rcParams restore (MP_RC), same MANIFEST renames -- and changes
ONLY the panel title strings.

  ssDNA / dsDNA          -> unchanged (nomenclature; "SsDNA" would be wrong)
  all (ssDNA + dsDNA)    -> All (ssDNA + dsDNA)

Axis labels ("selection round"), the colorbar label ("depth-normalized (z)"),
the "n = ..." annotation and the C1/C2/C3 strip labels are untouched.

PALETTE: this imports the optionB copy of make_panels, i.e. the CURRENT deck
palette C1 #E63946 / C2 #17375E / C3 #56C1B0. The canonical
``analysis/version 3/code/make_panels.py`` still carries the OLD palette
(C2 #457B9D / C3 #2A9D8F) -- do not render Supp 1B from there.

Pass --verify to render the original titles and md5-compare against the
deployed ``optionB/suppfig1B_heatmap_*.png``.

Run:  python "figures/panels/new color/optionB/code/render_recolored.py" [--verify]
"""
import hashlib
import shutil
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
PROJECT = HERE.parents[3]
OPTIONB_CODE = PROJECT / "figures" / "panels" / "new color" / "optionB" / "code"
DEPLOYED = OPTIONB_CODE.parent

sys.path.insert(0, str(OPTIONB_CODE))
import make_panels as mp  # noqa: E402

# make_panels' own rcParams -- the heatmaps were rendered under these
MP_RC = {"axes.titlesize": 16, "legend.fontsize": 11}

RENAME = {"enrichment_heatmap_ssDNA.png": "suppfig1B_heatmap_ssDNA.png",
          "enrichment_heatmap_dsDNA.png": "suppfig1B_heatmap_dsDNA.png",
          "enrichment_heatmap_pooled.png": "suppfig1B_heatmap_pooled.png"}

# --- the only difference from render_recolored.py ----------------------------
ORIG_LABEL = dict(mp.DATASET_LABEL)                       # all (ssDNA + dsDNA)
CAPS_LABEL = dict(ORIG_LABEL, pooled="All (ssDNA + dsDNA)")
# -----------------------------------------------------------------------------


def render(labels, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    mp.DATASET_LABEL = labels
    mp.PANELS = outdir                 # keep the optionB tree untouched
    mpl.rcParams.update(MP_RC)
    for ds in ("pooled", "ssDNA", "dsDNA"):
        mp.render_enrichment_heatmap(ds, info[ds])
    for native, keyed in RENAME.items():
        src = outdir / native
        shutil.move(src, outdir / keyed)
    return [outdir / k for k in RENAME.values()]


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


df, info = mp.load_and_cluster()

if "--verify" in sys.argv:
    vdir = OUT / "_verify"
    print("reproducing the DEPLOYED (original titles) heatmaps for md5 comparison:")
    for p in render(ORIG_LABEL, vdir):
        old = DEPLOYED / p.name
        ok = old.exists() and md5(p) == md5(old)
        print(f"  {p.name:32} {'MATCH' if ok else 'DIFFERS'}  {md5(p)[:10]}")
    shutil.rmtree(vdir, ignore_errors=True)   # ignore_errors: exFAT "._" forks
    print()

for p in render(CAPS_LABEL, OUT):
    print(f"wrote {p.name}")
print(f"\npalette in use: C1 {mp.C_CLUSTER[1]}  C2 {mp.C_CLUSTER[2]}  C3 {mp.C_CLUSTER[3]}")
