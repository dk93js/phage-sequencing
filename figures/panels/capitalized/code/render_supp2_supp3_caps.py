"""Supp Fig 2 (per-round scatters) + Supp Fig 3 A-D with capitalised titles.

Uses the optionB copies of make_panels / make_figures_sized, i.e. the CURRENT
deck palette C1 #E63946 / C2 #17375E / C3 #56C1B0. Title strings and the output
directory are the only changes; the panel functions are monkey-patched at the
title line, everything else runs untouched.

  Supp 2   "all - R{r}"                     -> "All - R{r}"          (5 panels)
  Supp 3A  "ss/ds overlap"                  -> UNCHANGED (see below)
  Supp 3B  "shared clonotypes: ss vs ds..." -> "Shared clonotypes: ss vs ds..."
  Supp 3C  "rank correlation by cluster"    -> "Rank correlation by cluster"
  Supp 3D  "frequency bias by cluster"      -> "Frequency bias by cluster"

Supp 3A is the one title that cannot take an initial capital: it begins with the
nomenclature "ss/ds", and "Ss/ds overlap" would be wrong. Left as published --
decide separately whether to reword it (e.g. "Abundance overlap (ss/ds)").

!! Supp 3B: the panel currently in the deck is the OLD palette (pixel-identical
!! to figures/panels/new color/_current_control/suppfig3B_cluster_confusion.png).
!! 3C and 3D in the deck are the new palette. This script renders 3B with the
!! new palette, so swapping it in also closes that gap. See README.md.

Pass --verify to render the original titles and md5-compare against the
published panels.

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

sys.path.insert(0, str(OPTIONB_CODE))
import make_panels as mp            # noqa: E402
import make_figures_sized as mfs    # noqa: E402

# published locations, for --verify
PUBLISHED = {
    **{f"suppfig2_R{r}.png": PROJECT / "figures" / "panels" / f"suppfig2_R{r}.png"
       for r in range(5)},
    "suppfig3A_morisita_horn.png": PROJECT / "figures" / "panels" / "suppfig3A_morisita_horn.png",
    "suppfig3B_cluster_confusion.png": OPTIONB_CODE.parent / "suppfig3B_cluster_confusion.png",
    "suppfig3C_spearman_rho.png": OPTIONB_CODE.parent / "suppfig3C_spearman_rho.png",
    "suppfig3D_loglog_slope.png": OPTIONB_CODE.parent / "suppfig3D_loglog_slope.png",
}
RENAME = {**{f"scatter_pooled_R{r}_sized.png": f"suppfig2_R{r}.png" for r in range(5)},
          "morisita_horn_overlap_sized.png": "suppfig3A_morisita_horn.png",
          "shared_clone_cluster_confusion_sized.png": "suppfig3B_cluster_confusion.png",
          "spearman_rho_by_cluster_sized.png": "suppfig3C_spearman_rho.png",
          "loglog_slope_by_cluster_sized.png": "suppfig3D_loglog_slope.png"}

# --- the only difference: title strings --------------------------------------
TITLES = {
    "caps": {"supp2": "All — R{r}",
             "supp3b_head": "Shared clonotypes: ss vs ds cluster",
             "supp3c": "Rank correlation by cluster",
             "supp3d": "Frequency bias by cluster"},
    "orig": {"supp2": "all — R{r}",
             "supp3b_head": "shared clonotypes: ss vs ds cluster",
             "supp3c": "rank correlation by cluster",
             "supp3d": "frequency bias by cluster"},
}
# -----------------------------------------------------------------------------

_orig_scatter = mfs.render_scatter
_orig_confusion = mfs.render_confusion


def render(mode, outdir):
    """Run the published Supp 2 / Supp 3 block with `mode` titles."""
    T = TITLES[mode]
    outdir.mkdir(parents=True, exist_ok=True)
    mfs.PANELS = outdir
    mp.PANELS = outdir

    # Faithful route: let the stock panel functions run untouched and remap the
    # title string at ax.set_title, so nothing else can drift.
    import matplotlib.axes as maxes
    real_set_title = maxes.Axes.set_title
    remap = {
        TITLES["orig"]["supp2"].format(r=r): T["supp2"].format(r=r) for r in range(5)
    }
    remap[TITLES["orig"]["supp3c"]] = T["supp3c"]
    remap[TITLES["orig"]["supp3d"]] = T["supp3d"]

    def patched(self, label="", *a, **kw):
        if isinstance(label, str):
            if label in remap:
                label = remap[label]
            elif label.startswith(TITLES["orig"]["supp3b_head"]):
                label = T["supp3b_head"] + label[len(TITLES["orig"]["supp3b_head"]):]
        return real_set_title(self, label, *a, **kw)

    maxes.Axes.set_title = patched
    try:
        for r in mp.ROUNDS:
            mfs.render_scatter(df, r)                                    # Supp 2
        mfs.render_morisita(df)                                          # Supp 3A
        _orig_confusion(df)                                              # Supp 3B
        mfs.render_line_by_cluster(stats, "spearman", "Spearman ρ (ds vs ss)",
                                   TITLES["orig"]["supp3c"], (0.6, 1.02),
                                   "spearman_rho_by_cluster_sized.png")  # Supp 3C
        mfs.render_line_by_cluster(stats, "slope_loglog", "log-log slope (ds vs ss)",
                                   TITLES["orig"]["supp3d"], (0.45, 1.1),
                                   "loglog_slope_by_cluster_sized.png")  # Supp 3D
    finally:
        maxes.Axes.set_title = real_set_title

    made = []
    for native, keyed in RENAME.items():
        src = outdir / native
        if src.exists():
            shutil.move(src, outdir / keyed)
            made.append(outdir / keyed)
    return made


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


df, info = mp.load_and_cluster()
stats = mp.compute_per_cluster_stats(df)

if "--verify" in sys.argv:
    vdir = OUT / "_verify"
    print("reproducing the PUBLISHED (original titles) panels for md5 comparison:")
    for p in render("orig", vdir):
        old = PUBLISHED[p.name]
        ok = old.exists() and md5(p) == md5(old)
        note = "" if p.name != "suppfig3B_cluster_confusion.png" else "  <- new palette; deck has OLD"
        print(f"  {p.name:32} {'MATCH' if ok else 'DIFFERS'}  {md5(p)[:10]}{note}")
    shutil.rmtree(vdir, ignore_errors=True)   # ignore_errors: exFAT "._" forks
    print()

for p in render("caps", OUT):
    print(f"wrote {p.name}")
print(f"\npalette in use: C1 {mp.C_CLUSTER[1]}  C2 {mp.C_CLUSTER[2]}  C3 {mp.C_CLUSTER[3]}")
print("Supp 3A title left as published ('ss/ds overlap') - see module docstring.")
