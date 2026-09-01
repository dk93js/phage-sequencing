"""Render every cluster-coloured panel with this option's palette.

Run:  python "figures/panels/new color/optionB/code/render_recolored.py"

Covers the 16 panels that carry a C1/C2/C3 colour. Panels that do not (Fig 1A,
Fig 2C, all of Fig 3, Supp 2, Supp 4) are deliberately NOT regenerated.
"""
import shutil
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_panels as mp
import make_panels_sized as mps
import make_figures_sized as mfs

OUT = HERE.parent          # option root — figure-keyed copies land here
RENAME = {'cluster_trajectory_pooled_sized.png': 'fig1B_cluster_trajectory_pooled.png', 'clonal_takeover_stack_pooled_sized.png': 'fig1C_clonal_composition.png', 'scatter_C1_enriching_R3_sized.png': 'fig2A_C1_enriching_R3.png', 'scatter_C2_intermediate_R3_sized.png': 'fig2A_C2_intermediate_R3.png', 'scatter_C3_depleted_R3_sized.png': 'fig2A_C3_depleted_R3.png', 'scatter_C1_enriching_R4_sized.png': 'fig2B_C1_enriching_R4.png', 'scatter_C2_intermediate_R4_sized.png': 'fig2B_C2_intermediate_R4.png', 'scatter_C3_depleted_R4_sized.png': 'fig2B_C3_depleted_R4.png', 'enrichment_heatmap_ssDNA.png': 'suppfig1B_heatmap_ssDNA.png', 'enrichment_heatmap_dsDNA.png': 'suppfig1B_heatmap_dsDNA.png', 'enrichment_heatmap_pooled.png': 'suppfig1B_heatmap_pooled.png', 'cluster_trajectory_ssDNA_sized.png': 'suppfig1C_trajectory_ssDNA.png', 'cluster_trajectory_dsDNA_sized.png': 'suppfig1C_trajectory_dsDNA.png', 'shared_clone_cluster_confusion_sized.png': 'suppfig3B_cluster_confusion.png', 'spearman_rho_by_cluster_sized.png': 'suppfig3C_spearman_rho.png', 'loglog_slope_by_cluster_sized.png': 'suppfig3D_loglog_slope.png'}

# make_panels' own rcParams (heatmaps were rendered under these); the *_sized
# modules overwrite them at import, so restore before the heatmaps.
MP_RC = {"axes.titlesize": 16, "legend.fontsize": 11}
SIZED_RC = {"axes.titlesize": 14, "legend.fontsize": 10}

df, info = mp.load_and_cluster()
stats = mp.compute_per_cluster_stats(df)

# ---- Supp Fig 1B: enrichment heatmaps (cluster annotation strip) ----
mpl.rcParams.update(MP_RC)
for ds in ("pooled", "ssDNA", "dsDNA"):
    mp.render_enrichment_heatmap(ds, info[ds])

# ---- Fig 1B/1C + Supp Fig 1C ----
mpl.rcParams.update(SIZED_RC)
mps.render_B_sized(info["pooled"])
mps.render_C_sized(df)
mps.render_supp_trajectory_sized("ssDNA", info["ssDNA"], short_labels=False)
mps.render_supp_trajectory_sized("dsDNA", info["dsDNA"], short_labels=True)

# ---- Fig 2A/2B scatters ----
for r in (3, 4):
    for cl in (1, 2, 3):
        mfs.render_scatter(df, r, cl)

# ---- Supp Fig 3 B/C/D ----
mfs.render_confusion(df)
mfs.render_line_by_cluster(stats, "spearman", "Spearman \u03c1 (ds vs ss)",
                           "rank correlation by cluster", (0.6, 1.02),
                           "spearman_rho_by_cluster_sized.png")
mfs.render_line_by_cluster(stats, "slope_loglog", "log-log slope (ds vs ss)",
                           "frequency bias by cluster", (0.45, 1.1),
                           "loglog_slope_by_cluster_sized.png")

# ---- publish under the MANIFEST figure-keyed names ----
missing = []
for native, keyed in RENAME.items():
    src = mp.PANELS / native
    if src.exists():
        shutil.copy2(src, OUT / keyed)
    else:
        missing.append(native)
print(f"published {len(RENAME) - len(missing)}/{len(RENAME)} panels to {OUT}")
if missing:
    sys.exit("MISSING: " + ", ".join(missing))
