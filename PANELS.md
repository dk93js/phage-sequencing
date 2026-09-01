# Which script makes which panel

Every panel of every figure in the paper, with the script that produces it and the file
it writes. Run `bash run_all_figures.sh` from the repository root to regenerate all of
them in dependency order; 42 of the 44 images in the submitted figure deck come out
byte-identical (the other two are the arrow icons drawn in PowerPoint).

Two rendering stages produce most panels. The first writes the panels with the final
cluster palette (C1 `#E63946`, C2 `#17375E`, C3 `#56C1B0`); the second re-renders the
subset whose panel titles were capitalised for the submitted version, changing the title
text only. `figures/panels/capitalized/code/README_capitalisation.md` records exactly
which strings changed and includes a `--verify` mode that re-renders with the original
titles and checks the result against the published files.

## Figure 1 — sequencing of the two templates

| panel | script | output |
|---|---|---|
| A | `analysis/version 3/code/panelA_prototypes.py`, then `figures/panels/capitalized/code/render_fig1A_caps.py` | `fig1A_smallmultiples_retitled.png` |
| B | `figures/panels/new color/optionB/code/render_fig1B_promoted.py`, then `.../capitalized/code/render_fig1B_caps.py` | `fig1B_trajectory_dsDNA.png`, `fig1B_trajectory_ssDNA.png`, `fig1B_trajectory_pooled.png` |
| C | `figures/panels/capitalized/code/render_fig1C_caps.py` | `clonal_takeover_stack_pooled_sized_inset.png` |

`render_fig1C_split.py` and `render_fig1C_split_top1pct.py` produce the alternative
per-library version of panel C (one stack for ssDNA and one for dsDNA); the top-1%
variant is the one whose numbers are quoted in the Results (C1 reaches 65.6% of ssDNA
reads at R3 and 36.7% of dsDNA reads).

Panel A's summary values (titer, Shannon H, top-1% fraction) are written into
`panelA_prototypes.py` as literals; they are reproduced from the read table by
`figures/panels/capitalized/code/render_fig1A_caps.py`, which prints them on each run.

## Figure 2 — frequency bias by cluster

| panel | script | output |
|---|---|---|
| A, B | `figures/panels/new color/optionB/code/make_panels.py` → `make_panels_sized.py` → `render_recolored.py` | `fig2A_C{1,2,3}_*_R3.png`, `fig2B_C{1,2,3}_*_R4.png` |
| C | `figures/panels/new color/optionB/code/make_panels_sized.py` | `panels/cdrh3_length_by_round_sized.png` |

## Figure 3 — antigen-reactive clonotypes

| panel | script | output |
|---|---|---|
| A | `analysis/fig3_mapping/make_fig3_sized.py` | `fig3A_scatter_elisa_R{0..4}_sized.png` |
| B | `analysis/fig3_mapping/make_logratio_panel.py` | `analysis/version 3/panels/logratio_by_elisa_R3_sized.png` |

## Supplementary Figure 1 — sequencing statistics, enrichment, choice of k

| panel | script | output |
|---|---|---|
| A | table typed from `data/read_depth_per_round.csv` | — |
| B | `figures/panels/capitalized/code/render_supp1B_caps.py` | `suppfig1B_heatmap_{dsDNA,ssDNA,pooled}.png` |
| C, D, E | `analysis/revision_2026-08/r2_02_k_selection.py` → `render_revision_panels.py` | `out/figures/suppfig1{C,D,E}_*.png` |

## Supplementary Figure 2 — per-round frequency correlation

| panel | script | output |
|---|---|---|
| R0–R4 | `figures/panels/capitalized/code/render_supp2_supp3_caps.py` | `suppfig2_R{0..4}.png` |

## Supplementary Figure 3 — overlap, correlation, and discordant clonotypes

| panel | script | output |
|---|---|---|
| A | `figures/panels/capitalized/code/render_supp2_supp3_caps.py` | `suppfig3A_morisita_horn.png` |
| B | 〃 | `suppfig3B_cluster_confusion.png` |
| C | 〃 | `suppfig3C_spearman_rho.png` |
| D | 〃 | `suppfig3D_loglog_slope.png` |
| E, F | `analysis/revision_2026-08/r2_03_discordant.py` → `render_revision_panels.py` | `out/figures/suppfig3{E,F}_*.png` |

Panel A's title is left lowercase (`ss/ds overlap`) because capitalising it would give
`Ss/ds overlap`; the same applies to the `ssDNA`/`dsDNA` panel titles in Figures 1B and
Supplementary 1B, which are unchanged from the published versions.

## Supplementary Figure 4 — log-ratio at every round

| panel | script | output |
|---|---|---|
| A–F | `analysis/fig3_mapping/supp5_logratio_rounds/make_supp5_logratio.py` | `suppfig4_logratio_R{0..4}.png`, `suppfig4_logratio_summary.png` |

## Supplementary Figure 5 — HCDR3 sequence features

| panel | script | output |
|---|---|---|
| A–F | `analysis/fig3_mapping/make_supp_sequence.py` | `analysis/fig3_mapping/supp_sequence_sized.png` |

The UMAP coordinates in panels E and F are read from `atchley_umap_coords.npz` and
`esm_umap_coords.npz`, which are included so the panel can be regenerated without
running ESM-2. The AUC values quoted in the caption (0.53 and 0.55) come from the
logistic-regression probe described in the Methods.

## Analyses that produce numbers rather than panels

| script | what it reports |
|---|---|
| `analysis/revision_2026-08/r2_02_k_selection.py` | validity indices and bootstrap stability for k = 2–8, and the key slope as a function of k |
| `analysis/revision_2026-08/r2_03_discordant.py` | the 215 cluster-discordant clonotypes, by transition class and by antigen reactivity |
| `analysis/revision_2026-08/r2_04_06_regression_and_thresholds.py` | OLS vs standardised major axis vs Deming slopes; PPM cut-off sensitivity |
| `analysis/revision_2026-08/r2_08_elisa_threshold.py` | the antigen-reactivity result under five ELISA calling rules |
| `analysis/version 3/code/check_cluster_discordance.py` | the ssDNA/dsDNA cluster-label confusion matrix |
| `analysis/version 3/code/map_yoo2020_validated_clones.py` | the validated binders of Yoo et al. 2020 on the ss-vs-ds map |
