# Template-choice bias in deep sequencing of antibody phage-display libraries

Analysis code for:

> Yang H, Lee HK, Ryu T, Yoo D-K\*, Chung J\*. *Quantifying template-choice biases in deep
> sequencing of antibody phage display libraries.* **Biomolecules** (2026), manuscript
> biomolecules-4527970.

The study compares two DNA templates prepared in parallel from the same c-Met bio-panning
campaign — intracellular replicative-form phagemid (**dsDNA**) and phagemid packaged inside
M13 particles (**ssDNA**) — across the input library and four rounds of selection.

Raw sequencing reads: NCBI SRA **PRJNA1507317**.

## Layout

```
data/       inputs (see "Data" below)
figures/    scripts that produce the published main and supplementary figure panels
revision/   analyses added in response to peer review
out/        everything the scripts write (created on first run)
```

## Requirements

Python >= 3.11 with `numpy pandas scipy matplotlib openpyxl scikit-learn`
(`umap-learn fair-esm torch` are needed only for `figures/make_supp_sequence.py`).

```bash
pip install -r requirements.txt
```

Run every script from the repository root, e.g. `python revision/r2_02_k_selection.py`.

## Data

| file | contents |
|---|---|
| `report_q40_ppm100_list.csv` | the analysis set: 2,145 clonotype-by-library rows (1,052 ssDNA + 1,093 dsDNA unique HCDR3 amino-acid sequences), read counts and depth-normalised abundance (PPM) for R0-R4, after Phred Q >= 40 and PPM >= 100 filtering |
| `clusters_pooled.csv`, `clusters_ssDNA.csv`, `clusters_dsDNA.csv` | Ward (k = 3) trajectory clusters and the z-scored trajectories they were built from |
| `per_cluster_stats.csv`, `cluster_size_summary.csv` | per-round Spearman rho and log-log slope by cluster; cluster sizes |
| `read_depth_per_round.csv` | reads per library and round before and after filtering |
| `phage_elisa_wells.xlsx` | phage-ELISA well-level calls for the retrieved clones (692 wells -> 568 unique HCDR3) |
| `elisa_calls_per_hcdr3.csv` | the same data collapsed to one call per HCDR3 (90 antigen-reactive / 478 non-reactive) |

## What the scripts do

**figures/**

- `make_panels.py`, `make_panels_sized.py`, `make_figures_sized.py`, `panelA_prototypes.py`,
  `panel_C_prototypes.py`, `make_fig1C_inset.py` — Figure 1 and 2 panels: diversity, Top-1%
  concentration, Ward clustering, per-round composition, ss-vs-ds scatters, HCDR3 length.
- `check_cluster_discordance.py` — the 3 x 3 ssDNA/dsDNA cluster-label confusion matrix.
- `make_fig3_sized.py`, `make_logratio_panel.py` — Figure 3: antigen-reactive clonotypes on the
  ss-vs-ds scatter and the log10(ssDNA/dsDNA) comparison.
- `make_supp5_logratio.py` — the same log-ratio comparison at every round (R0-R4).
- `make_supp_sequence.py` — HCDR3 sequence features, Atchley-factor and ESM-2 embeddings,
  logistic-regression probe and UMAP projections.
- `map_yoo2020_validated_clones.py` — maps the validated binders of Yoo et al. 2020 onto the map.

**revision/**

- `r2_02_k_selection.py` — silhouette / Calinski-Harabasz / Davies-Bouldin / gap statistic for
  k = 2-8, bootstrap cluster stability (Jaccard, 200 resamples) and the key slope as a function of k.
- `r2_03_discordant.py` — the 215 cluster-discordant clonotypes: transition classes, their
  log-ratios, and the antigen-reactivity of clonotypes the two templates disagree about.
- `r2_04_06_regression_and_thresholds.py` — standardised major axis and Deming regression
  alongside OLS; sensitivity of the whole pipeline to the PPM cut-off (100 / 200 / 500 / 1000).
- `r2_08_elisa_threshold.py` — sensitivity of the antigen-reactivity result to the ELISA
  calling rule.
- `render_revision_panels.py` — panels for the supplementary figures added at revision.

## Key values reproduced by this code

- analysis set 1,052 ssDNA + 1,093 dsDNA clonotypes, 1,142 unique, 1,003 shared
- log-log slope of the enriching cluster C1: 0.995 (R2), 0.574 (R3), 0.549 (R4); C2 and C3 stay at ~1.0
- Morisita-Horn overlap 0.997 / 0.998 / 0.994 / 0.501 / 0.379 for R0-R4
- 692 ELISA wells -> 568 unique HCDR3 -> 90 antigen-reactive; 146 present in both templates (35 +, 111 -)
- median log10(ssDNA/dsDNA) at R3: +0.444 (reactive) vs -0.174 (non-reactive), P = 6.3e-8

## License

MIT (see `LICENSE`).
