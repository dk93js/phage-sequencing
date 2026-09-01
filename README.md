# Template-choice bias in deep sequencing of antibody phage-display libraries

Analysis code for:

> Yang H, Lee HK, Ryu T, Yoo D-K\*, Chung J\*. *Quantifying template-choice biases in deep
> sequencing of antibody phage display libraries.* **Biomolecules** (2026), manuscript
> biomolecules-4527970.

The study compares two DNA templates prepared in parallel from the same c-Met bio-panning
campaign — intracellular replicative-form phagemid (**dsDNA**) and phagemid packaged inside
M13 particles (**ssDNA**) — across the input library and four rounds of selection.

Raw sequencing reads: NCBI SRA **PRJNA1507317**.

## Regenerating the figures

```bash
pip install -r requirements.txt
bash run_all_figures.sh          # or: PYTHON=python3 bash run_all_figures.sh
```

Run everything from the repository root. The script regenerates every panel in
dependency order; 42 of the 44 images in the submitted figure deck come out
byte-identical (the remaining two are arrow icons drawn in PowerPoint).
`PANELS.md` lists which script produces which panel.

## Layout

```
data/                                  inputs (below)
figures/panels/new color/optionB/code/ base panels in the final cluster palette
figures/panels/capitalized/code/       the same panels with capitalised titles, as submitted
analysis/version 3/code/               Figure 1A, the Figure 1C inset, clustering checks
analysis/fig3_mapping/                 Figure 3, Supplementary Figures 4 and 5
analysis/revision_2026-08/             analyses and panels added in response to peer review
out/                                   written by the revision scripts
```

The directory names are the ones the scripts were written against and are kept so that
they run unmodified; the only edits made for release were replacing absolute paths and
pointing the phage-ELISA reader at `data/phage_elisa_wells.xlsx`.

## Data

| file | contents |
|---|---|
| `report_q40_ppm100_list.csv` | the analysis set: 2,145 clonotype-by-library rows (1,052 ssDNA + 1,093 dsDNA unique HCDR3 amino-acid sequences), read counts and depth-normalised abundance (PPM) for R0-R4, after Phred Q >= 40 and PPM >= 100 filtering |
| `clusters_pooled.csv`, `clusters_ssDNA.csv`, `clusters_dsDNA.csv` | Ward (k = 3) trajectory clusters and the z-scored trajectories they were built from |
| `per_cluster_stats.csv`, `cluster_size_summary.csv` | per-round Spearman rho and log-log slope by cluster; cluster sizes |
| `read_depth_per_round.csv` | reads per library and round: those belonging to the analysis set, the same value as PPM of all quality-filtered reads, and the quality-filtered total |
| `phage_elisa_wells.xlsx` | phage-ELISA well-level calls for the retrieved clones (692 wells -> 568 unique HCDR3) |
| `elisa_calls_per_hcdr3.csv` | the same data collapsed to one call per HCDR3 (90 antigen-reactive / 478 non-reactive) |

## Requirements

Python >= 3.11 with `numpy pandas scipy matplotlib openpyxl scikit-learn pillow`.
`figures/panels/capitalized/code/` and the base panel scripts use Arial; on systems
without it Matplotlib substitutes a default sans-serif and the panels will differ
cosmetically from the published ones.

`analysis/fig3_mapping/make_supp_sequence.py` reads pre-computed UMAP coordinates from
the two `.npz` files in that directory, so it does not need `umap-learn`, `fair-esm` or
PyTorch. Those packages are needed only to recompute the embeddings from scratch.

## Key values reproduced by this code

- analysis set 1,052 ssDNA + 1,093 dsDNA clonotypes, 1,142 unique, 1,003 shared
- log-log slope of the enriching cluster C1: 0.995 (R2), 0.574 (R3), 0.549 (R4); C2 and C3 stay at ~1.0
- Morisita-Horn overlap 0.997 / 0.998 / 0.994 / 0.501 / 0.379 for R0-R4
- 692 ELISA wells -> 568 unique HCDR3 -> 90 antigen-reactive; 146 present in both templates (35 +, 111 -)
- median log10(ssDNA/dsDNA) at R3: +0.444 (reactive) vs -0.174 (non-reactive), P = 6.3e-8

## License

MIT (see `LICENSE`).
