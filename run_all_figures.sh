#!/bin/bash
# Regenerate every published panel, in dependency order.
# Run from the repository root:  bash run_all_figures.sh
set -e
PY="${PYTHON:-python}"

echo "== base panels (final palette) =="
$PY "figures/panels/new color/optionB/code/make_panels.py"
$PY "figures/panels/new color/optionB/code/make_panels_sized.py"
$PY "figures/panels/new color/optionB/code/make_figures_sized.py"
$PY "figures/panels/new color/optionB/code/render_recolored.py"
$PY "figures/panels/new color/optionB/code/render_fig1B_promoted.py"

echo "== capitalised panel titles (the versions in the submitted figures) =="
for s in figures/panels/capitalized/code/render_*.py; do $PY "$s"; done

echo "== Figure 1A =="
mkdir -p "analysis/version 3/panels"
$PY "analysis/version 3/code/panelA_prototypes.py"

echo "== Figure 3 and the sequence-feature supplementary figure =="
$PY analysis/fig3_mapping/make_fig3_sized.py
$PY analysis/fig3_mapping/make_logratio_panel.py
$PY analysis/fig3_mapping/make_supp_sequence.py

echo "== Supplementary Figure 4 (log-ratio at every round) =="
$PY analysis/fig3_mapping/supp5_logratio_rounds/make_supp5_logratio.py

echo "== analyses and panels added at revision =="
for s in analysis/revision_2026-08/r2_0*.py; do $PY "$s"; done
$PY analysis/revision_2026-08/render_revision_panels.py

echo "done"
