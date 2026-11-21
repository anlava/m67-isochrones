#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export MPLBACKEND=Agg

python src/plot_m67.py
python checks/check_comparison.py
python src/m67_residuals.py full
python src/m67_residuals_inrt.py
python experiment/check_inrt_residuals.py
python -m pytest -q checks
python manuscript/figures/plot_m67_figure1.py

cd manuscript
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
bibtex paper
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
pdflatex -interaction=nonstopmode -halt-on-error paper.tex
