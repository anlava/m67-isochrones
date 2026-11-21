# Colour residuals in M67

This repository contains the sources for the M67 paper, its five raw inputs and the numerical results used in the comparison. The two PARSEC reconstructions use the parameters of Cavallo et al. (2024) and Bossini et al. (2019). Each star contributes a pair of colour residuals at its observed Gaia G magnitude.

## Reproduce the paper

Use Python 3.12, Bash, and TeX Live with pdfLaTeX, BibTeX and English and Russian language support (the course style loads both). The build was tested with TeX Live 2026. From the repository root:

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-paper.txt
bash reproduce.sh
```

On macOS with MacTeX, put `/Library/TeX/texbin` first in your PATH if another TeX distribution is the default. After the dependencies are installed the analysis runs offline.

The script rebuilds the stellar sample and model curves, checks them against the raw files, calculates the baseline and tidal-radius subsets, runs the independent membership check and regression tests, renders Figure 1, and compiles [the paper](manuscript/paper.pdf). To compare the regenerated numerical outputs with the committed reference in a fresh clone, run:

```sh
git diff --exit-code -- results/m67_residuals_by_g.csv experiment/results/
```

The [original run log](results/week6_runs.csv) records the dry run, smoke test, synthetic sweep, timing measurements and full analysis reported in METHOD.md. Running the analysis appends a new timestamped row for the full run.

## Inputs and results

- [Membership catalogue](https://doi.org/10.1051/0004-6361/202346285) and the [Cavallo](https://doi.org/10.3847/1538-3881/ad07e5) and [Bossini](https://doi.org/10.1051/0004-6361/201834693) parameter extracts are in [data/raw](data/raw). Each TSV header records its exact VizieR query. The two CMD 3.9 model files preserve the PARSEC settings and primary references in their headers. Please cite the original data and model providers when reusing these inputs.
- [Method](METHOD.md), [analysis code](src), [independent checks and tests](checks), [tidal-radius check](experiment/check_inrt_residuals.py).
- [Table 1](experiment/results/m67_residuals_inrt_comparison.csv) has six rows: three G intervals for each sample. The expected counts are 235, 204 and 219 in the full sample, and 143, 95 and 81 within the tidal radius. [Membership IDs](experiment/results/m67_residuals_inrt_ids.csv) identify the stars used for each row.
- [Figure 1 source](manuscript/figures/plot_m67_figure1.py) reads the output CSVs. The numerical table in [the LaTeX source](manuscript/paper.tex) reports the CSV values rounded to four decimal places.
