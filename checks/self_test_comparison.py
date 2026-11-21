#!/usr/bin/env python3
"""Test the additional checker on explicit fixtures and deliberate errors.

Install beside check.py, self_test.py and check_comparison.py. All fixtures
are temporary; this does not modify or test the real production program.
"""

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
import self_test as fixtures


CHECK = Path(__file__).with_name("check_comparison.py")
MODEL = "results/m67_isochrone_bossini2019_apparent.csv"
RAW_MODEL = "data/raw/m67_parsec12s_bossini2019_av0115.dat"
PARAMS = "data/raw/m67_parameters_bossini2019.tsv"
COMMON_METADATA = (
    "# isochrones based on PARSEC release v1.2S + COLIBRI S_37 + S_35 + PR16\n"
    "# Photometric system: Gaia EDR3 (all Vegamags, Gaia passbands from ESA/Gaia website)\n"
    "# Using YBC version of bolometric corrections as in Chen et al. (2019)\n"
    "# O-rich circumstellar dust ignored\n"
    "# C-rich circumstellar dust ignored\n"
)


def make_fixture(root):
    fixtures.make_fixture(root)
    raw = root / "data" / "raw"
    results = root / "results"
    old_model = raw / "m67_parsec12s_av018.dat"
    old_text = old_model.read_text(encoding="utf-8").replace(
        "with coefficients derived star-by-star.",
        "with coefficients derived star-by-star, for Cardelli et al 89 + O'Donnell 94 Rv=3.1 extinction curve.")
    old_model.write_text(COMMON_METADATA + old_text, encoding="utf-8")
    fixtures.write_tsv(raw / Path(PARAMS).name,
                       ["Cluster", "logA", "Dist", "AV", "AGTO", "Fe/H"],
                       ["", "[yr]", "mag", "mag", "mag", ""],
                       [["NGC2682", "9.561", "9.726", "0.115", "0.112", "0.000"]])
    lines = [
        COMMON_METADATA.rstrip(),
        "# Synthetic test fixture: not a scientific isochrone",
        "# Attention: photometry includes extinction of Av=0.115, with coefficients derived star-by-star, "
        "for Cardelli et al 89 + O'Donnell 94 Rv=3.1 extinction curve.",
        "# " + " ".join(fixtures.MODEL_HEADER),
    ]
    for values in (
        {"Mini": "0.0900000036", "Gmag": "15.436", "G_BPmag": "18.409", "G_RPmag": "13.948"},
        {"Mini": "0.1000000015", "Gmag": "14.975", "G_BPmag": "17.935", "G_RPmag": "13.489"},
        # Explicit edge case: label=9 ALONE must not cause removal.
        {"Mini": "1.2", "label": "9", "logL": "-1", "Gmag": "2", "G_BPmag": "3", "G_RPmag": "1.5"},
        {"Mini": "1.3960039616", "label": "9", "logL": "-9.999",
         "Gmag": "30.087", "G_BPmag": "30.091", "G_RPmag": "30.066"},
    ):
        row = dict.fromkeys(fixtures.MODEL_HEADER, "0")
        row.update({"Zini": "0.01471", "MH": "0.00000", "logAge": "9.56100", "label": "1"})
        row.update(values)
        lines.append(" ".join(row[key] for key in fixtures.MODEL_HEADER))
    lines.append("#isochrone terminated")
    (root / RAW_MODEL).write_text("\n".join(lines) + "\n", encoding="utf-8")
    # Expected arithmetic is explicit, without using checker/production functions.
    fixtures.write_csv(root / MODEL, [
        ["Gmag", "G_BPmag", "G_RPmag", "G_app", "BP_RP", "Mini", "logAge", "MH", "label"],
        ["15.436", "18.409", "13.948", "25.162", "4.461", "0.0900000036", "9.561", "0", "1"],
        ["14.975", "17.935", "13.489", "24.701", "4.446", "0.1000000015", "9.561", "0", "1"],
        ["2", "3", "1.5", "11.726", "1.5", "1.2", "9.561", "0", "9"],
    ])
    (results / "m67_cmd_comparison.png").write_bytes((results / "m67_cmd.png").read_bytes())


def edit_rows(transform):
    def mutate(root):
        path = root / MODEL
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.reader(stream))
        fixtures.write_csv(path, transform(rows))
    return mutate


def duplicate_parameters(root):
    path = root / PARAMS
    text = path.read_text(encoding="utf-8")
    path.write_text(text + text.splitlines()[-1] + "\n", encoding="utf-8")


def swap_models(root):
    (root / MODEL).write_bytes((root / "results/m67_isochrone_apparent.csv").read_bytes())


def main():
    cases = [
        ("correct fixture, including faint points and non-terminal label=9", None, 0, "Bossini"),
        ("Cavallo distance modulus applied to Bossini",
         fixtures.change_cell(MODEL, 1, 3, "25.076"), 1, "G_app"),
        ("Av added a second time", fixtures.change_cell(MODEL, 1, 3, "25.277"), 1, "G_app"),
        ("wrong model colour", fixtures.change_cell(MODEL, 1, 4, "-4.461"), 1, "BP_RP"),
        ("Cavallo model saved under Bossini name", swap_models, 1, "rows, expected"),
        ("all label=9 points removed", edit_rows(lambda rows: rows[:3]), 1, "rows, expected"),
        ("terminal point retained", edit_rows(lambda rows: rows + [[
            "30.087", "30.091", "30.066", "39.813", "0.025", "1.3960039616", "9.561", "0", "9"
        ]]), 1, "rows, expected"),
        ("source model order changed", edit_rows(lambda rows: [rows[0], *reversed(rows[1:])]), 1, "Gmag"),
        ("source photometry changed", fixtures.change_cell(MODEL, 1, 1, "18.4"), 1, "G_BPmag"),
        ("AGTO mistaken for Av", fixtures.replace_text(RAW_MODEL, "Av=0.115", "Av=0.112"), 1, "Av versus AV"),
        ("wrong age throughout source model", fixtures.replace_text(RAW_MODEL, "9.56100", "9.40000"), 1, "logAge"),
        ("mixed model metallicities", fixtures.replace_text(RAW_MODEL, "0.00000", "0.32000", 1), 1, "multiple metallicities"),
        ("damaged model row", fixtures.replace_text(RAW_MODEL, "15.436 18.409 13.948", "15.436 18.409"), 1, "fields, expected"),
        ("duplicate parameter row", duplicate_parameters, 1, "expected exactly one NGC2682"),
        ("wrong cluster spelling", fixtures.replace_text(PARAMS, "NGC2682", "NGC_2682"), 1, "expected exactly one NGC2682"),
        ("malformed parameter row", fixtures.replace_text(PARAMS, "9.726\t0.115", "9.726"), 1, "fields, expected"),
        ("non-finite model output", fixtures.change_cell(MODEL, 1, 0, "nan"), 1, "non-finite"),
        ("Gaia DR2 used for second model", fixtures.replace_text(RAW_MODEL, "Gaia EDR3", "Gaia DR2"), 1, "Gaia EDR3"),
        ("Cavallo settings changed", fixtures.replace_text("data/raw/m67_parsec12s_av018.dat", "Gaia EDR3", "Gaia DR2"), 1, "Gaia EDR3"),
        ("previous result regressed", fixtures.change_cell("results/m67_isochrone_apparent.csv", 1, 3, "0"), 1, "G_app"),
        ("shared stellar sample changed", fixtures.change_cell("results/m67_members_selected.csv", 1, 1, "604904641373912999"), 1, "GaiaDR3"),
        ("comparison PNG absent", lambda root: (root / "results/m67_cmd_comparison.png").unlink(), 1, "m67_cmd_comparison.png"),
    ]
    failures = []
    with tempfile.TemporaryDirectory(prefix="m67-comparison-check-") as temporary:
        for index, (name, mutate, expected_code, expected_text) in enumerate(cases):
            root = Path(temporary) / str(index)
            make_fixture(root)
            if mutate is not None:
                mutate(root)
            run = subprocess.run([sys.executable, str(CHECK), "--project-dir", str(root)],
                                 capture_output=True, text=True, timeout=15)
            output = run.stdout + run.stderr
            if run.returncode != expected_code or expected_text not in output:
                failures.append(name)
                print(f"FAIL: {name}: exit={run.returncode}\n{output}", file=sys.stderr)
            else:
                print(f"OK: {name}: checker exit={run.returncode}")
    if failures:
        print(f"FAILED: {len(failures)} of {len(cases)} comparison checker tests", file=sys.stderr)
        return 1
    print(f"PASSED: {len(cases)} comparison checker tests; correct result -> 0, damaged copies -> 1.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
