#!/usr/bin/env python3
"""Test check.py using small, explicit fixtures and deliberately wrong copies.

This validates the checker, not the scientific model or src/plot_m67.py.
Fixtures live in a temporary directory; real data/results are never touched.
"""

import base64
import csv
import subprocess
import sys
import tempfile
from pathlib import Path


CHECK = Path(__file__).resolve().with_name("check.py")
MODEL_HEADER = (
    "Zini MH logAge Mini int_IMF Mass logL logTe logg label McoreTP C_O "
    "period0 period1 period2 period3 period4 pmode Mloss tau1m X Y Xc Xn "
    "Xo Cexcess Z mbolmag Gmag G_BPmag G_RPmag"
).split()


def write_tsv(path, header, units, rows):
    lines = [
        "# Synthetic format/arithmetic fixture, not scientific M67 data",
        "\t".join(header), "\t".join(units), "\t".join(["------"] * len(header)),
    ]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        csv.writer(stream).writerows(rows)


def make_fixture(root):
    raw = root / "data" / "raw"
    result = root / "results"
    raw.mkdir(parents=True)
    result.mkdir()
    write_tsv(
        raw / "m67_parameters_cavallo2024.tsv",
        ["Cluster", "logAge50", "[Fe/H]50", "AV50", "dMod50"],
        ["", "[yr]", "[-]", "mag", "mag"],
        [["NGC_2682", "9.40", "0.32", "0.18", "9.64"]],
    )
    write_tsv(
        raw / "m67_members_hr2023.tsv",
        ["Name", "GaiaDR3", "Prob", "inrt", "Gmag", "BPmag", "RPmag"],
        ["", "", "", "", "mag", "mag", "mag"],
        [
            ["NGC_2682", "604904641373912960", "0.9", "True", "15", "16.3", "14.7"],
            ["NGC_2682", "604904641373912961", "0.8", "False", "18", "19.2", "17.6"],
            ["NGC_2682", "604904641373912962", "0.9", "True", "18.001", "19.5", "17.7"],
            ["OTHER", "604904641373912963", "0.9", "True", "12", "13", "11"],
            ["NGC_2682", "604904641373912964", "0.9", "True", "17", "nan", "16"],
            ["NGC_2682", "604904641373912965", "0.9", "True", "-inf", "17", "16"],
            ["NGC_2682", "604904641373912966", "0.9", "True", "17", "18", ""],
        ],
    )
    lines = [
        "# Synthetic format/arithmetic fixture, not a complete physical isochrone",
        "# Attention: photometry includes extinction of Av=0.18, with coefficients derived star-by-star.",
        "# " + " ".join(MODEL_HEADER),
    ]
    for values in (
        {"Mini": "0.0991749018", "Gmag": "15.232", "G_BPmag": "18.076", "G_RPmag": "13.740"},
        {"Mini": "0.1039476618", "Gmag": "15.083", "G_BPmag": "17.923", "G_RPmag": "13.592"},
        {"Mini": "1.6751031876", "label": "9", "logL": "-9.999",
         "Gmag": "30.128", "G_BPmag": "30.141", "G_RPmag": "30.089"},
    ):
        row = dict.fromkeys(MODEL_HEADER, "0")
        row.update({"Zini": "0.029013", "MH": "0.32000", "logAge": "9.40000"})
        row.update(values)
        lines.append(" ".join(row[key] for key in MODEL_HEADER))
    lines.append("#isochrone terminated")
    (raw / "m67_parsec12s_av018.dat").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Expected numbers are explicit. They are not computed with checker functions.
    write_csv(result / "m67_isochrone_apparent.csv", [
        ["Gmag", "G_BPmag", "G_RPmag", "G_app", "BP_RP", "Mini", "logAge", "MH", "label"],
        ["15.232", "18.076", "13.740", "24.872", "4.336", "0.0991749018", "9.4", "0.32", "0"],
        ["15.083", "17.923", "13.592", "24.723", "4.331", "0.1039476618", "9.4", "0.32", "0"],
    ])
    write_csv(result / "m67_members_selected.csv", [
        ["Name", "GaiaDR3", "Gmag", "BPmag", "RPmag", "BP_RP"],
        ["NGC_2682", "604904641373912960", "15", "16.3", "14.7", "1.6"],
        ["NGC_2682", "604904641373912961", "18", "19.2", "17.6", "1.6"],
    ])
    (result / "m67_cmd.png").write_bytes(base64.b64decode(
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+j3ioAAAAASUVORK5CYII="
    ))


def change_cell(relative_path, row_index, column_index, value):
    def mutate(root):
        path = root / relative_path
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.reader(stream))
        rows[row_index][column_index] = value
        write_csv(path, rows)
    return mutate


def replace_text(relative_path, before, after, count=-1):
    def mutate(root):
        path = root / relative_path
        text = path.read_text(encoding="utf-8")
        if before not in text:
            raise RuntimeError(f"Test fixture no longer contains {before!r}")
        path.write_text(text.replace(before, after, count), encoding="utf-8")
    return mutate


def remove_model_row(root):
    path = root / "results/m67_isochrone_apparent.csv"
    path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")


def reverse_model_order(root):
    path = root / "results/m67_isochrone_apparent.csv"
    lines = path.read_text().splitlines()
    path.write_text("\n".join([lines[0], *reversed(lines[1:])]) + "\n")


def duplicate_parameters(root):
    path = root / "data/raw/m67_parameters_cavallo2024.tsv"
    text = path.read_text()
    path.write_text(text + text.splitlines()[-1] + "\n")


def main():
    model = "results/m67_isochrone_apparent.csv"
    members = "results/m67_members_selected.csv"
    raw_model = "data/raw/m67_parsec12s_av018.dat"
    raw_members = "data/raw/m67_members_hr2023.tsv"
    raw_params = "data/raw/m67_parameters_cavallo2024.tsv"
    cases = [
        ("correct fixture, including boundary G=18 and missing photometry", None, 0, "принято:"),
        ("extra extinction in G_app", change_cell(model, 1, 3, "25.052"), 1, "G_app"),
        ("incorrect model colour", change_cell(model, 1, 4, "4.0"), 1, "BP_RP"),
        ("missing model point", remove_model_row, 1, "rows, expected"),
        ("model order changed", reverse_model_order, 1, "Gmag"),
        ("one wrong age throughout the model",
         replace_text(raw_model, "9.40000", "9.50000"), 1, "logAge"),
        ("mixed model metallicities",
         replace_text(raw_model, "0.32000", "0.33000", 1), 1, "multiple metallicities"),
        ("wrong Av in metadata",
         replace_text(raw_model, "Av=0.18", "Av=0.28"), 1, "Av versus AV50"),
        ("malformed model row",
         replace_text(raw_model, "15.232 18.076 13.740", "15.232 18.076"), 1, "fields, expected"),
        ("rounded Gaia source ID",
         change_cell(members, 2, 1, "604904641373912960"), 1, "GaiaDR3"),
        ("wrong cluster in selected output",
         change_cell(members, 1, 0, "OTHER"), 1, ":Name"),
        ("non-finite selected photometry",
         change_cell(members, 1, 3, "nan"), 1, "non-finite"),
        ("malformed member row",
         replace_text(raw_members, "\t15\t16.3\t14.7", "\t15\t16.3"), 1, "fields, expected"),
        ("invalid photometry token",
         replace_text(raw_members, "\t15\t16.3\t14.7", "\tbad\t16.3\t14.7"), 1, "not a number"),
        ("malformed parameter row",
         replace_text(raw_params, "\t0.18\t9.64", "\t0.18"), 1, "fields, expected"),
        ("ambiguous parameter row", duplicate_parameters, 1, "exactly one NGC_2682"),
    ]
    failures = []
    with tempfile.TemporaryDirectory(prefix="m67-check-self-test-") as temporary:
        for index, (name, mutate, expected_code, expected_text) in enumerate(cases):
            root = Path(temporary) / str(index)
            make_fixture(root)
            if mutate:
                mutate(root)
            run = subprocess.run(
                [sys.executable, str(CHECK), "--project-dir", str(root)],
                capture_output=True, text=True, timeout=15,
            )
            output = run.stdout + run.stderr
            if run.returncode != expected_code or expected_text not in output:
                failures.append(name)
                print(f"FAIL: {name}: exit={run.returncode}\n{output}", file=sys.stderr)
            else:
                print(f"OK: {name}: checker exit={run.returncode}")
    if failures:
        print(f"FAILED: {len(failures)} of {len(cases)} checker tests", file=sys.stderr)
        return 1
    print(f"PASSED: {len(cases)} checker tests; correct result -> 0, damaged copies -> 1.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
