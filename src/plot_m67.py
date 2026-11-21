#!/usr/bin/env python3
"""
Plot observed CMD of M67 (NGC_2682) and overlay PARSEC isochrones
from Cavallo et al. (2024) and Bossini et al. (2019).

Reads:
  - data/raw/m67_members_hr2023.tsv  (Hunt & Reffert 2023 members)
  - data/raw/m67_parameters_cavallo2024.tsv  (Cavallo et al. 2024 parameters)
  - data/raw/m67_parsec12s_av018.dat  (PARSEC 1.2S isochrone, Av=0.18 already included)
  - data/raw/m67_parameters_bossini2019.tsv  (Bossini et al. 2019 parameters)
  - data/raw/m67_parsec12s_bossini2019_av0115.dat  (PARSEC 1.2S isochrone, Av=0.115 already included)

Outputs:
  - results/m67_members_selected.csv
  - results/m67_cmd.png
  - results/m67_isochrone_apparent.csv (Cavallo model)
  - results/m67_isochrone_bossini2019_apparent.csv (Bossini model)
  - results/m67_cmd_comparison.png (both models)
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Paths
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, "data", "raw")
RES = os.path.join(BASE, "results")
os.makedirs(RES, exist_ok=True)

MEMBERS_FILE = os.path.join(RAW, "m67_members_hr2023.tsv")

# Cavallo et al. (2024)
CAVALLO_PARAMS_FILE = os.path.join(RAW, "m67_parameters_cavallo2024.tsv")
CAVALLO_ISO_FILE = os.path.join(RAW, "m67_parsec12s_av018.dat")

# Bossini et al. (2019)
BOSSINI_PARAMS_FILE = os.path.join(RAW, "m67_parameters_bossini2019.tsv")
BOSSINI_ISO_FILE = os.path.join(RAW, "m67_parsec12s_bossini2019_av0115.dat")

CLUSTER_CAVALLO = "NGC_2682"
CLUSTER_BOSSINI = "NGC2682"


def _error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


# ============================================================
# Common utilities
# ============================================================
def read_members_tsv(path):
    """Read Hunt & Reffert (2023) members TSV.
    Returns list of dicts with keys: Name, GaiaDR3, Prob, inrt, Gmag, BPmag, RPmag.

    Данные: data/raw/m67_members_hr2023.tsv — таблица членов M67 с фотометрией Gaia DR3.
    Проверка: checks/check_comparison.py (через checks/check.py) — сверяет все строки
      выходного CSV results/m67_members_selected.csv с допуском 1e-8, без импорта src/.
    Используется в: src/plot_m67.py::main.
    """
    print(f"\nReading member stars from {os.path.basename(path)}...")
    with open(path, "r", encoding="utf-8-sig") as f:
        raw_lines = f.readlines()

    # Filter out comment lines
    data_lines = [l for l in raw_lines if not l.lstrip().startswith("#")]

    # Find header line
    header_idx = None
    for i, line in enumerate(data_lines):
        if line.startswith("Name\t") and "GaiaDR3" in line:
            header_idx = i
            break

    if header_idx is None:
        _error(f"Could not find header in {os.path.basename(path)}")

    # Parse header
    header_line = data_lines[header_idx].rstrip("\r\n")
    header_fields = header_line.split("\t")
    if len(header_fields) != 7:
        _error(f"Header has {len(header_fields)} fields, expected 7")

    # Find data start (skip units and separator lines)
    data_start = header_idx + 1
    for i in range(header_idx + 1, len(data_lines)):
        stripped = data_lines[i].strip()
        if stripped == "":
            data_start = i + 1
            continue
        sep_fields = data_lines[i].rstrip("\r\n").split("\t")
        if len(sep_fields) == 7 and all(c.strip() and all(ch == "-" for ch in c.strip()) for c in sep_fields):
            data_start = i + 1
            break
        if len(sep_fields) == 7:
            if (sep_fields[0].strip() == "" and sep_fields[1].strip() == ""
                    and sep_fields[2].strip() == "" and sep_fields[3].strip() == ""
                    and sep_fields[4].strip() == "mag" and sep_fields[5].strip() == "mag"
                    and sep_fields[6].strip() == "mag"):
                data_start = i + 1
                continue

    # Parse data rows
    total_read = 0
    members = []

    for row_num, line in enumerate(data_lines[data_start:], start=data_start + 1):
        if not line.strip():
            continue
        line_stripped = line.rstrip("\r\n")
        cols = line_stripped.split("\t")

        if len(cols) != 7:
            _error(f"Line {row_num}: expected 7 columns, got {len(cols)}: {line_stripped!r}")

        total_read += 1
        name_s, gaia_s, prob_s, inrt_s, gmag_s, bpmag_s, rpmag_s = cols

        # Validate Name
        name = name_s.strip()
        if not name:
            _error(f"Line {row_num}: empty cluster name")

        # Validate GaiaDR3: must be decimal digits only
        gaia = gaia_s.strip()
        if not gaia or not all(c in "0123456789" for c in gaia):
            _error(f"Line {row_num}: GaiaDR3 must contain decimal digits only, got {gaia!r}")

        # Validate Prob (can be empty = missing)
        prob_str = prob_s.strip()
        if prob_str:
            try:
                prob_val = float(prob_str)
            except ValueError:
                _error(f"Line {row_num}: Prob is not a number: {prob_str!r}")
        else:
            prob_val = math.nan

        # Validate inrt
        inrt_str = inrt_s.strip()
        if inrt_str not in ("", "True", "False"):
            _error(f"Line {row_num}: unexpected inrt value {inrt_str!r}")

        # Validate magnitudes (can be empty = missing)
        def parse_mag(val_s, field_name):
            v = val_s.strip()
            if not v:
                return math.nan
            try:
                fval = float(v)
            except ValueError:
                _error(f"Line {row_num}: {field_name} is not a number: {v!r}")
            if not math.isfinite(fval):
                _error(f"Line {row_num}: {field_name} is non-finite: {v!r}")
            return fval

        gmag = parse_mag(gmag_s, "Gmag")
        bpmag = parse_mag(bpmag_s, "BPmag")
        rpmag = parse_mag(rpmag_s, "RPmag")

        members.append({
            "Name": name,
            "GaiaDR3": gaia,
            "Prob": prob_val,
            "inrt": inrt_str,
            "Gmag": gmag,
            "BPmag": bpmag,
            "RPmag": rpmag,
        })

    print(f"  Total rows read: {total_read}")
    return members


def select_members(members, cluster_name):
    """Apply selection: Name == cluster, finite Gmag/BPmag/RPmag, Gmag <= 18."""
    selected = [
        m for m in members
        if m["Name"] == cluster_name
        and math.isfinite(m["Gmag"])
        and math.isfinite(m["BPmag"])
        and math.isfinite(m["RPmag"])
        and m["Gmag"] <= 18.0
    ]
    return selected


def save_members_csv(selected, out_path):
    """Save selected members to CSV with exact header."""
    print(f"\nSaving members CSV to {out_path}...")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Name", "GaiaDR3", "Gmag", "BPmag", "RPmag", "BP_RP"])
        for m in selected:
            bp_rp = m["BPmag"] - m["RPmag"]
            writer.writerow([
                m["Name"],
                m["GaiaDR3"],
                m["Gmag"],
                m["BPmag"],
                m["RPmag"],
                bp_rp,
            ])
    print(f"  Saved {len(selected)} selected stars")


def read_params_tsv(path, required_columns, cluster_name):
    """Read parameter TSV and return dict with required columns for given cluster."""
    print(f"\nReading cluster parameters from {os.path.basename(path)}...")
    with open(path, "r", encoding="utf-8-sig") as f:
        lines = f.readlines()

    header_cols = None
    data_cols = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Cluster\t") and required_columns[0] in stripped:
            reader = csv.reader([line], delimiter="\t")
            for row in reader:
                header_cols = [c.strip() for c in row]
        if stripped.startswith(cluster_name):
            reader = csv.reader([line], delimiter="\t")
            for row in reader:
                data_cols = [c.strip() for c in row]

    if header_cols is None or data_cols is None:
        _error(f"Could not find {cluster_name} data in {os.path.basename(path)}")

    params = {}
    for h, v in zip(header_cols, data_cols):
        params[h] = v

    result = {}
    for col in required_columns[1:]:
        if col not in params:
            _error(f"Missing required column {col} in {os.path.basename(path)}")
        result[col] = float(params[col])

    print(f"  {', '.join(f'{k} = {v}' for k, v in result.items())}")
    return result


def read_isochrone_dat(path, expected_log_age, expected_mh, expected_av, label_suffix):
    """Read PARSEC isochrone .dat file.
    Returns list of dicts with model points (excluding termination point).
    """
    print(f"\nReading PARSEC isochrone ({label_suffix})...")
    isochrone = []
    with open(path, "r", encoding="utf-8-sig") as f:
        iso_lines = f.readlines()

    # Find data start (after comment header)
    iso_data_start = None
    for i, line in enumerate(iso_lines):
        if not line.startswith("#") and line.strip() != "":
            iso_data_start = i
            break

    if iso_data_start is None:
        _error(f"Could not find data in {os.path.basename(path)}")

    total_iso = 0
    excluded_termination = 0

    for i in range(iso_data_start, len(iso_lines)):
        line = iso_lines[i].strip()
        if line == "" or line.startswith("#"):
            continue

        cols = line.split()
        if len(cols) != 31:
            _error(f"Isochrone line {i+1}: expected 31 columns, got {len(cols)}")

        total_iso += 1

        try:
            mh = float(cols[1])
            log_age = float(cols[2])
            mini = float(cols[3])
            label = int(cols[9])
            log_l = float(cols[6])
            gmag_model = float(cols[28])
            g_bp_model = float(cols[29])
            g_rp_model = float(cols[30])
        except (ValueError, IndexError) as e:
            _error(f"Isochrone line {i+1}: parse error: {e}")

        # Skip termination point (label=9, logL=-9.999)
        if label == 9 and log_l == -9.999:
            excluded_termination += 1
            continue

        isochrone.append({
            "Mini": mini,
            "logAge": log_age,
            "MH": mh,
            "Gmag": gmag_model,
            "G_BPmag": g_bp_model,
            "G_RPmag": g_rp_model,
            "label": label,
            "logL": log_l,
        })

    print(f"  Total isochrone points read: {total_iso}")
    print(f"  Termination points excluded: {excluded_termination}")
    print(f"  Isochrone points kept: {len(isochrone)}")

    # Verify isochrone parameters match expected
    unique_ages = set(p["logAge"] for p in isochrone)
    unique_mh = set(p["MH"] for p in isochrone)
    print(f"  Unique logAge values: {unique_ages}")
    print(f"  Unique MH values: {unique_mh}")

    if len(unique_ages) > 1 or len(unique_mh) > 1:
        _error(f"{label_suffix}: isochrone contains multiple ages or metallicities!")

    # Verify against expected parameters
    if not math.isclose(next(iter(unique_ages)), expected_log_age, rel_tol=0, abs_tol=1e-8):
        _error(f"{label_suffix}: isochrone logAge {next(iter(unique_ages))} != expected {expected_log_age}")
    if not math.isclose(next(iter(unique_mh)), expected_mh, rel_tol=0, abs_tol=1e-8):
        _error(f"{label_suffix}: isochrone MH {next(iter(unique_mh))} != expected {expected_mh}")

    return isochrone, total_iso, excluded_termination


def compute_apparent_magnitudes(isochrone, dist_mod):
    """Compute G_app = Gmag + dist_mod, BP_RP = G_BPmag - G_RPmag.

    Данные: isochrone — список точек из src/plot_m67.py::read_isochrone_dat
      (data/raw/m67_parsec12s_av018.dat или data/raw/m67_parsec12s_bossini2019_av0115.dat);
      dist_mod — dMod50 из data/raw/m67_parameters_cavallo2024.tsv или Dist из
      data/raw/m67_parameters_bossini2019.tsv (через src/plot_m67.py::read_params_tsv).
    Проверка: checks/check_comparison.py — независимо пересчитывает G_app и BP_RP
      для обеих моделей и сравнивает с results/m67_isochrone_apparent.csv и
      results/m67_isochrone_bossini2019_apparent.csv (допуск 1e-8).
    Используется в: src/plot_m67.py::main (для Cavallo и Bossini перед сохранением
      CSV и построением графиков).
    """
    for pt in isochrone:
        pt["G_app"] = pt["Gmag"] + dist_mod
        pt["BP_RP"] = pt["G_BPmag"] - pt["G_RPmag"]


def save_isochrone_csv(isochrone, out_path):
    """Save isochrone to CSV with exact header."""
    print(f"\nSaving isochrone CSV to {out_path}...")
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Gmag", "G_BPmag", "G_RPmag", "G_app", "BP_RP",
                          "Mini", "logAge", "MH", "label"])
        for pt in isochrone:
            writer.writerow([
                pt["Gmag"], pt["G_BPmag"], pt["G_RPmag"],
                pt["G_app"], pt["BP_RP"],
                pt["Mini"], pt["logAge"], pt["MH"], pt["label"],
            ])
    print(f"  Saved {len(isochrone)} model points")


def split_isochrone_by_label(isochrone):
    """Split isochrone into continuous blocks of same label.
    Returns list of (label, list_of_points) preserving original order.
    Does NOT merge separated blocks with same label.
    """
    if not isochrone:
        return []

    blocks = []
    current_label = isochrone[0]["label"]
    current_block = []

    for pt in isochrone:
        if pt["label"] != current_label:
            # Label changed - save current block and start new
            blocks.append((current_label, current_block))
            current_label = pt["label"]
            current_block = [pt]
        else:
            current_block.append(pt)

    # Save last block
    if current_block:
        blocks.append((current_label, current_block))

    return blocks


def _setup_cmd_axes(ax, obs_g, obs_bp_rp):
    """Common CMD axes setup: invert Y, labels, grid, limits."""
    g_min, g_max = obs_g.min(), obs_g.max()
    bp_rp_min, bp_rp_max = obs_bp_rp.min(), obs_bp_rp.max()

    g_margin = (g_max - g_min) * 0.03
    bp_rp_margin = (bp_rp_max - bp_rp_min) * 0.03

    ax.invert_yaxis()
    ax.set_xlabel("G_BP − G_RP (mag)")
    ax.set_ylabel("G (mag)")
    ax.grid(True, alpha=0.3)

    ax.set_xlim(bp_rp_min - bp_rp_margin, bp_rp_max + bp_rp_margin)
    ax.set_ylim(g_max + g_margin, g_min - g_margin)


def plot_cmd_single(obs_g, obs_bp_rp, model_isochrone, out_path, title):
    """Plot single-model CMD with isochrone split by label blocks."""
    print(f"\nPlotting {title}...")

    fig, ax = plt.subplots(figsize=(8, 8))

    # Observed stars: gray points
    ax.scatter(obs_bp_rp, obs_g, c="gray", s=12, alpha=0.6,
               edgecolors="none", label="Члены M67 (G ≤ 18)", zorder=2)

    # Isochrone: plot each label block as separate line segment
    blocks = split_isochrone_by_label(model_isochrone)
    for i, (label, block) in enumerate(blocks):
        bp_rp = np.array([pt["BP_RP"] for pt in block])
        g_app = np.array([pt["G_app"] for pt in block])
        # Only first block gets legend entry
        lbl = "Cavallo et al. (2024)" if i == 0 else None
        ax.plot(bp_rp, g_app, c="red", linewidth=1.5,
                label=lbl, zorder=3)

    _setup_cmd_axes(ax, obs_g, obs_bp_rp)

    # Title: one line, via suptitle
    fig.suptitle("M67 (NGC 2682): изохрона Cavallo", y=0.98, fontsize=14)

    # Legend below title, 2 columns
    fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.94),
               ncol=2, fontsize=11, frameon=False, columnspacing=1.5)

    # Make room for title + legend
    fig.subplots_adjust(top=0.88)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Report label blocks for diagnostics
    label_seq = [lbl for lbl, _ in blocks]
    print(f"  Isochrone label blocks: {label_seq}")
    print(f"  Saved plot to {out_path}")


def plot_cmd_comparison(obs_g, obs_bp_rp,
                        cav_isochrone, bos_isochrone,
                        out_path, title):
    """Plot comparison CMD with both models split by label blocks.

    Данные: obs_g, obs_bp_rp — массивы из общей выборки (src/plot_m67.py::select_members);
      cav_isochrone, bos_isochrone — списки точек с G_app, BP_RP после
      src/plot_m67.py::compute_apparent_magnitudes.
    Проверка: checks/check_comparison.py сверяет CSV с исходными данными
      и проверяет заголовки PNG. Оси, легенду, видимость обеих моделей
      и соответствие изображения CSV проверяет человек.
    Используется в: src/plot_m67.py::main (создаёт results/m67_cmd_comparison.png).
      Параллельно src/plot_m67.py::plot_cmd_single вызывается из main для
      results/m67_cmd.png (только Cavallo).
    """
    print(f"\nPlotting {title}...")

    fig, ax = plt.subplots(figsize=(8, 8))

    # Observed stars: gray points
    ax.scatter(obs_bp_rp, obs_g, c="gray", s=12, alpha=0.6,
               edgecolors="none", label="Члены M67 (G ≤ 18)", zorder=2)

    # Cavallo (2024): red line, split by label blocks
    cav_blocks = split_isochrone_by_label(cav_isochrone)
    for i, (label, block) in enumerate(cav_blocks):
        bp_rp = np.array([pt["BP_RP"] for pt in block])
        g_app = np.array([pt["G_app"] for pt in block])
        lbl = "Cavallo et al. (2024)" if i == 0 else None
        ax.plot(bp_rp, g_app, c="red", linewidth=1.5,
                label=lbl, zorder=3)

    # Bossini (2019): blue line, split by label blocks
    bos_blocks = split_isochrone_by_label(bos_isochrone)
    for i, (label, block) in enumerate(bos_blocks):
        bp_rp = np.array([pt["BP_RP"] for pt in block])
        g_app = np.array([pt["G_app"] for pt in block])
        lbl = "Bossini et al. (2019)" if i == 0 else None
        ax.plot(bp_rp, g_app, c="blue", linewidth=1.5,
                label=lbl, zorder=3)

    _setup_cmd_axes(ax, obs_g, obs_bp_rp)

    # Title: one line, via suptitle
    fig.suptitle("M67 (NGC 2682): сравнение изохрон", y=0.98, fontsize=14)

    # Legend below title, 3 columns
    fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.94),
               ncol=3, fontsize=11, frameon=False, columnspacing=1.5)

    # Make room for title + legend
    fig.subplots_adjust(top=0.88)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Report label blocks for diagnostics
    cav_labels = [lbl for lbl, _ in cav_blocks]
    bos_labels = [lbl for lbl, _ in bos_blocks]
    print(f"  Cavallo label blocks: {cav_labels}")
    print(f"  Bossini label blocks: {bos_labels}")
    print(f"  Saved plot to {out_path}")


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    # --------------------------------------------------------
    # 1. Read member stars (common for both models)
    # --------------------------------------------------------
    members = read_members_tsv(MEMBERS_FILE)

    # --------------------------------------------------------
    # 2. Select common sample (Name == NGC_2682, Gmag <= 18)
    # --------------------------------------------------------
    selected = select_members(members, CLUSTER_CAVALLO)
    print(f"\nSelection: Name == {CLUSTER_CAVALLO}, finite photometry, Gmag <= 18")
    print(f"  Stars before cut: {len(members)}")
    print(f"  Stars after cut:  {len(selected)}")

    # --------------------------------------------------------
    # 3. Save common members CSV
    # --------------------------------------------------------
    members_csv_path = os.path.join(RES, "m67_members_selected.csv")
    save_members_csv(selected, members_csv_path)

    # --------------------------------------------------------
    # 4. CAVALLO et al. (2024) model
    # --------------------------------------------------------
    cavallo_params = read_params_tsv(
        CAVALLO_PARAMS_FILE,
        ["Cluster", "logAge50", "[Fe/H]50", "AV50", "dMod50"],
        CLUSTER_CAVALLO
    )
    logAge50 = cavallo_params["logAge50"]
    feh50 = cavallo_params["[Fe/H]50"]
    av50 = cavallo_params["AV50"]
    dmod50 = cavallo_params["dMod50"]

    cavallo_iso, cav_total_iso, cav_excl = read_isochrone_dat(
        CAVALLO_ISO_FILE, logAge50, feh50, av50, "Cavallo et al. 2024"
    )

    compute_apparent_magnitudes(cavallo_iso, dmod50)

    cav_csv_path = os.path.join(RES, "m67_isochrone_apparent.csv")
    save_isochrone_csv(cavallo_iso, cav_csv_path)

    # Prepare arrays for plotting
    cav_g_app = np.array([pt["G_app"] for pt in cavallo_iso])
    cav_bp_rp = np.array([pt["BP_RP"] for pt in cavallo_iso])

    # --------------------------------------------------------
    # 5. BOSSINI et al. (2019) model
    # --------------------------------------------------------
    bossini_params = read_params_tsv(
        BOSSINI_PARAMS_FILE,
        ["Cluster", "logA", "Dist", "AV", "Fe/H"],
        CLUSTER_BOSSINI
    )
    logA = bossini_params["logA"]
    dist = bossini_params["Dist"]
    av = bossini_params["AV"]
    feh = bossini_params["Fe/H"]

    bossini_iso, bos_total_iso, bos_excl = read_isochrone_dat(
        BOSSINI_ISO_FILE, logA, feh, av, "Bossini et al. 2019"
    )

    compute_apparent_magnitudes(bossini_iso, dist)

    bos_csv_path = os.path.join(RES, "m67_isochrone_bossini2019_apparent.csv")
    save_isochrone_csv(bossini_iso, bos_csv_path)

    # Prepare arrays for plotting
    bos_g_app = np.array([pt["G_app"] for pt in bossini_iso])
    bos_bp_rp = np.array([pt["BP_RP"] for pt in bossini_iso])

    # --------------------------------------------------------
    # 6. Plot single-model CMD (Cavallo only, as before)
    # --------------------------------------------------------
    obs_g = np.array([m["Gmag"] for m in selected])
    obs_bp_rp = np.array([m["BPmag"] - m["RPmag"] for m in selected])

    png_path = os.path.join(RES, "m67_cmd.png")
    plot_cmd_single(obs_g, obs_bp_rp, cavallo_iso, png_path,
                    "M67 (NGC 2682) — Observed CMD +\nPARSEC Isochrone (Cavallo et al. 2024)")

    # --------------------------------------------------------
    # 7. Plot comparison CMD (both models)
    # --------------------------------------------------------
    comparison_png = os.path.join(RES, "m67_cmd_comparison.png")
    plot_cmd_comparison(obs_g, obs_bp_rp,
                        cavallo_iso, bossini_iso,
                        comparison_png,
                        "M67 (NGC 2682) — Reconstruction comparison:\nCavallo et al. (2024) vs Bossini et al. (2019)")

    # --------------------------------------------------------
    # 8. Summary
    # --------------------------------------------------------
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Cluster: M67 (NGC_2682 / NGC2682)")
    print(f"\nParameters (Cavallo et al. 2024):")
    print(f"  logAge50 = {logAge50}")
    print(f"  [Fe/H]50 = {feh50}")
    print(f"  AV50     = {av50}")
    print(f"  dMod50   = {dmod50}")
    print(f"\nParameters (Bossini et al. 2019):")
    print(f"  logA     = {logA}")
    print(f"  [Fe/H]   = {feh}")
    print(f"  AV       = {av}")
    print(f"  Dist     = {dist}")
    print(f"\nObservations (common sample):")
    print(f"  Total rows read:     {len(members)}")
    print(f"  Stars selected:      {len(selected)} (Gmag <= 18)")
    print(f"\nIsochrone (Cavallo 2024):")
    print(f"  Total points read:   {cav_total_iso}")
    print(f"  Termination excluded: {cav_excl}")
    print(f"  Points saved:        {len(cavallo_iso)}")
    print(f"\nIsochrone (Bossini 2019):")
    print(f"  Total points read:   {bos_total_iso}")
    print(f"  Termination excluded: {bos_excl}")
    print(f"  Points saved:        {len(bossini_iso)}")
    print(f"\nOutput files:")
    print(f"  {members_csv_path}")
    print(f"  {cav_csv_path}")
    print(f"  {bos_csv_path}")
    print(f"  {png_path}")
    print(f"  {comparison_png}")
    print("=" * 60)