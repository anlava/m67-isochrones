#!/usr/bin/env python3
"""Independent check for inrt=True residuals.

Reads raw TSV directly, filters for base and inrt=True, and compares GaiaDR3 sets
with experiment/results/m67_residuals_inrt_ids.csv.

Checks:
  - Each GaiaDR3 appears exactly once per (sample, g_bin) pair in the ID CSV.
  - Presence in both base and inrt_true is expected (overlap allowed).
  - n in each of the 6 rows of the residuals table matches the ID count.
"""

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_TSV = ROOT / "data" / "raw" / "m67_members_hr2023.tsv"
IDS_CSV = ROOT / "experiment" / "results" / "m67_residuals_inrt_ids.csv"
TABLE_CSV = ROOT / "experiment" / "results" / "m67_residuals_inrt_comparison.csv"
BINS = ("[15,16)", "[16,17)", "[17,18]")


def bin_index(g):
    if 15 <= g < 16:
        return 0
    if 16 <= g < 17:
        return 1
    if 17 <= g <= 18:
        return 2
    return None


def read_raw_tsv(path):
    """Read raw TSV and return list of (GaiaDR3, Gmag, inrt) for NGC_2682."""
    records = []
    with path.open(encoding="utf-8-sig", newline="") as f:
        lines = f.readlines()

    # Find header
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("Name\t") and "GaiaDR3" in line:
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("Could not find header in raw TSV")

    # Find data start
    data_start = header_idx + 1
    for i in range(header_idx + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped == "":
            data_start = i + 1
            continue
        sep_fields = lines[i].rstrip("\r\n").split("\t")
        if len(sep_fields) == 7 and all(
            c.strip() and all(ch == "-" for ch in c.strip()) for c in sep_fields
        ):
            data_start = i + 1
            break
        if len(sep_fields) == 7:
            if (
                sep_fields[0].strip() == ""
                and sep_fields[1].strip() == ""
                and sep_fields[2].strip() == ""
                and sep_fields[3].strip() == ""
                and sep_fields[4].strip() == "mag"
                and sep_fields[5].strip() == "mag"
                and sep_fields[6].strip() == "mag"
            ):
                data_start = i + 1
                continue

    # Parse data
    for row_num, line in enumerate(lines[data_start:], start=data_start + 1):
        if not line.strip():
            continue
        cols = line.rstrip("\r\n").split("\t")
        if len(cols) != 7:
            continue

        name, gaia, prob, inrt, gmag, bpmag, rpmag = cols
        name = name.strip()
        gaia = gaia.strip()
        inrt = inrt.strip()
        gmag_str = gmag.strip()
        bpmag_str = bpmag.strip()
        rpmag_str = rpmag.strip()

        if name != "NGC_2682":
            continue

        # Parse magnitudes
        try:
            gmag_val = float(gmag_str) if gmag_str else float("nan")
            bpmag_val = float(bpmag_str) if bpmag_str else float("nan")
            rpmag_val = float(rpmag_str) if rpmag_str else float("nan")
        except ValueError:
            continue

        # Filter: finite photometry, Gmag <= 18
        if not (
            math.isfinite(gmag_val)
            and math.isfinite(bpmag_val)
            and math.isfinite(rpmag_val)
            and gmag_val <= 18.0
        ):
            continue

        records.append((gaia, gmag_val, inrt))

    return records


def main():
    print(f"ROOT: {ROOT}")
    print(f"RAW_TSV: {RAW_TSV}")
    print(f"IDS_CSV: {IDS_CSV}")
    print(f"TABLE_CSV: {TABLE_CSV}")

    # 1. Read raw TSV
    print("\n1. Чтение raw TSV...")
    records = read_raw_tsv(RAW_TSV)
    print(f"  Всего записей (конечные, G<=18): {len(records)}")

    # 2. Group by G bin for base (all) and inrt_true (inrt=="True")
    print("\n2. Группировка по интервалам G...")
    base_groups = {i: [] for i in range(3)}
    inrt_groups = {i: [] for i in range(3)}
    for gaia, gmag, inrt in records:
        i = bin_index(gmag)
        if i is None:
            continue
        base_groups[i].append(gaia)
        if inrt == "True":
            inrt_groups[i].append(gaia)

    for i, g_bin in enumerate(BINS):
        print(f"  {g_bin}: base={len(base_groups[i])}, inrt={len(inrt_groups[i])}")

    # 3. Read IDs CSV
    print("\n3. Чтение IDs CSV...")
    csv_base = {i: [] for i in range(3)}
    csv_inrt = {i: [] for i in range(3)}
    with IDS_CSV.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["sample"]
            i = BINS.index(row["g_bin"])
            gaia = row["GaiaDR3"]
            if sample == "base":
                csv_base[i].append(gaia)
            elif sample == "inrt_true":
                csv_inrt[i].append(gaia)

    # 4. Check multiplicities in IDs CSV
    print("\n4. Проверка уникальности ID в каждой паре (sample, g_bin)...")
    errors = []
    for i, g_bin in enumerate(BINS):
        for label, ids in [("base", csv_base[i]), ("inrt_true", csv_inrt[i])]:
            seen = {}
            for gaia in ids:
                seen[gaia] = seen.get(gaia, 0) + 1
            for gaia, cnt in seen.items():
                if cnt > 1:
                    errors.append(
                        f"  ОШИБКА: {g_bin}.{label}: GaiaDR3 {gaia} встречается {cnt} раз(а)"
                    )

    if errors:
        for e in errors:
            print(e)
        print("\nПроверка не пройдена: дубликаты ID.")
        return 1

    # 5. Compare sets (order-independent)
    print("\n5. Сравнение наборов ID...")
    all_ok = True
    for i, g_bin in enumerate(BINS):
        raw_base = set(base_groups[i])
        raw_inrt = set(inrt_groups[i])
        csv_base_set = set(csv_base[i])
        csv_inrt_set = set(csv_inrt[i])

        for label, raw_set, csv_set in [
            ("base", raw_base, csv_base_set),
            ("inrt_true", raw_inrt, csv_inrt_set),
        ]:
            if raw_set == csv_set:
                print(f"  {g_bin}.{label}: СОВПАДАЕТ ({len(raw_set)} звёзд)")
            else:
                print(f"  {g_bin}.{label}: НЕ СОВПАДАЕТ!")
                print(f"    В raw, но не в CSV: {raw_set - csv_set}")
                print(f"    В CSV, но не в raw: {csv_set - raw_set}")
                all_ok = False

    # 6. Check n in comparison table
    print("\n6. Проверка n в итоговой таблице...")
    table_n = {}
    with TABLE_CSV.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row["sample"], row["g_bin"])
            table_n[key] = int(row["n"])

    for i, g_bin in enumerate(BINS):
        for label, ids in [("base", csv_base[i]), ("inrt_true", csv_inrt[i])]:
            key = (label, g_bin)
            expected_n = len(ids)
            actual_n = table_n.get(key)
            if actual_n is None:
                print(f"  ОШИБКА: {key}: строка отсутствует в таблице")
                all_ok = False
            elif actual_n != expected_n:
                print(
                    f"  ОШИБКА: {key}: ожидаемое n={expected_n}, фактическое n={actual_n}"
                )
                all_ok = False
            else:
                print(f"  {key}: n={actual_n} OK")

    if all_ok:
        print("\nПроверка пройдена: наборы ID и n совпадают.")
        return 0
    else:
        print("\nПроверка не пройдена: расхождения.")
        return 1


if __name__ == "__main__":
    sys.exit(main())