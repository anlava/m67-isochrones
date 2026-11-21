#!/usr/bin/env python3
"""Calculate colour residuals for inrt=True subset of M67 members.

Reuses validated_inputs, summarise_bin, result_row from m67_residuals.py.
Reads raw TSV to map GaiaDR3 -> inrt.
Outputs:
  - experiment/results/m67_residuals_inrt_comparison.csv (6 rows)
  - experiment/results/m67_residuals_inrt_ids.csv (sample, g_bin, GaiaDR3)
"""

import csv
import math
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent))

from m67_residuals import (
    validated_inputs,
    summarise_bin,
    result_row,
    SCIENCE_FIELDS,
    BINS,
    MIN_N,
)

ROOT = Path(__file__).resolve().parents[1]
RAW_TSV = ROOT / "data" / "raw" / "m67_members_hr2023.tsv"
OUTPUT_DIR = ROOT / "experiment" / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def read_gaia_to_inrt(path):
    """Read raw TSV and return GaiaDR3 -> inrt mapping.
    Returns dict with string keys and string values ('True', 'False', '').
    """
    mapping = {}
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
            continue  # Skip malformed lines

        name, gaia, prob, inrt, gmag, bpmag, rpmag = cols
        name = name.strip()
        gaia = gaia.strip()
        inrt = inrt.strip()

        # Only store for NGC_2682
        if name == "NGC_2682" and gaia:
            mapping[gaia] = inrt

    return mapping


def filter_groups_by_inrt(groups, gaia_to_inrt):
    """Filter each group to keep only stars with inrt == 'True'.
    Preserves original order within each group.
    """
    filtered = []
    for group in groups:
        filtered_group = [
            (gaia, g, colour)
            for gaia, g, colour in group
            if gaia_to_inrt.get(gaia) == "True"
        ]
        filtered.append(filtered_group)
    return filtered


def verify_base_against_csv(rows_base, output_path):
    """Verify base results against existing CSV.

    Compares all SCIENCE_FIELDS for each g_bin with absolute tolerance 1e-10.
    n and status are compared exactly (string).
    For numeric fields: checks math.isfinite of both values, then compares.
    CSV values are "expected", rows_base values are "actual".
    Exits with code 1 on any mismatch before writing output CSVs.
    """
    existing = {}
    with output_path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gb = row["g_bin"]
            if gb in existing:
                print(f"  ОШИБКА: дубликат g_bin '{gb}' в CSV")
                sys.exit(1)
            existing[gb] = row

    # Exact match of bin sets
    if set(existing.keys()) != set(BINS):
        extra = set(existing.keys()) - set(BINS)
        missing = set(BINS) - set(existing.keys())
        if extra:
            print(f"  ОШИБКА: лишние интервалы в CSV: {extra}")
        if missing:
            print(f"  ОШИБКА: отсутствуют интервалы в CSV: {missing}")
        sys.exit(1)

    print("\nВерификация базовых результатов:")
    for row in rows_base:
        g_bin = row["g_bin"]
        csv_row = existing[g_bin]
        for field in SCIENCE_FIELDS:
            if field == "g_bin":
                continue
            csv_val = csv_row[field]
            row_val = str(row[field]) if row[field] is not None else ""
            if field in ("n", "status"):
                if csv_val != row_val:
                    print(
                        f"  ОШИБКА: {g_bin}.{field}: ожидаемое={csv_val!r}, фактическое={row_val!r}"
                    )
                    sys.exit(1)
            else:
                # Numeric field
                if csv_val == "" and row_val == "":
                    continue  # both empty (n < MIN_N case)
                try:
                    exp_val = float(csv_val)
                except ValueError:
                    print(
                        f"  ОШИБКА: {g_bin}.{field}: ожидаемое не число '{csv_val}'"
                    )
                    sys.exit(1)
                try:
                    act_val = float(row_val)
                except ValueError:
                    print(
                        f"  ОШИБКА: {g_bin}.{field}: фактическое не число '{row_val}'"
                    )
                    sys.exit(1)
                if not math.isfinite(exp_val):
                    print(
                        f"  ОШИБКА: {g_bin}.{field}: ожидаемое не конечно ({exp_val})"
                    )
                    sys.exit(1)
                if not math.isfinite(act_val):
                    print(
                        f"  ОШИБКА: {g_bin}.{field}: фактическое не конечно ({act_val})"
                    )
                    sys.exit(1)
                if abs(exp_val - act_val) > 1e-10:
                    print(
                        f"  ОШИБКА: {g_bin}.{field}: ожидаемое={exp_val}, фактическое={act_val}"
                    )
                    sys.exit(1)
        print(f"  {g_bin}: OK")
    print("  Все поля совпадают.")


def main():
    print(f"Python: {sys.version}")
    print(f"ROOT: {ROOT}")

    # 1. Get base groups and curves
    print("\n1. Вызов validated_inputs...")
    groups_base, curves = validated_inputs(ROOT)
    print(f"  Базовые группы: {[len(g) for g in groups_base]}")

    # 2. Read GaiaDR3 -> inrt mapping
    print("\n2. Чтение mapping GaiaDR3 -> inrt...")
    gaia_to_inrt = read_gaia_to_inrt(RAW_TSV)
    print(f"  Всего записей в mapping: {len(gaia_to_inrt)}")

    # Check uniqueness and coverage
    unique_gaia = set(gaia_to_inrt.keys())
    print(f"  Уникальных GaiaDR3: {len(unique_gaia)}")

    # Check for empty inrt values
    empty_inrt = sum(1 for v in gaia_to_inrt.values() if v == "")
    print(f"  Пустые inrt: {empty_inrt}")

    # 3. Filter for inrt=True
    print("\n3. Фильтрация для inrt=True...")
    groups_inrt = filter_groups_by_inrt(groups_base, gaia_to_inrt)
    print(f"  inrt=True группы: {[len(g) for g in groups_inrt]}")

    # 4. Calculate residuals for base
    print("\n4. Расчёт остатков для base...")
    summaries_base = [summarise_bin(group, curves) for group in groups_base]
    rows_base = [result_row(g_bin, s) for g_bin, s in zip(BINS, summaries_base)]

    # 5. Calculate residuals for inrt=True
    print("\n5. Расчёт остатков для inrt=True...")
    summaries_inrt = [summarise_bin(group, curves) for group in groups_inrt]
    rows_inrt = [result_row(g_bin, s) for g_bin, s in zip(BINS, summaries_inrt)]

    # 6. Verify base against existing CSV
    print("\n6. Верификация base против results/m67_residuals_by_g.csv...")
    verify_base_against_csv(rows_base, ROOT / "results" / "m67_residuals_by_g.csv")

    # 7. Save comparison CSV
    print("\n7. Сохранение experiment/results/m67_residuals_inrt_comparison.csv...")
    comparison_path = OUTPUT_DIR / "m67_residuals_inrt_comparison.csv"
    with comparison_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, ["sample"] + list(SCIENCE_FIELDS))
        writer.writeheader()
        for row in rows_base:
            writer.writerow({"sample": "base", **row})
        for row in rows_inrt:
            writer.writerow({"sample": "inrt_true", **row})
    print(f"  Сохранено {len(rows_base) + len(rows_inrt)} строк")

    # 8. Save IDs CSV
    print("\n8. Сохранение experiment/results/m67_residuals_inrt_ids.csv...")
    ids_path = OUTPUT_DIR / "m67_residuals_inrt_ids.csv"
    with ids_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "g_bin", "GaiaDR3"])
        for i, g_bin in enumerate(BINS):
            for sample, groups in [("base", groups_base), ("inrt_true", groups_inrt)]:
                for gaia, g, colour in groups[i]:
                    writer.writerow([sample, g_bin, gaia])
    print(f"  Сохранено ID для {sum(len(g) for g in groups_base) + sum(len(g) for g in groups_inrt)} звёзд")

    # 9. Print table
    print("\n" + "=" * 100)
    print("ТАБЛИЦА ОСТАТКОВ")
    print("=" * 100)
    print(f"{'sample':<12} {'g_bin':<10} {'n':>4} {'median_r_cav':>14} {'median_r_bos':>14} {'D':>10} {'ci_low':>10} {'ci_high':>10} {'status'}")
    print("-" * 100)
    for sample, rows in [("base", rows_base), ("inrt_true", rows_inrt)]:
        for row in rows:
            n = row["n"]
            if n < MIN_N:
                print(f"{sample:<12} {row['g_bin']:<10} {n:>4}  {'недостаточно данных':<50}")
            else:
                print(
                    f"{sample:<12} {row['g_bin']:<10} {n:>4} "
                    f"{row['median_r_cav']:>14.8f} {row['median_r_bos']:>14.8f} "
                    f"{row['D']:>10.8f} {row['ci_low']:>10.8f} {row['ci_high']:>10.8f} "
                    f"{row['status']}"
                )
    print("=" * 100)

    # 10. Compare D sign and CI zero-crossing
    print("\nСРАВНЕНИЕ D:")
    for i, g_bin in enumerate(BINS):
        base_row = rows_base[i]
        inrt_row = rows_inrt[i]
        base_d = base_row.get("D")
        inrt_d = inrt_row.get("D")
        base_ci = (base_row.get("ci_low"), base_row.get("ci_high"))
        inrt_ci = (inrt_row.get("ci_low"), inrt_row.get("ci_high"))

        print(f"\n  {g_bin}:")
        print(f"    base: D={base_d}, CI={base_ci}, status={base_row['status']}")
        print(f"    inrt: D={inrt_d}, CI={inrt_ci}, status={inrt_row['status']}")

        if base_d is not None and inrt_d is not None:
            base_sign = "positive" if base_d > 0 else "negative" if base_d < 0 else "zero"
            inrt_sign = "positive" if inrt_d > 0 else "negative" if inrt_d < 0 else "zero"
            print(f"    Знак D: base={base_sign}, inrt={inrt_sign}")
            base_zero = base_ci[0] <= 0 <= base_ci[1]
            inrt_zero = inrt_ci[0] <= 0 <= inrt_ci[1]
            print(f"    CI включает ноль: base={base_zero}, inrt={inrt_zero}")

    print("\nГотово.")
    return 0


if __name__ == "__main__":
    sys.exit(main())