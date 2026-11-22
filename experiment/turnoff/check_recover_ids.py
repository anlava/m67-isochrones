#!/usr/bin/env python3
"""Проверка 1: восстановить набор GaiaDR3 субгигантов из исходного TSV.

Читает data/raw/m67_members_hr2023.tsv независимым парсером,
применяет фильтры (Name == NGC_2682, конечные G/BP/RP, 11.5 ≤ G < 13.5)
и правило S1 (скользящая медиана + 0.10 mag), затем сравнивает
полученный набор ID с experiment/turnoff/subgiant_candidates_ids.csv.

Ожидание: совпадение множеств GaiaDR3.
"""
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW_TSV = ROOT / "data" / "raw" / "m67_members_hr2023.tsv"
CANDIDATES_CSV = ROOT / "experiment" / "turnoff" / "subgiant_candidates_ids.csv"


def parse_mag(val_str, label):
    """Парсит величину; возвращает nan при пустом или некорректном значении."""
    s = val_str.strip()
    if not s:
        return math.nan
    try:
        v = float(s)
    except ValueError:
        return math.nan
    if not math.isfinite(v):
        return math.nan
    return v


def read_raw_tsv(path):
    """Независимый парсер TSV (не использует plot_m67.read_members_tsv)."""
    members = []
    with open(path, encoding="utf-8-sig") as f:
        lines = f.readlines()
    # Убираем комментарии
    data = [l for l in lines if not l.lstrip().startswith("#")]
    # Ищем заголовок
    header_idx = None
    for i, line in enumerate(data):
        if line.startswith("Name\t") and "GaiaDR3" in line:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Заголовок не найден")
    # Пропускаем строки единиц и разделителей
    start = header_idx + 1
    while start < len(data):
        stripped = data[start].strip()
        if not stripped:
            start += 1
            continue
        cols = data[start].rstrip("\r\n").split("\t")
        if len(cols) == 7 and all(c.strip() == "mag" for c in cols[4:7]):
            start += 1
            continue
        if len(cols) == 7 and all(set(c.strip()).issubset({"-"}) for c in cols):
            start += 1
            continue
        break
    for line in data[start:]:
        if not line.strip():
            continue
        cols = line.rstrip("\r\n").split("\t")
        if len(cols) != 7:
            continue
        name = cols[0].strip()
        gaia = cols[1].strip()
        gmag = parse_mag(cols[4], "Gmag")
        bpmag = parse_mag(cols[5], "BPmag")
        rpmag = parse_mag(cols[6], "RPmag")
        members.append({
            "Name": name,
            "GaiaDR3": gaia,
            "Gmag": gmag,
            "BPmag": bpmag,
            "RPmag": rpmag,
        })
    return members


def select_turnoff(members):
    """Фильтр: Name == NGC_2682, конечные величины, 11.5 ≤ G < 13.5."""
    result = []
    for m in members:
        if m["Name"] not in ("NGC_2682", "NGC2682"):
            continue
        if not (math.isfinite(m["Gmag"]) and math.isfinite(m["BPmag"])
                and math.isfinite(m["RPmag"])):
            continue
        if not (11.5 <= m["Gmag"] < 13.5):
            continue
        result.append(m)
    return result


def rule_S1(stars):
    """Правило S1: скользящая медиана + 0.10 mag.

    Окно [G-0.25, G+0.25] (обе границы включены).
    Медиана — numpy.median (среднее двух центральных при чётном N).
    """
    import numpy as np  # импорт здесь, чтобы не тянуть numpy на весь скрипт

    if not stars:
        return set()
    sorted_stars = sorted(stars, key=lambda s: s["Gmag"])
    g_vals = np.array([s["Gmag"] for s in sorted_stars])
    bp_rp_vals = np.array([s["BPmag"] - s["RPmag"] for s in sorted_stars])

    selected = set()
    window = 0.5  # mag

    for i, s in enumerate(sorted_stars):
        g = s["Gmag"]
        lo, hi = g - window / 2, g + window / 2
        mask = (g_vals >= lo) & (g_vals <= hi)
        if not np.any(mask):
            continue
        median = float(np.median(bp_rp_vals[mask]))
        bp_rp = s["BPmag"] - s["RPmag"]
        if bp_rp > median + 0.10:
            selected.add(s["GaiaDR3"])
    return selected


def read_candidates_csv(path):
    """Читает CSV кандидатов и возвращает множество GaiaDR3."""
    ids = set()
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids.add(row["GaiaDR3"].strip())
    return ids


def main():
    print("=" * 60)
    print("Проверка 1: восстановление набора ID из исходного TSV")
    print("=" * 60)

    # Читаем исходный TSV
    try:
        members = read_raw_tsv(RAW_TSV)
        print(f"Прочитано строк из TSV: {len(members)}")
    except (OSError, ValueError) as exc:
        print(f"ОШИБКА: не удалось прочитать TSV: {exc}")
        sys.exit(1)

    # Отбираем область поворота
    turnoff = select_turnoff(members)
    print(f"Звёзд с 11.5 ≤ G < 13.5: {len(turnoff)}")

    # Применяем правило S1
    recovered = rule_S1(turnoff)
    print(f"Кандидатов по правилу S1: {len(recovered)}")

    # Читаем эталонный набор
    if not CANDIDATES_CSV.exists():
        print(f"ОШИБКА: эталонный файл не найден: {CANDIDATES_CSV}")
        sys.exit(1)
    reference = read_candidates_csv(CANDIDATES_CSV)
    print(f"Эталонных ID из CSV: {len(reference)}")

    # Сравнение
    if recovered == reference:
        print("\nРЕЗУЛЬТАТ: PASS — набор ID полностью восстановлен.")
        sys.exit(0)
    else:
        only_recovered = recovered - reference
        only_ref = reference - recovered
        print(f"\nРЕЗУЛЬТАТ: FAIL — наборы не совпадают.")
        if only_recovered:
            print(f"  В восстановленном, но не в эталоне ({len(only_recovered)}):")
            for gaia in sorted(only_recovered)[:5]:
                print(f"    {gaia}")
            if len(only_recovered) > 5:
                print(f"    ... и ещё {len(only_recovered)-5}")
        if only_ref:
            print(f"  В эталоне, но не в восстановленном ({len(only_ref)}):")
            for gaia in sorted(only_ref)[:5]:
                print(f"    {gaia}")
            if len(only_ref) > 5:
                print(f"    ... и ещё {len(only_ref)-5}")
        sys.exit(1)


if __name__ == "__main__":
    main()