#!/usr/bin/env python3
"""Проверка покрытия label=2 для 35 кандидатов S1.

Для каждого кандидата (из subgiant_candidates_ids.csv) проверяет,
сколько значений цвета даёт label=2 каждой модели (Cavallo, Bossini)
при G = Gmag кандидата без экстраполяции.

Выход: таблица по модели: число случаев с 0, 1, >1 цветами.
"""
import csv
import sys
from fractions import Fraction
# Python 3.9 compatibility: pairwise added in 3.10
try:
    from itertools import pairwise
except ImportError:
    def pairwise(iterable):
        it = iter(iterable)
        try:
            prev = next(it)
        except StopIteration:
            return
        for curr in it:
            yield (prev, curr)
            prev = curr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CANDIDATES_CSV = ROOT / "experiment" / "turnoff" / "subgiant_candidates_ids.csv"
CAVALLO_CSV = ROOT / "results" / "m67_isochrone_apparent.csv"
BOSSINI_CSV = ROOT / "results" / "m67_isochrone_bossini2019_apparent.csv"


def read_candidates(path):
    """Читает список кандидатов S1, возвращает список (GaiaDR3, Gmag, BP_RP)."""
    result = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            result.append((
                row["GaiaDR3"].strip(),
                Fraction(row["Gmag"]),
                Fraction(row["BP_RP"]),
            ))
    return result


def read_label2_segments(path):
    """Читает CSV изохроны, возвращает список отрезков label=2.

    Каждый отрезок: ((G_left, BP_RP_left), (G_right, BP_RP_right)).
    Соседние точки с label=2 образуют отрезок.
    """
    segments = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        points = [(Fraction(r["G_app"]), Fraction(r["BP_RP"]), int(r["label"]))
                  for r in reader]

    # Собираем непрерывные блоки label=2
    block = []
    for g, c, lbl in points:
        if lbl == 2:
            block.append((g, c))
        else:
            if len(block) >= 2:
                for (g1, c1), (g2, c2) in pairwise(block):
                    segments.append(((g1, c1), (g2, c2)))
            block = []
    if len(block) >= 2:
        for (g1, c1), (g2, c2) in pairwise(block):
            segments.append(((g1, c1), (g2, c2)))
    return segments


def colours_at_g(g, segments):
    """Сколько разных цветов дают отрезки label=2 при данной G.

    Возвращает (число_цветов, множество_цветов).
    """
    colours = set()
    for (g1, c1), (g2, c2) in segments:
        left, right = (g1, g2) if g1 <= g2 else (g2, g1)
        # Проверяем: G внутри отрезка (включая границы)
        if not (left <= g <= right):
            continue
        # Экстраполяция запрещена
        if g < min(g1, g2) or g > max(g1, g2):
            continue
        if g1 == g2:
            colours.add(c1)
            colours.add(c2)
        else:
            # Линейная интерполяция
            c = c1 + (c2 - c1) * (g - g1) / (g2 - g1)
            colours.add(c)
    return len(colours), colours


def summarize(segments, candidates, model_name):
    """Для всех кандидатов подсчитывает случаи 0, 1, >1 цветов."""
    counts = {0: 0, 1: 0, "multi": 0}
    details = []
    for gaia, g, obs_colour in candidates:
        n, cols = colours_at_g(g, segments)
        key = n if n <= 1 else "multi"
        counts[key] += 1
        details.append((gaia, g, n, cols))
    return counts, details


def main():
    print("=" * 70)
    print("Проверка покрытия label=2 для 35 кандидатов S1")
    print("=" * 70)

    candidates = read_candidates(CANDIDATES_CSV)
    print(f"\nКандидатов S1: {len(candidates)}")
    print(f"  Gmag: [{float(min(c[1] for c in candidates)):.3f}, "
          f"{float(max(c[1] for c in candidates)):.3f}]")
    print()

    for model_name, csv_path in [("Cavallo", CAVALLO_CSV), ("Bossini", BOSSINI_CSV)]:
        segments = read_label2_segments(csv_path)
        print(f"  {model_name}: {len(segments)} отрезков label=2")

        # Покажем диапазон покрытия label=2
        if segments:
            all_g = [g for seg in segments for (g, _) in seg]
            all_c = [c for seg in segments for (_, c) in seg]
            print(f"    G_app label=2: [{float(min(all_g)):.3f}, {float(max(all_g)):.3f}]")
            print(f"    BP_RP label=2: [{float(min(all_c)):.3f}, {float(max(all_c)):.3f}]")

        counts, details = summarize(segments, candidates, model_name)
        print(f"\n  --- {model_name} ---")
        print(f"  Цветов = 0 (нет label=2 при G): {counts[0]}")
        print(f"  Цветов = 1 (однозначно):         {counts[1]}")
        print(f"  Цветов > 1 (неоднозначно):       {counts['multi']}")

        # Подробности для случаев != 1
        problem_stars = [(g, n, cols) for _, g, n, cols in details if n != 1]
        if counts[0] > 0:
            print(f"\n  Звёзды без покрытия label=2:")
            for gaia, g, n, cols in details:
                if n == 0:
                    g_float = float(g)
                    print(f"    GaiaDR3={gaia}  G={g_float:.3f}")
        if counts["multi"] > 0:
            print(f"\n  Звёзды с несколькими цветами label=2:")
            for gaia, g, n, cols in details:
                if n > 1:
                    g_float = float(g)
                    cols_str = ", ".join(f"{float(c):.4f}" for c in cols)
                    print(f"    GaiaDR3={gaia}  G={g_float:.3f}  цветов={n}: [{cols_str}]")
        print()

    print("Сравнение преимущества Cavallo/Bossini не производилось.")


if __name__ == "__main__":
    main()