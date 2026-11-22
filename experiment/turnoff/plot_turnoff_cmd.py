#!/usr/bin/env python3
"""CMD области поворота M67 (11.5 ≤ G < 13.5) и отбор кандидатов в субгиганты.

Строит диаграмму цвет–величина только по Gaia DR3, без моделей.
Применяет правило S1 (скользящая медиана + 0.10 mag) и S2 (явное сечение).

Выходы:
  - turnoff_cmd_gaia_only.png  — CMD
  - subgiant_candidates_ids.csv — кандидаты по S1
"""
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MEMBERS_CSV = ROOT / "results" / "m67_members_selected.csv"
OUT_DIR = ROOT / "experiment" / "turnoff"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PNG_PATH = OUT_DIR / "turnoff_cmd_gaia_only.png"
CSV_PATH = OUT_DIR / "subgiant_candidates_ids.csv"


def read_members_csv(path):
    """Читает results/m67_members_selected.csv."""
    stars = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            g = float(row["Gmag"])
            if not (11.5 <= g < 13.5):
                continue
            stars.append({
                "GaiaDR3": row["GaiaDR3"].strip(),
                "Gmag": g,
                "BP_RP": float(row["BP_RP"]),
            })
    return stars


def rolling_median_with_threshold(stars, window=0.5, threshold=0.10):
    """Правило S1: скользящая медиана + порог.

    Для каждой звезды: если BP_RP > медиана(BP_RP) в окне [G-window/2, G+window/2]
    + threshold, то кандидат.
    """
    sorted_stars = sorted(stars, key=lambda s: s["Gmag"])
    g_vals = np.array([s["Gmag"] for s in sorted_stars])
    bp_rp_vals = np.array([s["BP_RP"] for s in sorted_stars])

    selected = []
    for i, s in enumerate(sorted_stars):
        g = s["Gmag"]
        lo, hi = g - window / 2, g + window / 2
        mask = (g_vals >= lo) & (g_vals <= hi)
        if not np.any(mask):
            continue
        median = float(np.median(bp_rp_vals[mask]))
        if s["BP_RP"] > median + threshold:
            selected.append(s)
    return selected


def explicit_cut_rule(stars):
    """Правило S2: явное сечение по BP_RP.

    BP_RP > 1.00 для G < 12.5
    BP_RP > 0.90 для 12.5 ≤ G < 13.5
    """
    selected = []
    for s in stars:
        if s["Gmag"] < 12.5 and s["BP_RP"] > 1.00:
            selected.append(s)
        elif 12.5 <= s["Gmag"] < 13.5 and s["BP_RP"] > 0.90:
            selected.append(s)
    return selected


def plot_cmd(stars, s1_selected, s2_selected, out_path):
    """Строит CMD области поворота с выделением субгигантов."""
    fig, ax = plt.subplots(figsize=(8, 9))

    # Все звёзды в области (серые)
    g_all = np.array([s["Gmag"] for s in stars])
    bp_rp_all = np.array([s["BP_RP"] for s in stars])
    ax.scatter(bp_rp_all, g_all, c="gray", s=10, alpha=0.5,
               edgecolors="none", label=f"Все звёзды ({len(stars)})", zorder=2)

    # S2 (явное сечение) — светло-оранжевые
    if s2_selected:
        g_s2 = np.array([s["Gmag"] for s in s2_selected])
        bp_rp_s2 = np.array([s["BP_RP"] for s in s2_selected])
        ax.scatter(bp_rp_s2, g_s2, c="orange", s=15, alpha=0.6,
                   edgecolors="none", label=f"Правило S2 ({len(s2_selected)})", zorder=3)

    # S1 (скользящая медиана) — красные
    if s1_selected:
        g_s1 = np.array([s["Gmag"] for s in s1_selected])
        bp_rp_s1 = np.array([s["BP_RP"] for s in s1_selected])
        ax.scatter(bp_rp_s1, g_s1, c="red", s=20, alpha=0.8,
                   edgecolors="black", linewidths=0.3,
                   label=f"Субгиганты S1 ({len(s1_selected)})", zorder=4)

    # Оформление
    ax.invert_yaxis()
    ax.set_xlabel("G_BP − G_RP (mag)", fontsize=12)
    ax.set_ylabel("G (mag)", fontsize=12)
    ax.set_title("M67 (NGC 2682) — область поворота\n11.5 ≤ G < 13.5, только Gaia DR3", fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  CMD сохранён: {out_path}")


def save_candidates_csv(candidates, out_path):
    """Сохраняет список кандидатов в CSV."""
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["GaiaDR3", "Gmag", "BP_RP", "rule"])
        for s in candidates:
            writer.writerow([s["GaiaDR3"], f"{s['Gmag']:.8f}",
                             f"{s['BP_RP']:.8f}", "S1"])
    print(f"  Кандидаты сохранены: {out_path} ({len(candidates)} звёзд)")


def main():
    print("=" * 60)
    print("CMD области поворота M67 и отбор субгигантов")
    print("=" * 60)

    # Чтение
    stars = read_members_csv(MEMBERS_CSV)
    print(f"Звёзд с 11.5 ≤ G < 13.5: {len(stars)}")

    if not stars:
        print("ОШИБКА: нет звёзд в заданном диапазоне G.")
        sys.exit(1)

    # Правила отбора
    s1_selected = rolling_median_with_threshold(stars)
    print(f"Правило S1 (скользящая медиана + 0.10): {len(s1_selected)} кандидатов")

    s2_selected = explicit_cut_rule(stars)
    print(f"Правило S2 (явное сечение): {len(s2_selected)} кандидатов")

    # Пересечение
    s1_ids = {s["GaiaDR3"] for s in s1_selected}
    s2_ids = {s["GaiaDR3"] for s in s2_selected}
    common = s1_ids & s2_ids
    print(f"Пересечение S1 ∩ S2: {len(common)}")
    print(f"Только S1: {len(s1_ids - s2_ids)}")
    print(f"Только S2: {len(s2_ids - s1_ids)}")

    # Проверка достаточности
    if len(s1_selected) < 10:
        print("ПРЕДУПРЕЖДЕНИЕ: кандидатов < 10 — недостаточно для сравнения моделей.")

    # Построение CMD
    print("\nПостроение CMD...")
    plot_cmd(stars, s1_selected, s2_selected, PNG_PATH)

    # Сохранение
    save_candidates_csv(s1_selected, CSV_PATH)

    print("\nГотово. Сравнение моделей Cavallo/Bossini не производилось.")
    print("=" * 60)


if __name__ == "__main__":
    main()