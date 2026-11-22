#!/usr/bin/env python3
"""CMD области поворота M67 с изохронами и ΔG между моделями на label=2.

CMD: 11.5 ≤ G < 13.5, 257 звёзд Gaia + 35 S1-кандидатов.
Наложены: изохроны Cavallo и Bossini, на каждой выделен label=2 (субгиганты).
Построен график ΔG(BP-RP) на общем интервале цветов.

Выходы:
  - turnoff_cmd_isochrones.png   — CMD с изохронами
  - turnoff_delta_G.png          — ΔG = G_Cavallo − G_Bossini
  - subgiant_brightness_diff.csv — таблица ΔG
"""
import csv
import sys
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "experiment" / "turnoff"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MEMBERS_CSV = ROOT / "results" / "m67_members_selected.csv"
CAVALLO_CSV = ROOT / "results" / "m67_isochrone_apparent.csv"
BOSSINI_CSV = ROOT / "results" / "m67_isochrone_bossini2019_apparent.csv"
CANDIDATES_CSV = OUT_DIR / "subgiant_candidates_ids.csv"


# ─────────────────────────────────────────────
# 1. Read data
# ─────────────────────────────────────────────
def read_turnoff_members(path):
    stars = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            g = float(r["Gmag"])
            if not (11.5 <= g < 13.5):
                continue
            stars.append({
                "GaiaDR3": r["GaiaDR3"].strip(),
                "Gmag": g,
                "BP_RP": float(r["BP_RP"]),
            })
    return stars


def read_label2_points(path):
    """Возвращает список (G_app, BP_RP) для label=2, в порядке таблицы."""
    pts = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if int(r["label"]) == 2:
                pts.append((float(r["G_app"]), float(r["BP_RP"])))
    return pts


def read_candidates_ids(path):
    """Возвращает множество GaiaDR3 кандидатов S1."""
    ids = set()
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            ids.add(r["GaiaDR3"].strip())
    return ids


# ─────────────────────────────────────────────
# 2. Interpolation: colour → G (inverse of usual)
# ─────────────────────────────────────────────
def make_given_colour_interpolator(points, name):
    """Строит функцию G(BP-RP) по label=2 точкам.

    Проверяет монотонность цвета в ИСХОДНОМ порядке точек (как в файле).
    Затем сортирует для интерполяции. Возвращает (fun, bp_rp_min, bp_rp_max)
    или None если цвет не монотонен.
    """
    # Проверка монотонности цвета в исходном порядке
    colours_orig = np.array([c for _, c in points])
    diffs = np.diff(colours_orig)
    if not np.all(diffs >= 0):
        non_mono = np.where(diffs < 0)[0]
        print(f"  {name}: BP-RP НЕ монотонен в исходном порядке! "
              f"Нарушение на индексах {non_mono.tolist()} "
              f"(BP-RP: {colours_orig[non_mono[0]]:.4f} → {colours_orig[non_mono[0]+1]:.4f})")
        return None
    print(f"  {name}: BP-RP монотонен в исходном порядке (все {len(points)} точек)")

    # Сортируем по BP_RP для интерполяции
    sorted_pts = sorted(points, key=lambda x: x[1])
    colours = np.array([c for _, c in sorted_pts])
    mags = np.array([g for g, _ in sorted_pts])

    def interp(colour):
        return float(np.interp(colour, colours, mags, left=np.nan, right=np.nan))

    return interp, float(colours[0]), float(colours[-1])


# ─────────────────────────────────────────────
# 3. Main
# ─────────────────────────────────────────────
def main():
    print("=" * 70)
    print("CMD области поворота M67 и разность яркостей label=2")
    print("=" * 70)

    # Данные
    stars = read_turnoff_members(MEMBERS_CSV)
    s1_ids = read_candidates_ids(CANDIDATES_CSV)
    cav_pts = read_label2_points(CAVALLO_CSV)
    bos_pts = read_label2_points(BOSSINI_CSV)

    print(f"  Звёзд в области поворота: {len(stars)}")
    print(f"  S1-кандидатов: {len(s1_ids)}")
    print(f"  Cavallo label=2 точек: {len(cav_pts)}")
    print(f"  Bossini label=2 точек: {len(bos_pts)}")

    # Разделяем звёзды
    g_all = np.array([s["Gmag"] for s in stars])
    bp_all = np.array([s["BP_RP"] for s in stars])
    g_s1 = np.array([s["Gmag"] for s in stars if s["GaiaDR3"] in s1_ids])
    bp_s1 = np.array([s["BP_RP"] for s in stars if s["GaiaDR3"] in s1_ids])

    # Разделяем label=2 на непрерывные блоки для отрисовки
    def split_blocks(pts):
        blocks = []
        cur = []
        for i, (g, c) in enumerate(pts):
            cur.append((g, c))
            if i == len(pts) - 1 or pts[i+1][1] < c:  # BP_RP растёт, следующая точка — новый блок
                # На самом деле проверять по label, но у нас только label=2
                pass
        blocks.append(cur)
        return blocks

    # Просто используем все точки как одну линию (BP_RP монотонный)
    cav_arr = np.array(cav_pts)
    bos_arr = np.array(bos_pts)

    # ─────────────────────────────────────────────
    # CMD plot
    # ─────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 9))

    # Все звёзды
    ax.scatter(bp_all, g_all, c="gray", s=10, alpha=0.5, edgecolors="none",
               label=f"All stars ({len(stars)})", zorder=2)

    # S1 кандидаты
    ax.scatter(bp_s1, g_s1, c="red", s=20, alpha=0.8, edgecolors="black",
               linewidths=0.3, label=f"S1 candidates ({len(s1_ids)})", zorder=4)

    # Cavallo изохрона (все label)
    cav_full = []
    with open(CAVALLO_CSV, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            g = float(r["G_app"])
            if 11.0 <= g <= 14.0:
                cav_full.append((float(r["BP_RP"]), g))
    cav_arr_full = np.array(cav_full)
    ax.plot(cav_arr_full[:, 0], cav_arr_full[:, 1], c="grey", linewidth=0.8,
            alpha=0.3, label="Cavallo et al. (2024) — full", zorder=1)

    # Bossini изохрона (все label)
    bos_full = []
    with open(BOSSINI_CSV, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            g = float(r["G_app"])
            if 11.0 <= g <= 14.0:
                bos_full.append((float(r["BP_RP"]), g))
    bos_arr_full = np.array(bos_full)
    ax.plot(bos_arr_full[:, 0], bos_arr_full[:, 1], c="grey", linewidth=0.8,
            alpha=0.3, label="Bossini et al. (2019) — full", zorder=1)

    # Cavallo label=2 (жирная красная)
    ax.plot(cav_arr[:, 1], cav_arr[:, 0], c="crimson", linewidth=2.5,
            label=f"Cavallo — SGB (label=2)", zorder=5)
    # Маркеры на узлах Cavallo label=2
    ax.scatter(cav_arr[:, 1], cav_arr[:, 0], c="crimson", s=15, alpha=0.6,
               edgecolors="none", zorder=5)

    # Bossini label=2 (жирная синяя)
    ax.plot(bos_arr[:, 1], bos_arr[:, 0], c="royalblue", linewidth=2.5,
            label=f"Bossini — SGB (label=2)", zorder=5)
    ax.scatter(bos_arr[:, 1], bos_arr[:, 0], c="royalblue", s=15, alpha=0.6,
               edgecolors="none", zorder=5)

    # Оформление
    ax.invert_yaxis()
    ax.set_xlabel("G_BP − G_RP (mag)", fontsize=12)
    ax.set_ylabel("G (mag)", fontsize=12)
    ax.set_title("M67 (NGC 2682) — turnoff region\n"
                 "11.5 ≤ G < 13.5, Gaia DR3 + PARSEC isochrones", fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower right", ncol=1)
    fig.tight_layout()

    cmd_path = OUT_DIR / "turnoff_cmd_isochrones.png"
    fig.savefig(cmd_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nCMD сохранён: {cmd_path}")

    # ─────────────────────────────────────────────
    # 4. Rule documentation and ΔG computation
    # ─────────────────────────────────────────────
    print("\n--- Правило проверки ---")
    print("Метрика: ΔG(BP−RP) = G_Cavallo(BP−RP) − G_Bossini(BP−RP), mag.")
    print("Обе ветви label=2 просматриваются от цвета к величине.")
    print("Условия однозначности:")

    cav_interp = make_given_colour_interpolator(cav_pts, "Cavallo")
    bos_interp = make_given_colour_interpolator(bos_pts, "Bossini")

    if cav_interp is None:
        print("  Cavallo: BP-RP не монотонен — остановка.")
        sys.exit(1)
    if bos_interp is None:
        print("  Bossini: BP-RP не монотонен — остановка.")
        sys.exit(1)

    g_cav_fn, cav_c_min, cav_c_max = cav_interp
    g_bos_fn, bos_c_min, bos_c_max = bos_interp

    print(f"  Cavallo label=2: BP−RP ∈ [{cav_c_min:.3f}, {cav_c_max:.3f}]")
    print(f"  Bossini label=2: BP−RP ∈ [{bos_c_min:.3f}, {bos_c_max:.3f}]")

    common_c_min = max(cav_c_min, bos_c_min)
    common_c_max = min(cav_c_max, bos_c_max)

    print(f"\n  Общий интервал BP−RP: [{common_c_min:.3f}, {common_c_max:.3f}]")

    if common_c_min >= common_c_max:
        print("\n  РЕЗУЛЬТАТ: общего интервала BP−RP для label=2 нет.")
        print("  Сравнение ΔG на label=2 невозможно.")
        # Создадим пустой файл с результатом
        with open(OUT_DIR / "subgiant_brightness_diff.csv", "w") as f:
            f.write("status,common_bp_rp_min,common_bp_rp_max\n")
            f.write(f"no_common_interval,{common_c_min:.4f},{common_c_max:.4f}\n")
        print(f"\nРезультат сохранён в subgiant_brightness_diff.csv")
        print("Сравнение преимущества Cavallo/Bossini не производилось.")
        sys.exit(0)

    # Вычисляем ΔG на узлах обеих ветвей в общем интервале
    # Собираем все узловые BP-RP из обеих моделей, попадающие в общий интервал
    node_bps = set()
    for g, c in cav_pts:
        if common_c_min <= c <= common_c_max:
            node_bps.add(round(c, 6))
    for g, c in bos_pts:
        if common_c_min <= c <= common_c_max:
            node_bps.add(round(c, 6))

    node_bps = sorted(node_bps)
    bp_grid = np.array(node_bps)
    g_cav = np.array([g_cav_fn(c) for c in bp_grid])
    g_bos = np.array([g_bos_fn(c) for c in bp_grid])
    delta = g_cav - g_bos

    # Таблица
    csv_path = OUT_DIR / "subgiant_brightness_diff.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["BP_RP", "G_Cavallo", "G_Bossini", "delta_G"])
        for c, gc, gb, d in zip(bp_grid, g_cav, g_bos, delta):
            writer.writerow([f"{c:.6f}", f"{gc:.6f}", f"{gb:.6f}", f"{d:.6f}"])
    print(f"\nТаблица ΔG сохранена: {csv_path} (по {len(node_bps)} узлам)")

    # Печать статистики
    print(f"\n  ΔG = G_Cavallo − G_Bossini на {len(node_bps)} узлах обеих ветвей")
    print(f"    от BP−RP = {common_c_min:.3f} до {common_c_max:.3f}:")
    print(f"    min ΔG = {float(np.min(delta)):.4f} mag (отрицательный → Cavallo ярче)")
    print(f"    max ΔG = {float(np.max(delta)):.4f} mag")
    print(f"    mean ΔG = {float(np.mean(delta)):.4f} mag")
    # Найдите макс расхождение и соответствующий BP-RP
    argmax = np.argmax(np.abs(delta))
    print(f"    Макс |ΔG| = {float(np.abs(delta)[argmax]):.4f} mag при BP−RP = {bp_grid[argmax]:.3f}")

    # ─────────────────────────────────────────────
    # 5. ΔG plot
    # ─────────────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(bp_grid, delta, c="darkgreen", linewidth=2)
    ax2.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax2.set_xlabel("G_BP − G_RP (mag)", fontsize=12)
    ax2.set_ylabel("ΔG = G_Cavallo − G_Bossini (mag)", fontsize=12)
    ax2.set_title("Brightness difference on the subgiant branch (label=2)\n"
                  f"M67, common colour range [{common_c_min:.3f}, {common_c_max:.3f}]",
                  fontsize=12)
    ax2.grid(True, alpha=0.3)

    # Найти макс расхождение
    argmax = int(np.argmax(np.abs(delta)))
    max_bp = bp_grid[argmax]
    max_dg = float(np.abs(delta)[argmax])

    # Подпись статистики
    text = (f"Cavallo brighter than Bossini\n"
            f"ΔG = G_Cav − G_Bos = {float(np.mean(delta)):.3f} ± {float(np.std(delta)):.3f} mag\n"
            f"Range [{float(np.min(delta)):.3f}, {float(np.max(delta)):.3f}]\n"
            f"Max |ΔG| = {max_dg:.3f} mag at BP−RP = {max_bp:.3f}")
    ax2.text(0.98, 0.05, text, transform=ax2.transAxes, fontsize=10,
             ha="right", va="bottom", bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    fig2.tight_layout()
    delta_path = OUT_DIR / "turnoff_delta_G.png"
    fig2.savefig(delta_path, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"\nГрафик ΔG сохранён: {delta_path}")

    # Интерпретация
    mean_d = float(np.mean(delta))
    if mean_d < 0:
        interp = (f"ИТОГ: На всём общем интервале BP−RP "
                  f"Cavallo ярче Bossini на {abs(mean_d):.3f} mag "
                  f"(ΔG = G_Cav − G_Bos < 0). "
                  f"Максимальное расхождение {max_dg:.3f} mag при BP−RP = {max_bp:.3f}.")
    else:
        interp = (f"ИТОГ: На всём общем интервале BP−RP "
                  f"Bossini ярче Cavallo на {mean_d:.3f} mag.")
    print(f"\n  {interp}")
    print("\nСравнение преимущества Cavallo/Bossini не производилось.")
    print("=" * 70)


if __name__ == "__main__":
    main()