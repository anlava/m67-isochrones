#!/usr/bin/env python3
"""
synthetic_test.py — синтетический тест по §4 протокола AGE_FIT_PROTOCOL.md.

Генерирует 257 систем из готовых изохрон (synth_*), подбирает возраст
по 9 узлам сетки (grid_*). Выводит таблицу PASS/FAIL для 4 пар
(сетка A/B × истинный возраст 2.5/3.6 Грд).

Использование:
    python3 experiment/age_fit/synthetic_test.py
"""

import csv
import math
import random
import sys
from pathlib import Path

# ── конфигурация ────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]
AGE_DIR = ROOT / "experiment" / "age_fit"

GRIDS = {
    "Cavallo": {
        "grid_csv": AGE_DIR / "grid_cavallo2024_av018.csv",
        "synth": {
            9.398: AGE_DIR / "synth_cavallo2024_av018_logage9.398.csv",
            9.556: AGE_DIR / "synth_cavallo2024_av018_logage9.556.csv",
        },
    },
    "Bossini": {
        "grid_csv": AGE_DIR / "grid_bossini2019_av0115.csv",
        "synth": {
            9.398: AGE_DIR / "synth_bossini2019_av0115_logage9.398.csv",
            9.556: AGE_DIR / "synth_bossini2019_av0115_logage9.556.csv",
        },
    },
}

K = 100
N_SYN = 257
TOLERANCE = 0.05
PASS_THRESHOLD = 90

SIGMA_G = 0.005
SIGMA_BP = 0.008
SIGMA_RP = 0.008
BINARY_FRAC = 0.15
G_MIN, G_MAX = 11.5, 13.5
M_MIN, M_MAX = 0.08, 3.0


# ── изохроны ────────────────────────────────────────────────────
def load_isochrone(csv_path, label_set={0, 1}):
    """Читает CSV, возвращает список (Mini, G_app, BP_RP, G_BPmag, G_RPmag)."""
    rows = []
    with open(csv_path, encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            label = int(rec["label"])
            if label not in label_set:
                continue
            rows.append((
                float(rec["Mini"]),
                float(rec["G_app"]),
                float(rec["BP_RP"]),
                float(rec["G_BPmag"]),
                float(rec["G_RPmag"]),
            ))
    rows.sort(key=lambda x: x[0])
    return rows


def interpolate(mass, iso):
    """Линейная интерполяция по массе. Возвращает (G_app, BP_RP, G_BP, G_RP)."""
    masses = [r[0] for r in iso]
    if mass < masses[0] or mass > masses[-1]:
        raise ValueError("mass out of range")
    lo, hi = 0, len(masses) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mass < masses[mid]:
            hi = mid
        else:
            lo = mid
    m_lo, g_lo, c_lo, bp_lo, rp_lo = iso[lo]
    m_hi, g_hi, c_hi, bp_hi, rp_hi = iso[hi]
    if m_hi == m_lo:
        return g_lo, c_lo, bp_lo, rp_lo
    t = (mass - m_lo) / (m_hi - m_lo)
    return (
        g_lo + t * (g_hi - g_lo),
        c_lo + t * (c_hi - c_lo),
        bp_lo + t * (bp_hi - bp_lo),
        rp_lo + t * (rp_hi - rp_lo),
    )


def load_grid(csv_path):
    """Читает сетку, возвращает {logAge: [(Mini, G_app, BP_RP), ...]}."""
    ages = {}
    with open(csv_path, encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            label = int(rec["label"])
            if label not in (0, 1):
                continue
            la = round(float(rec["logAge"]), 2)
            if la not in ages:
                ages[la] = []
            ages[la].append((
                float(rec["Mini"]),
                float(rec["G_app"]),
                float(rec["BP_RP"]),
            ))
    for la in ages:
        ages[la].sort(key=lambda x: x[0])
    return ages


# ── генерация ───────────────────────────────────────────────────
def kroupa_mass(rng):
    """Kroupa (2001) IMF, rejection sampling на [M_MIN, M_MAX]."""
    log_lo, log_hi = math.log10(M_MIN), math.log10(M_MAX)
    while True:
        m = 10 ** (log_lo + rng.random() * (log_hi - log_lo))
        w_target = m ** (-1.3) if m < 0.5 else m ** (-2.3)
        if rng.random() < w_target / (m ** (-1.3)):
            return m


def generate_realization(seed, synth_iso):
    """
    Одна реализация: ровно N_SYN систем с (G, BP-RP).
    """
    rng = random.Random(seed)
    systems = []
    while len(systems) < N_SYN:
        m1 = kroupa_mass(rng)
        try:
            g1, c1, bp1, rp1 = interpolate(m1, synth_iso)
        except ValueError:
            continue
        is_binary = rng.random() < BINARY_FRAC
        if is_binary:
            q = 0.5 + rng.random() * 0.5
            m2 = q * m1
            try:
                g2, c2, bp2, rp2 = interpolate(m2, synth_iso)
            except ValueError:
                is_binary = False
        if is_binary:
            # Сложение потоков в G, BP, RP независимо
            f_g = 10 ** (-0.4 * g1) + 10 ** (-0.4 * g2)
            f_bp = 10 ** (-0.4 * bp1) + 10 ** (-0.4 * bp2)
            f_rp = 10 ** (-0.4 * rp1) + 10 ** (-0.4 * rp2)
            g_comb = -2.5 * math.log10(f_g)
            bp_comb = -2.5 * math.log10(f_bp)
            rp_comb = -2.5 * math.log10(f_rp)
        else:
            g_comb, bp_comb, rp_comb = g1, bp1, rp1
        # Шум
        g_obs = g_comb + rng.gauss(0, SIGMA_G)
        bp_obs = bp_comb + rng.gauss(0, SIGMA_BP)
        rp_obs = rp_comb + rng.gauss(0, SIGMA_RP)
        # Фильтр G
        if G_MIN <= g_obs < G_MAX:
            systems.append((g_obs, bp_obs - rp_obs))
    return systems[:N_SYN]


# ── подбор возраста ─────────────────────────────────────────────
def fit_age(stars, grid_ages):
    """
    По набору звёзд stars = [(G, BP-RP)] и словарю grid_ages
    {logAge: [(Mini, G_app, BP_RP), ...]} возвращает logAge с
    минимальной медианой расстояний.
    """
    best_la = None
    best_m = float("inf")
    for la, iso in grid_ages.items():
        dists = []
        for g_star, c_star in stars:
            # Ищем ближайшую точку на изохроне
            min_d = float("inf")
            for _, g_mod, c_mod in iso:
                d = math.sqrt((g_star - g_mod) ** 2 + (c_star - c_mod) ** 2)
                if d < min_d:
                    min_d = d
            dists.append(min_d)
        dists.sort()
        m = dists[len(dists) // 2]  # медиана
        if m < best_m:
            best_m = m
            best_la = la
    return best_la


# ── main ────────────────────────────────────────────────────────
def main():
    print("=" * 72)
    print("Синтетический тест — AGE_FIT_PROTOCOL.md §4")
    print("=" * 72)

    # 1. Загрузка данных
    print("\n1. Загрузка изохрон...")
    grid_isos = {}
    synth_isos = {}
    for gname, cfg in GRIDS.items():
        grid_isos[gname] = load_grid(cfg["grid_csv"])
        ages = sorted(grid_isos[gname].keys())
        print(f"  {gname}: сетка {len(ages)} возрастов: "
              f"{ages[0]:.2f}–{ages[-1]:.2f}")
        synth_isos[gname] = {}
        for true_la, path in cfg["synth"].items():
            iso = load_isochrone(path)
            synth_isos[gname][true_la] = iso
            print(f"    synth {true_la:.3f}: {len(iso)} точек")

    # 2. Прогон
    print(f"\n2. Прогон {4 * K} реализаций...")
    all_rows = []
    summary = []

    # Для воспроизводимости: детерминированная последовательность seed
    base_seed = 20261122

    pair_idx = 0
    for gname in ["Cavallo", "Bossini"]:
        for true_la in sorted(GRIDS[gname]["synth"].keys()):
            synth_iso = synth_isos[gname][true_la]
            grid = grid_isos[gname]

            estimates = []
            for k in range(K):
                seed = base_seed + pair_idx * 1000 + k
                stars = generate_realization(seed, synth_iso)
                best_la = fit_age(stars, grid)
                estimates.append(best_la)
                all_rows.append({
                    "grid": gname,
                    "true_logAge": f"{true_la:.3f}",
                    "realization": k,
                    "seed": seed,
                    "estimated_logAge": f"{best_la:.2f}",
                })

            within = sum(1 for la in estimates if abs(la - true_la) <= TOLERANCE)
            passed = within >= PASS_THRESHOLD
            status = "PASS" if passed else "FAIL"
            summary.append({
                "grid": gname,
                "true_logAge": f"{true_la:.3f}",
                "within": within,
                "threshold": PASS_THRESHOLD,
                "status": status,
            })
            print(f"  {gname} true={true_la:.3f}: {within}/{K} "
                  f"(≥{PASS_THRESHOLD}) → {status}")
            pair_idx += 1

    # 3. Таблица результатов
    print("\n" + "=" * 72)
    print("ТАБЛИЦА РЕЗУЛЬТАТОВ")
    print("=" * 72)
    print(f"{'Сетка':<12} {'true logAge':>12} {'в допуске':>10} "
          f"{'порог':>8} {'статус':>8}")
    print("-" * 56)
    for s in summary:
        print(f"{s['grid']:<12} {s['true_logAge']:>12} "
              f"{s['within']:>4}/100    {PASS_THRESHOLD:>3}/100   "
              f"{s['status']:>8}")
    print("-" * 56)
    overall = all(s["status"] == "PASS" for s in summary)
    print(f"\nИТОГ: {'PASS' if overall else 'FAIL'}")
    print()

    # 4. Сохранение CSV
    out_csv = AGE_DIR / "synthetic_test_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "grid", "true_logAge", "realization", "seed", "estimated_logAge"])
        w.writeheader()
        w.writerows(all_rows)
    print(f"Детальные результаты: {out_csv} ({len(all_rows)} строк)")

    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())