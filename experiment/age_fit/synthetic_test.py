#!/usr/bin/env python3
"""
synthetic_test.py — синтетический тест по §4 протокола AGE_FIT_PROTOCOL.md.

Исправленная версия:
  - Kroupa IMF: inverse-CDF (проверен)
  - label=0,1,2 (PMS, MS, SGB) — область поворота
  - Двойные: при выходе массы за диапазон обе массы отбрасываются
  - Расстояние до отрезка (не вершины)
  - Фиксированные seed для воспроизводимости
"""

import csv, math, random, sys
from pathlib import Path

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

SIGMA_G, SIGMA_BP, SIGMA_RP = 0.005, 0.008, 0.008
BINARY_FRAC = 0.15
G_MIN, G_MAX = 11.5, 13.5
LO_MASS, MID_MASS, HI_MASS = 0.08, 0.5, 3.0
ALPHA_LO, ALPHA_HI = 1.3, 2.3

# ── Kroupa IMF (inverse CDF) ───────────────────────────────────
def pow_int(a, b, alpha):
    if alpha == 1.0:
        return math.log(b / a)
    return (b ** (1 - alpha) - a ** (1 - alpha)) / (1 - alpha)

_norm_lo = pow_int(LO_MASS, MID_MASS, ALPHA_LO)
_norm_hi = pow_int(MID_MASS, HI_MASS, ALPHA_HI) * MID_MASS
_total_norm = _norm_lo + _norm_hi
_P_SEG_LO = _norm_lo / _total_norm

def kroupa_mass(rng):
    if rng.random() < _P_SEG_LO:
        alpha = ALPHA_LO
        a, b = LO_MASS, MID_MASS
    else:
        alpha = ALPHA_HI
        a, b = MID_MASS, HI_MASS
    u = rng.random()
    pow1 = 1 - alpha
    return (u * (b ** pow1 - a ** pow1) + a ** pow1) ** (1.0 / pow1)

# ── Изохроны ───────────────────────────────────────────────────
def load_isochrone(csv_path, labels={0, 1, 2}):
    rows = []
    with open(csv_path, encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            if int(rec["label"]) not in labels:
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
    masses = [r[0] for r in iso]
    if not (masses[0] <= mass <= masses[-1]):
        raise ValueError("mass out of range")
    lo, hi = 0, len(masses) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if mass >= masses[mid] else (lo, mid)
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

def load_grid(csv_path, labels={0, 1, 2}):
    ages = {}
    with open(csv_path, encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            if int(rec["label"]) not in labels:
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

# ── Расстояние до отрезка изохроны ────────────────────────────
def distance_to_isochrone(g_star, c_star, iso):
    """Минимальное евклидово расстояние от звезды до ломаной изохроны."""
    best = float("inf")
    for i in range(len(iso) - 1):
        _, g1, c1 = iso[i]
        _, g2, c2 = iso[i + 1]
        dx, dy = g2 - g1, c2 - c1
        if dx == 0 and dy == 0:
            d = math.hypot(g_star - g1, c_star - c1)
        else:
            t = ((g_star - g1) * dx + (c_star - c1) * dy) / (dx*dx + dy*dy)
            if t < 0:
                d = math.hypot(g_star - g1, c_star - c1)
            elif t > 1:
                d = math.hypot(g_star - g2, c_star - c2)
            else:
                ix, iy = g1 + t * dx, c1 + t * dy
                d = math.hypot(g_star - ix, c_star - iy)
        if d < best:
            best = d
    return best

# ── Генерация одной реализации ─────────────────────────────────
def generate_realization(seed, synth_iso):
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
                # §4.2.3: масса вне диапазона — отбрасываем обе и выбираем заново
                continue
            # Сложение потоков
            f_g = 10 ** (-0.4 * g1) + 10 ** (-0.4 * g2)
            f_bp = 10 ** (-0.4 * bp1) + 10 ** (-0.4 * bp2)
            f_rp = 10 ** (-0.4 * rp1) + 10 ** (-0.4 * rp2)
            g_obs_comb = -2.5 * math.log10(f_g)
            bp_comb = -2.5 * math.log10(f_bp)
            rp_comb = -2.5 * math.log10(f_rp)
        else:
            g_obs_comb, bp_comb, rp_comb = g1, bp1, rp1
        # Шум Gaia
        g_obs = g_obs_comb + rng.gauss(0, SIGMA_G)
        bp_obs = bp_comb + rng.gauss(0, SIGMA_BP)
        rp_obs = rp_comb + rng.gauss(0, SIGMA_RP)
        if G_MIN <= g_obs < G_MAX:
            systems.append((g_obs, bp_obs - rp_obs))
    return systems[:N_SYN]

# ── Подбор возраста ────────────────────────────────────────────
def fit_age(stars, grid_ages):
    best_la = None
    best_m = float("inf")
    for la, iso in grid_ages.items():
        dists = [distance_to_isochrone(g, c, iso) for g, c in stars]
        dists.sort()
        m = dists[len(dists) // 2]
        if m < best_m:
            best_m = m
            best_la = la
    return best_la

# ── Main ───────────────────────────────────────────────────────
def main():
    print("=" * 72)
    print("Синтетический тест (исправленный) — AGE_FIT_PROTOCOL.md §4")
    print("=" * 72)

    # 1. Загрузка
    print("\n1. Загрузка изохрон (label=0,1,2)...")
    grid_isos, synth_isos = {}, {}
    for gname in ["Cavallo", "Bossini"]:
        cfg = GRIDS[gname]
        grid_isos[gname] = load_grid(cfg["grid_csv"])
        ages = sorted(grid_isos[gname].keys())
        print(f"  {gname}: сетка {len(ages)} возрастов: {ages[0]:.2f}–{ages[-1]:.2f}")
        synth_isos[gname] = {}
        for true_la, path in cfg["synth"].items():
            iso = load_isochrone(path)
            synth_isos[gname][true_la] = iso
            n_ms = sum(1 for r in iso if 0 <= r[3] <= 40)
            print(f"    synth {true_la:.3f}: {len(iso)} точек (вкл. subgiants)")

    # 2. Прогон
    print(f"\n2. Прогон {4 * K} реализаций (seed base=20261122)...")
    all_rows, summary, prev_fails = [], [], {}
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
                    "grid": gname, "true_logAge": f"{true_la:.3f}",
                    "realization": k, "seed": seed,
                    "estimated_logAge": f"{best_la:.2f}",
                })
            within = sum(1 for la in estimates if abs(la - true_la) <= TOLERANCE)
            passed = within >= PASS_THRESHOLD
            status = "PASS" if passed else "FAIL"
            summary.append({"grid": gname, "true_logAge": f"{true_la:.3f}",
                            "within": within, "threshold": PASS_THRESHOLD, "status": status})
            # Запоминаем предыдущие FAIL для Bossini 9.556 (первые 100)
            if gname == "Bossini" and abs(true_la - 9.556) < 0.001:
                prev_fails = {k for k, la in enumerate(estimates) if abs(la - true_la) > TOLERANCE}
                prev_ids = sorted(prev_fails)
                print(f"  {gname} true={true_la:.3f}: {within}/{K} (≥{PASS_THRESHOLD}) → {status}")
                if prev_ids:
                    print(f"    Номера реализаций вне допуска: {prev_ids}")
            else:
                print(f"  {gname} true={true_la:.3f}: {within}/{K} (≥{PASS_THRESHOLD}) → {status}")
            pair_idx += 1

    # 3. Таблица
    print("\n" + "=" * 72)
    print("ТАБЛИЦА РЕЗУЛЬТАТОВ")
    print("=" * 72)
    print(f"{'Сетка':<12} {'true logAge':>12} {'в допуске':>10} {'порог':>8} {'статус':>8}")
    print("-" * 56)
    for s in summary:
        print(f"{s['grid']:<12} {s['true_logAge']:>12} {s['within']:>4}/100    "
              f"{PASS_THRESHOLD:>3}/100   {s['status']:>8}")
    print("-" * 56)
    overall = all(s["status"] == "PASS" for s in summary)
    print(f"\nИТОГ: {'PASS' if overall else 'FAIL'}")

    # 4. Анализ предыдущих 27 неудач Bossini 9.556
    if prev_fails:
        print(f"\n3. Анализ {len(prev_fails)} неудачных реализаций Bossini (предыдущий запуск):")
        fail_rows = [r for r in all_rows if r["grid"] == "Bossini"
                     and abs(float(r["true_logAge"]) - 9.556) < 0.001
                     and int(r["realization"]) in prev_fails]
        # Группировка по оценённому возрасту
        from collections import Counter
        counter = Counter(r["estimated_logAge"] for r in fail_rows)
        print(f"  Распределение оценённых logAge:")
        for la, cnt in sorted(counter.items()):
            bar = "#" * cnt
            print(f"    logAge={la}: {cnt:>3} реализаций {bar}")
        print(f"  Ранее было 27 вне допуска (9.55, 9.60, 9.65); теперь — {len(prev_fails)}")

    # 5. Сохранение
    out_csv = AGE_DIR / "synthetic_test_results_v2.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "grid", "true_logAge", "realization", "seed", "estimated_logAge"])
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nДетальные результаты: {out_csv} ({len(all_rows)} строк)")

    return 0 if overall else 1

if __name__ == "__main__":
    sys.exit(main())