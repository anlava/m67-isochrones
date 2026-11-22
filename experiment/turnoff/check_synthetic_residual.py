#!/usr/bin/env python3
"""Проверка 2: знак и величина вертикального остатка на синтетической ветви.

Синтетическая линия: C_ref(G) = 0.70 + 0.20 × (G − 11.5) для 11.5 ≤ G < 13.5.
Генерируем 20 точек: G внутри интервала, цвет = C_ref(G) + s,
где s ∈ {−0.05, 0, +0.05} mag.

Проверяем: |r − s| <= 1e-8, где r = C_obs − C_ref(G).

Опция --flip-sign: обращает знак остатка (имитация ошибки).
"""
import argparse
import math
import sys
from fractions import Fraction

TOL = Fraction(1, 10**8)  # 1e-8 mag


# Синтетическая опорная линия: C_ref(G) = 0.70 + 0.20 × (G − 11.5)
# Используем целые аргументы для совместимости с Python 3.9
def reference_colour(g):
    return Fraction(70, 100) + Fraction(20, 100) * (g - Fraction(115, 10))


def residual(observed_colour, g, colour_at, flip_sign=False):
    """r = C_obs − C_ref(G). Если flip_sign, r = -(C_obs − C_ref(G))."""
    r = observed_colour - colour_at(g)
    if flip_sign:
        r = -r
    return r


def make_test_points():
    """Генерирует (G, s, ожидаемый r) для 20 тестовых точек."""
    points = []
    # shifts as strings to avoid TypeError in Python <3.11 Fraction(str)
    shift_strs = ["-5/100", "0", "5/100"]
    shifts = [Fraction(s) for s in shift_strs]

    # Точки G внутри [11.5, 13.5)
    g_strs = [
        "115/10",     # 11.5
        "1176/100",   # 11.76
        "1202/100",   # 12.02
        "1228/100",   # 12.28
        "1254/100",   # 12.54
        "128/10",     # 12.80
        "135/10",     # 13.5
        "1163/100",   # 11.63
        "119/10",     # 11.90
        "1215/100",   # 12.15
    ]
    all_g = [Fraction(s) for s in g_strs]

    # Каждый сдвиг для нескольких G
    for g in all_g[:3]:
        for s in shifts:
            points.append((g, s))
    for g in all_g[3:6]:
        for s in shifts[:2]:
            points.append((g, s))
    # Добавим остальные до 20
    idx = 0
    while len(points) < 20:
        g = all_g[6 + idx % (len(all_g) - 6)]
        s = shifts[idx % len(shifts)]
        points.append((g, s))
        idx += 1

    if len(points) > 20:
        points = points[:20]
    return points


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flip-sign", action="store_true",
                        help="Обратить знак остатка (имитация ошибки)")
    args = parser.parse_args()

    colour_at = reference_colour
    test_points = make_test_points()
    n = len(test_points)
    passed = 0
    failed = 0

    print(f"Синтетическая проверка: {n} точек{' (знак обращён)' if args.flip_sign else ''}")
    print(f"{'G':>8} {'C_ref':>8} {'сдвиг s':>8} {'r':>8}  {'|r-s|':>8}  статус")

    for g, s in test_points:
        c_ref = colour_at(g)
        obs = c_ref + s
        r = residual(obs, g, colour_at, flip_sign=args.flip_sign)
        err = abs(r - s)
        status = "OK" if err <= TOL else "FAIL"
        print(f"{float(g):8.3f} {float(c_ref):8.3f} {float(s):+8.2f} {float(r):+8.2f} {float(err):8.2e}  {status}")
        if status == "OK":
            passed += 1
        else:
            failed += 1

    print()
    if failed == 0:
        print(f"РЕЗУЛЬТАТ: PASS — все {n}/{n} точек прошли (max |r-s| <= 1e-8 mag)")
        sys.exit(0)
    else:
        print(f"РЕЗУЛЬТАТ: FAIL — {failed}/{n} точек не прошли")
        sys.exit(1)


if __name__ == "__main__":
    main()