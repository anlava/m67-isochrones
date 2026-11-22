#!/usr/bin/env python3
"""Independent verification of synthetic_test.py components."""
import csv, math, random, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AGE_DIR = ROOT / "experiment" / "age_fit"

# ── Correct Kroupa sampler using inverse CDF ──────────────────

# Kroupa 2001: xi(m) = A * m^-1.3 for 0.08<=m<0.5, B * m^-2.3 for 0.5<=m<=3.0
# Continuity: B = A * 0.5^(-1.3+2.3) = A * 0.5
# Power-law CDF: F(m|a,b,alpha) = [m^(1-alpha)-a^(1-alpha)]/[b^(1-alpha)-a^(1-alpha)]
# Inverse: m = [F*(b^(1-alpha)-a^(1-alpha)) + a^(1-alpha)]^(1/(1-alpha))

LO_MASS, MID_MASS, HI_MASS = 0.08, 0.5, 3.0
ALPHA_LO, ALPHA_HI = 1.3, 2.3

# Precompute normalization factors
def pow_int(a, b, alpha):
    """Integral of m^(-alpha) dm from a to b."""
    if alpha == 1.0:
        return math.log(b / a)
    return (b ** (1 - alpha) - a ** (1 - alpha)) / (1 - alpha)

norm_lo = pow_int(LO_MASS, MID_MASS, ALPHA_LO)
norm_hi = pow_int(MID_MASS, HI_MASS, ALPHA_HI) * MID_MASS  # continuity factor 0.5
total_norm = norm_lo + norm_hi
p_seg_low = norm_lo / total_norm  # probability of sampling from low-mass segment

def kroupa_mass_correct(rng):
    """Sample from Kroupa (2001) IMF using inverse CDF."""
    if rng.random() < p_seg_low:
        # Sample from m^(-1.3) on [LO, MID]
        alpha = ALPHA_LO
        a, b = LO_MASS, MID_MASS
    else:
        # Sample from m^(-2.3) on [MID, HI] (with continuity factor)
        alpha = ALPHA_HI
        a, b = MID_MASS, HI_MASS
    u = rng.random()
    # Inverse CDF: m = [u*(b^(1-alpha)-a^(1-alpha)) + a^(1-alpha)]^(1/(1-alpha))
    pow1 = 1 - alpha
    m = (u * (b ** pow1 - a ** pow1) + a ** pow1) ** (1.0 / pow1)
    return m

def test_kroupa():
    rng = random.Random(42)
    N = 100000
    masses = [kroupa_mass_correct(rng) for _ in range(N)]

    # Theoretical fractions
    F_lo = lambda x: x ** (1 - ALPHA_LO) / (1 - ALPHA_LO)
    F_hi = lambda x: x ** (1 - ALPHA_HI) / (1 - ALPHA_HI)

    bins = [(0.08, 0.15), (0.15, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 3.0)]
    print("Test 1: Kroupa IMF inverse-CDF (N=%d)" % N)
    print("%-20s %12s %12s %7s" % ("Range (Msun)", "Observed %", "Expected %", "Dev"))
    ok = True
    for lo, hi in bins:
        obs = sum(1 for m in masses if lo <= m < hi) / N * 100
        if hi <= MID_MASS:
            exp = (F_lo(hi) - F_lo(lo)) / total_norm * 100
        elif lo >= MID_MASS:
            exp = (F_hi(hi) - F_hi(lo)) * MID_MASS / total_norm * 100
        else:
            exp = ((F_lo(MID_MASS) - F_lo(lo)) + (F_hi(hi) - F_hi(MID_MASS)) * MID_MASS) / total_norm * 100
        dev = abs(obs - exp)
        flag = " <--" if dev > 1.5 else ""
        if dev > 1.5:
            ok = False
        print("%.2f-%.2f Msun       %10.2f%%  %10.2f%%  %6.2f%%%s" % (lo, hi, obs, exp, dev, flag))
    print("  Kroupa: %s" % ("PASS" if ok else "FAIL"))
    return ok

# ── Distance to segment ────────────────────────────────────────
def point_to_segment(px, py, x1, y1, x2, y2):
    """Min Euclidean distance from point (px,py) to line segment (x1,y1)-(x2,y2)."""
    dx, dy = x2 - x1, y2 - y1
    if dx == 0 and dy == 0:
        return math.hypot(px - x1, py - y1)
    t = ((px - x1) * dx + (py - y1) * dy) / (dx*dx + dy*dy)
    if t < 0:
        return math.hypot(px - x1, py - y1)
    elif t > 1:
        return math.hypot(px - x2, py - y2)
    else:
        ix, iy = x1 + t * dx, y1 + t * dy
        return math.hypot(px - ix, py - iy)

def test_distance():
    """Point exactly midway between two vertices should have d = 0."""
    # Segment from (12, 0.8) to (13, 1.0). Midpoint is (12.5, 0.9).
    d = point_to_segment(12.5, 0.9, 12, 0.8, 13, 1.0)
    ok1 = abs(d) < 1e-12
    print("\nTest 2: Distance to segment")
    print("  Midpoint (12.5,0.9) on segment: d=%.2e %s" % (d, "PASS" if ok1 else "FAIL"))
    # Point at (12.5, 1.0): projection gives t=0.519, dist ≈ 0.0978
    d2 = point_to_segment(12.5, 1.0, 12, 0.8, 13, 1.0)
    expected = math.hypot(0.5 - 0.519*1, 1.0 - (0.8 + 0.519*0.2))
    ok2 = abs(d2 - expected) < 1e-6
    print("  Point (12.5,1.0) off segment: d=%.10f (expected %.10f) %s" % (d2, expected, "PASS" if ok2 else "FAIL"))
    # Point on vertex
    d3 = point_to_segment(12, 0.8, 12, 0.8, 13, 1.0)
    ok3 = abs(d3) < 1e-12
    print("  At vertex (12,0.8): d=%.2e %s" % (d3, "PASS" if ok3 else "FAIL"))
    return ok1 and ok2 and ok3

if __name__ == "__main__":
    ok = True
    ok &= test_kroupa()
    ok &= test_distance()
    print("\n%s" % ("All PASS" if ok else "SOME FAILED"))
    sys.exit(0 if ok else 1)