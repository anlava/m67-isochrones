#!/usr/bin/env python3
"""Fit M67 age using the 257 turnoff stars against both PARSEC grids.

Uses distance_to_isochrone from synthetic_test.py (отрезки ломаной, label=0,1,2).
Сверяет BP-RP двух источников: исходный TSV и results/m67_members_selected.csv.
"""
import csv, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiment" / "age_fit"))

from synthetic_test import distance_to_isochrone, load_grid

DATA_TSV = ROOT / "data" / "raw" / "m67_members_hr2023.tsv"
SEL_CSV  = ROOT / "results" / "m67_members_selected.csv"
AGE_DIR  = ROOT / "experiment" / "age_fit"

GRIDS = {
    "Cavallo": AGE_DIR / "grid_cavallo2024_av018.csv",
    "Bossini": AGE_DIR / "grid_bossini2019_av0115.csv",
}
DMOBS = {"Cavallo": 9.64, "Bossini": 9.726}

# ── Cross-check colours ────────────────────────────────────────
def crosscheck_colours():
    """Сверяет BP-RP исходного TSV и selected CSV по GaiaDR3 для 257 звёзд."""
    errors = []

    # Read TSV
    tsv = {}
    with open(DATA_TSV, encoding="utf-8-sig") as f:
        for line in f:
            if line.startswith("Name\t"):
                break
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 7 or cols[0] == "--------------------":
                continue
            try:
                g = float(cols[4].strip())
                bp = float(cols[5].strip())
                rp = float(cols[6].strip())
            except ValueError:
                continue
            if 11.5 <= g < 13.5:
                tsv[cols[1].strip()] = {"g": g, "bp_rp": bp - rp, "bp": bp, "rp": rp}

    # Read selected CSV
    sel = {}
    with open(SEL_CSV) as f:
        for r in csv.DictReader(f):
            sel[r["GaiaDR3"].strip()] = float(r["BP_RP"])

    common = set(tsv.keys()) & set(sel.keys())
    for gaia in common:
        t = tsv[gaia]["bp_rp"]
        s = sel[gaia]
        if abs(t - s) > 1e-6:
            errors.append((gaia, t, s))

    if errors:
        print(f"⚠  BP-RP mismatches: {len(errors)}/{len(common)}")
        for g, t, s in errors[:5]:
            print(f"    {g}: TSV={t:.6f} CSV={s:.6f}")
        return False
    else:
        print(f"✓  BP-RP cross-check: all {len(common)} stars match")
        return True

# ── Load 257 stars ─────────────────────────────────────────────
def load_257_stars():
    """Читает 257 звёзд, возвращает [(G_app, BP_RP), ...]."""
    stars = []
    with open(DATA_TSV, encoding="utf-8-sig") as f:
        for line in f:
            if line.startswith("Name\t"):
                break
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 7 or cols[0] == "--------------------":
                continue
            try:
                g = float(cols[4].strip())
                bp = float(cols[5].strip())
                rp = float(cols[6].strip())
            except ValueError:
                continue
            if 11.5 <= g < 13.5:
                stars.append((g, bp - rp, cols[1].strip()))
    assert len(stars) == 257, f"Expected 257 stars, got {len(stars)}"
    return stars

# ── Main ───────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("M67 turnoff fit — 257 stars, 2 grids")
    print("=" * 60)

    # 1. Cross-check
    ok = crosscheck_colours()
    if not ok:
        print("ABORT: colour mismatch between TSV and CSV")
        return 1

    # 2. Load stars
    stars = load_257_stars()
    print(f"\nStars: {len(stars)} (G 11.5–13.5)")

    # 3. Fit
    results = []
    best_per_grid = {}
    for gname in ["Cavallo", "Bossini"]:
        grid = load_grid(GRIDS[gname])
        ages = sorted(grid.keys())
        print(f"\n{gname}:")
        grid_res = []
        for la in ages:
            iso = grid[la]
            dists = [distance_to_isochrone(s[0], s[1], iso) for s in stars]
            dists.sort()
            m = dists[len(dists) // 2]
            grid_res.append((la, m))
            results.append({"grid": gname, "logAge": la, "median_dist": m})

        best_la, best_m = min(grid_res, key=lambda x: x[1])
        best_idx = ages.index(best_la)
        age_gyr = 10 ** (best_la - 9)
        print(f"  Best logAge={best_la:.2f}  ({age_gyr:.2f} Gyr)  M={best_m:.6f}")
        if best_idx == 0:
            print("  ⚠  BOUNDARY: minimum at youngest node")
        elif best_idx == len(ages) - 1:
            print("  ⚠  BOUNDARY: minimum at oldest node")
        else:
            print(f"  Neighbor -0.05: logAge={grid_res[best_idx-1][0]:.2f}, M={grid_res[best_idx-1][1]:.6f}")
            print(f"  Neighbor +0.05: logAge={grid_res[best_idx+1][0]:.2f}, M={grid_res[best_idx+1][1]:.6f}")
        best_per_grid[gname] = (best_la, best_m, age_gyr)

    # 4. Save CSV
    out_csv = AGE_DIR / "m67_fit_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["grid", "logAge", "median_dist"])
        w.writeheader()
        w.writerows(results)
    print(f"\nSaved: {out_csv}")

    # 5. Table
    print(f"\n{'Grid':<12} {'Best logAge':>14} {'Age (Gyr)':>12} {'M (mag)':>12}")
    print("-" * 50)
    for g in ["Cavallo", "Bossini"]:
        la, m, gyr = best_per_grid[g]
        flag = " ⚠ boundary" if (la == 9.30 or la == 9.70) else ""
        print(f"{g:<12} {la:>14.2f} {gyr:>12.2f} {m:>12.6f}{flag}")
    print("-" * 50)

    return 0

if __name__ == "__main__":
    sys.exit(main())