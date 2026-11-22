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
    """
    Сверяет полные наборы из 257 уникальных GaiaDR3 и BP-RP между TSV и CSV.
    Выборка из CSV делается тем же фильтром 11.5 ≤ G < 13.5.
    """
    # Read TSV — 257 stars with 11.5 ≤ G < 13.5
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
                gaia = cols[1].strip()
                if gaia in tsv:
                    print(f"    DUPLICATE GaiaDR3 in TSV: {gaia}")
                    return False
                tsv[gaia] = {"g": g, "bp_rp": bp - rp}

    if len(tsv) != 257:
        print(f"✗  TSV: expected 257 unique GaiaDR3, got {len(tsv)}")
        return False

    # Read CSV — select same filter 11.5 ≤ G < 13.5
    sel = {}
    with open(SEL_CSV) as f:
        for r in csv.DictReader(f):
            g = float(r["Gmag"].strip())
            if 11.5 <= g < 13.5:
                gaia = r["GaiaDR3"].strip()
                if gaia in sel:
                    print(f"    DUPLICATE GaiaDR3 in CSV: {gaia}")
                    return False
                sel[gaia] = float(r["BP_RP"])

    if len(sel) != 257:
        print(f"✗  CSV: expected 257 unique GaiaDR3, got {len(sel)}")
        return False

    # Same IDs?
    tsv_ids = set(tsv.keys())
    sel_ids = set(sel.keys())
    if tsv_ids != sel_ids:
        missing = tsv_ids - sel_ids
        extra = sel_ids - tsv_ids
        print(f"✗  ID MISMATCH: missing={len(missing)}, extra={len(extra)}")
        if missing:
            print(f"    Missing from CSV: {list(missing)[:3]}")
        if extra:
            print(f"    Extra in CSV: {list(extra)[:3]}")
        return False

    # BP-RP values match?
    errors = []
    for gaia in tsv:
        t = tsv[gaia]["bp_rp"]
        s = sel[gaia]
        if abs(t - s) > 1e-6:
            errors.append((gaia, t, s))

    if errors:
        print(f"✗  BP-RP mismatches: {len(errors)}/257")
        for g, t, s in errors[:5]:
            print(f"    {g}: TSV={t:.6f} CSV={s:.6f}")
        return False
    print(f"✓  Colour cross-check: 257/257 GaiaDR3 IDs and BP-RP match")
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