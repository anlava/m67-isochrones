#!/usr/bin/env python3
"""Fit M67 age using the 257 turnoff stars against both PARSEC grids."""
import csv, math, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AGE_DIR = ROOT / "experiment" / "age_fit"

DATA_TSV = ROOT / "data" / "raw" / "m67_members_hr2023.tsv"
GRIDS = {
    "Cavallo": AGE_DIR / "grid_cavallo2024_av018.csv",
    "Bossini": AGE_DIR / "grid_bossini2019_av0115.csv",
}
DMOBS = {"Cavallo": 9.64, "Bossini": 9.726}

def load_257_stars():
    """Load exactly 257 stars with 11.5 ≤ G < 13.5 from the Hunt & Reffert 2023 TSV."""
    stars = []
    with open(DATA_TSV, encoding="utf-8-sig") as f:
        for line in f:
            if line.startswith("Name\t") and "GaiaDR3" in line:
                break
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 7:
                continue
            try:
                g = float(cols[4].strip())
                bp_rp = float(cols[6].strip())
            except ValueError:
                continue
            if 11.5 <= g < 13.5:
                stars.append((g, bp_rp))
    assert len(stars) == 257, f"Expected 257 stars, got {len(stars)}"
    return stars

def load_grid(csv_path):
    """Return {logAge: [(Mini, G, BP_RP), ...]} for labels 0,1,2."""
    ages = {}
    with open(csv_path, encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            if int(rec["label"]) not in (0, 1, 2):
                continue
            la = round(float(rec["logAge"]), 2)
            if la not in ages:
                ages[la] = []
            ages[la].append((float(rec["Mini"]), float(rec["G_app"]), float(rec["BP_RP"])))
    for la in ages:
        ages[la].sort(key=lambda x: x[0])
    return ages

def distance_to_isochrone(g_star, c_star, iso):
    """Min Euclidean distance from star to isochrone polyline segment."""
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

def main():
    stars = load_257_stars()
    print(f"Loaded {len(stars)} M67 turnoff stars.\n")

    results = []  # (grid, logAge, median_dist)
    best_per_grid = {}

    for grid_name in ["Cavallo", "Bossini"]:
        grid = load_grid(GRIDS[grid_name])
        ages = sorted(grid.keys())
        print(f"  {grid_name} grid")
        grid_results = []

        for la in ages:
            iso = grid[la]
            dists = [distance_to_isochrone(g, c, iso) for g, c in stars]
            dists.sort()
            m = dists[len(dists) // 2]
            grid_results.append((la, m))
            results.append({"grid": grid_name, "logAge": la, "median_dist": m})

        # Find best
        best_la, best_m = min(grid_results, key=lambda x: x[1])
        best_idx = ages.index(best_la)
        age_gyr = 10 ** (best_la - 9)

        print(f"    Best: logAge={best_la:.2f} ({age_gyr:.1f} Gyr), M={best_m:.6f} mag")
        if best_idx == 0:
            print("    ⚠ BOUNDARY: minimum at youngest node")
        elif best_idx == len(ages) - 1:
            print("    ⚠ BOUNDARY: minimum at oldest node")
        if best_idx > 0:
            prev_la, prev_m = grid_results[best_idx - 1]
            print(f"    Neighbor -0.05 dex: logAge={prev_la:.2f}, M={prev_m:.6f}")
        if best_idx < len(ages) - 1:
            next_la, next_m = grid_results[best_idx + 1]
            print(f"    Neighbor +0.05 dex: logAge={next_la:.2f}, M={next_m:.6f}")
        print()

        best_per_grid[grid_name] = (best_la, best_m, age_gyr)

    # Save CSV
    out_csv = AGE_DIR / "m67_fit_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["grid", "logAge", "median_dist"])
        w.writeheader()
        w.writerows(results)
    print(f"Saved: {out_csv}")

    # Summary table
    print("\n" + "=" * 56)
    print("FIT RESULTS SUMMARY")
    print("=" * 56)
    print(f"{'Grid':<12} {'Best logAge':>14} {'Age (Gyr)':>12} {'M (mag)':>12}")
    print("-" * 50)
    for g in ["Cavallo", "Bossini"]:
        la, m, gyr = best_per_grid[g]
        print(f"{g:<12} {la:>14.2f} {gyr:>12.2f} {m:>12.6f}")
    print("-" * 50)

    return 0

if __name__ == "__main__":
    sys.exit(main())