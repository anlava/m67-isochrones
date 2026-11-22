#!/usr/bin/env python3
"""Plot Figure 2 (age fit) from m67_fit_results.csv.

Usage:
    python3 experiment/age_fit/plot_age_fit.py
    # Output: experiment/age_fit/m67_cmd_fit.png
    # Copy to manuscript: cp experiment/age_fit/m67_cmd_fit.png manuscript/figures/figure3_age_fit.png
"""
import csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AGE_DIR = ROOT / "experiment" / "age_fit"
CSV_PATH = AGE_DIR / "m67_fit_results.csv"
OUT_PNG = AGE_DIR / "m67_cmd_fit.png"

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib is required. Install with: pip install matplotlib")
    sys.exit(1)

def main():
    with open(CSV_PATH) as f:
        rows = list(csv.DictReader(f))

    fig, ax = plt.subplots(figsize=(8, 5))

    for grid, color, marker in [('Cavallo', '#1f77b4', 'o'), ('Bossini', '#d62728', 's')]:
        ages = sorted(set(float(r['logAge']) for r in rows if r['grid'] == grid))
        vals = [float([r for r in rows if r['grid'] == grid and r['logAge'] == str(la)][0]['median_dist'])
                for la in ages]
        ax.plot(ages, vals, color=color, marker=marker, linewidth=1.5,
                markersize=6, label=f'{grid} PARSEC 1.2S')
        best_la = ages[vals.index(min(vals))]
        best_m = min(vals)
        ax.annotate(f'logAge={best_la:.2f}\nM={best_m:.4f}', xy=(best_la, best_m),
                    xytext=(best_la + 0.10, best_m - 0.005),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1),
                    fontsize=9, color=color, zorder=5)

    ax.set_xlabel('logAge (dex)', fontsize=12)
    ax.set_ylabel('Median distance M (mag)', fontsize=12)
    ax.set_title('M67 turnoff fit $-$ median CMD distance vs age', fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(9.25, 9.75)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Saved: {OUT_PNG}  ({Path(OUT_PNG).stat().st_size} bytes)")

if __name__ == "__main__":
    main()
