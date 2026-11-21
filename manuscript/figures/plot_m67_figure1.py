import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
MEMBERS_CSV = os.path.join(BASE_DIR, "results", "m67_members_selected.csv")
ISO_CAV_CSV = os.path.join(BASE_DIR, "results", "m67_isochrone_apparent.csv")
ISO_BOS_CSV = os.path.join(BASE_DIR, "results", "m67_isochrone_bossini2019_apparent.csv")
RESIDUALS_CSV = os.path.join(BASE_DIR, "experiment", "results", "m67_residuals_inrt_comparison.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "manuscript", "figures")

def load_data():
    members = pd.read_csv(MEMBERS_CSV)
    iso_cav = pd.read_csv(ISO_CAV_CSV)
    iso_bos = pd.read_csv(ISO_BOS_CSV)
    residuals = pd.read_csv(RESIDUALS_CSV)
    return members, iso_cav, iso_bos, residuals

def split_isochrone_by_label(iso_df):
    """Split isochrone into continuous blocks of same label."""
    if iso_df.empty:
        return []
    
    blocks = []
    current_label = iso_df.iloc[0]["label"]
    current_block = []
    
    for _, row in iso_df.iterrows():
        if row["label"] != current_label:
            blocks.append((current_label, current_block))
            current_label = row["label"]
            current_block = [row]
        else:
            current_block.append(row)
            
    if current_block:
        blocks.append((current_label, current_block))
        
    return blocks

def plot_panel_a(ax, members, iso_cav, iso_bos):
    # Observed stars
    obs_bp_rp = members["BP_RP"].values
    obs_g = members["Gmag"].values
    
    ax.scatter(obs_bp_rp, obs_g, s=5, color="grey", label="M67 members", zorder=3)
    
    # Shaded region 15 <= G <= 18
    g_min, g_max = obs_g.min(), obs_g.max()
    g_margin = (g_max - g_min) * 0.03
    bp_rp_min, bp_rp_max = obs_bp_rp.min(), obs_bp_rp.max()
    bp_rp_margin = (bp_rp_max - bp_rp_min) * 0.03
    
    ax.axhspan(15, 18, color="lightgray", alpha=0.3, zorder=1)
    
    # Cavallo (red)
    cav_blocks = split_isochrone_by_label(iso_cav)
    for i, (label, block) in enumerate(cav_blocks):
        if label == 9 and len(block) > 0 and block[0].get("logL") == -9.999:
            continue
        bp_rp = np.array([pt["BP_RP"] for pt in block])
        g_app = np.array([pt["G_app"] for pt in block])
        lbl = "Cavallo et al. (2024)" if i == 0 else None
        ax.plot(bp_rp, g_app, c="red", linewidth=1.5, label=lbl, zorder=2)
        
    # Bossini (blue)
    bos_blocks = split_isochrone_by_label(iso_bos)
    for i, (label, block) in enumerate(bos_blocks):
        if label == 9 and len(block) > 0 and block[0].get("logL") == -9.999:
            continue
        bp_rp = np.array([pt["BP_RP"] for pt in block])
        g_app = np.array([pt["G_app"] for pt in block])
        lbl = "Bossini et al. (2019)" if i == 0 else None
        ax.plot(bp_rp, g_app, c="blue", linewidth=1.5, label=lbl, zorder=2)
    
    # Axes setup
    ax.invert_yaxis()
    ax.set_xlabel("G$_{\\mathrm{BP}}$ − G$_{\\mathrm{RP}}$ (mag)")
    ax.set_ylabel("G (mag)")
    ax.set_xlim(bp_rp_min - bp_rp_margin, bp_rp_max + bp_rp_margin)
    ax.set_ylim(g_max + g_margin, g_min - g_margin)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

def plot_panel_b(ax, residuals):
    samples = ["base", "inrt_true"]
    labels_map = {"base": "Full sample", "inrt_true": "Within tidal radius"}
    markers = ["o", "s"]
    colors = ["red", "blue"]
    
    bins = ["[15,16)", "[16,17)", "[17,18]"]
    y_positions = {bin_val: 2.5 - i * 0.8 for i, bin_val in enumerate(bins)}
    
    # Смещение по y для двух выборок
    y_offsets = {"base": 0.10, "inrt_true": -0.10}
    
    for i, sample in enumerate(samples):
        sub = residuals[residuals["sample"] == sample]
        D_vals = sub["D"].tolist()
        ci_low = sub["ci_low"].tolist()
        ci_high = sub["ci_high"].tolist()
        
        lower_mag = [D - cl for D, cl in zip(D_vals, ci_low)]
        upper_mag = [ch - D for D, ch in zip(D_vals, ci_high)]
        
        # Смещаем y координаты для маркеров и интервалов
        y_coords = [y_positions[b] + y_offsets[sample] for b in sub["g_bin"]]
        
        ax.errorbar(D_vals, y_coords, xerr=np.array([lower_mag, upper_mag]), 
                    fmt=markers[i], color=colors[i], capsize=3, ecolor='black',
                    label=labels_map[sample], zorder=3)
        
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("D (mag)")
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(bins)
    ax.set_ylim(-0.5, 3.5)
    ax.legend(loc="upper left")

def main():
    members, iso_cav, iso_bos, residuals = load_data()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    plot_panel_a(ax1, members, iso_cav, iso_bos)
    plot_panel_b(ax2, residuals)
    
    # Убираем общий заголовок
    # fig.suptitle("Figure 1", fontsize=14, fontweight="bold")
    
    # Добавляем метки панелей (a) и (b)
    ax1.text(0.5, 1.05, "(a)", transform=ax1.transAxes, fontsize=14, fontweight="bold", ha="center")
    ax2.text(0.5, 1.05, "(b)", transform=ax2.transAxes, fontsize=14, fontweight="bold", ha="center")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    
    pdf_path = os.path.join(OUTPUT_DIR, "figure1.pdf")
    png_path = os.path.join(OUTPUT_DIR, "figure1.png")
    fig.savefig(pdf_path, format="pdf")
    fig.savefig(png_path, format="png", dpi=300)
    print(f"Saved {pdf_path} and {png_path}")

if __name__ == "__main__":
    main()