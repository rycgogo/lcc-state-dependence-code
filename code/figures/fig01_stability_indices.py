# ==========================================================
# Figure 1
# This script plots the relationships between EIS and two
# stability metrics, LTS and NSS, over global oceans.
#
# Panels:
# (a) EIS-LTS in polar regions
# (b) EIS-LTS in mid-low latitude regions
# (c) EIS-NSS in polar regions
# (d) EIS-NSS in mid-low latitude regions
#
# Modifications:
# 1. Panel labels are placed at the beginning of subplot titles.
# 2. In-panel label boxes are removed.
# 3. Colorbar thickness is unified across all panels.
# ==========================================================

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
from scipy.stats import linregress
import os
import pandas as pd
from pathlib import Path

# ==========================================================
# 1. Global settings
# ==========================================================
project_dir = Path(__file__).resolve().parents[2]
data_file = project_dir / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc"
output_dir = project_dir / "outputs"
data_save_dir = os.path.join(output_dir, "Sensitivity_Data")

output_filename = "Figure1_Stability_Scatter_Final.png"

summary_csv_filename = "Figure1_Stability_Scatter_Stats.csv"
detailed_csv_filename = "Figure1_Detailed_Bin_Data.csv"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(data_save_dir):
    os.makedirs(data_save_dir)

# --- Font settings ---
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

title_fontsize = 24
label_fontsize = 26
tick_fontsize = 24
cbar_tick_fontsize = 18
font_weight = "bold"

# unified colorbar style
cbar_fraction = 0.046
cbar_pad = 0.04

# ==========================================================
# 2. Data loading and processing
# ==========================================================
print("Loading data...")

try:
    ds = xr.open_dataset(data_file)

    # Coordinate processing
    rename_dict = {}

    if "latitude" in ds.coords:
        rename_dict["latitude"] = "lat"

    if "longitude" in ds.coords:
        rename_dict["longitude"] = "lon"

    if rename_dict:
        ds = ds.rename(rename_dict)

    # The archived SIC field is used as the ocean/valid-data mask, so no
    # separate land-sea-mask file is required.
    ocean_mask = ds["sic"].mean(dim="time").notnull()

    def get_region_data(ds, mask_array, lat_condition):
        ds_ocean = ds.where(mask_array)
        lat_mask = lat_condition(ds_ocean["lat"])
        ds_region = ds_ocean.where(lat_mask)

        eis = ds_region["eis"].values.flatten()
        lts = ds_region["lts"].values.flatten()
        nss = ds_region["nss"].values.flatten()

        mask_valid = ~np.isnan(eis) & ~np.isnan(lts) & ~np.isnan(nss)

        return {
            "eis": eis[mask_valid],
            "lts": lts[mask_valid],
            "nss": nss[mask_valid],
        }

    print("Extracting data...")

    cond_polar = lambda lat: np.abs(lat) >= 60
    cond_midlow = lambda lat: np.abs(lat) < 60

    data_polar = get_region_data(ds, ocean_mask, cond_polar)
    data_midlow = get_region_data(ds, ocean_mask, cond_midlow)

    all_vals = np.concatenate([
        data_polar["eis"],
        data_polar["lts"],
        data_polar["nss"],
        data_midlow["eis"],
        data_midlow["lts"],
        data_midlow["nss"],
    ])

    global_min = np.percentile(all_vals, 0.5)
    global_max = np.percentile(all_vals, 99.5)

    span = global_max - global_min

    plot_lim_min = global_min - 0.05 * span
    plot_lim_max = global_max + 0.05 * span

    # ======================================================
    # 3. Density calculations
    # ======================================================
    bins = 80
    range_lims = [
        [plot_lim_min, plot_lim_max],
        [plot_lim_min, plot_lim_max],
    ]

    hist_p1, _, _ = np.histogram2d(
        data_polar["lts"],
        data_polar["eis"],
        bins=bins,
        range=range_lims,
    )

    hist_p2, _, _ = np.histogram2d(
        data_polar["nss"],
        data_polar["eis"],
        bins=bins,
        range=range_lims,
    )

    max_polar = max(hist_p1.max(), hist_p2.max())

    hist_m1, _, _ = np.histogram2d(
        data_midlow["lts"],
        data_midlow["eis"],
        bins=bins,
        range=range_lims,
    )

    hist_m2, _, _ = np.histogram2d(
        data_midlow["nss"],
        data_midlow["eis"],
        bins=bins,
        range=range_lims,
    )

    max_midlow = max(hist_m1.max(), hist_m2.max())

    stats_list = []
    detailed_list = []

    # ======================================================
    # 4. Plotting function
    # ======================================================
    def plot_density_scatter_fit(
        ax,
        x,
        y,
        xlabel,
        ylabel,
        title,
        label_text,
        vmax_val,
        region_name,
        pair_name,
    ):
        counts, xedges, yedges, im = ax.hist2d(
            x,
            y,
            bins=bins,
            range=range_lims,
            cmap="jet",
            norm=LogNorm(vmin=1, vmax=vmax_val),
            cmin=1,
        )

        nx, ny = counts.shape
        subplot_label = label_text.replace("(", "").replace(")", "")

        for i in range(nx):
            for j in range(ny):
                if counts[i, j] > 0:
                    detailed_list.append({
                        "Subplot": subplot_label,
                        "Region": region_name,
                        "X_Var": pair_name.split("-")[1],
                        "Y_Var": "EIS",
                        "X_Bin_Center": round((xedges[i] + xedges[i + 1]) / 2, 3),
                        "Y_Bin_Center": round((yedges[j] + yedges[j + 1]) / 2, 3),
                        "Count": int(counts[i, j]),
                    })

        stats_entry = {
            "Region": region_name,
            "Pair": pair_name,
            "Label": subplot_label,
            "Sample_Size": len(x),
        }

        if len(x) > 0:
            slope, intercept, r_value, p_value, std_err = linregress(x, y)

            stats_entry.update({
                "Slope": slope,
                "Intercept": intercept,
                "R_squared": r_value ** 2,
                "P_value": p_value,
            })

            line_x = np.array([plot_lim_min, plot_lim_max])
            ax.plot(line_x, slope * line_x + intercept, "r-", linewidth=2.5)

            stats_text = (
                f"$y = {slope:.2f}x {'+' if intercept >= 0 else '-'} {abs(intercept):.2f}$\n"
                f"$R^2 = {r_value ** 2:.2f}$"
            )

            ax.text(
                0.95,
                0.05,
                stats_text,
                transform=ax.transAxes,
                fontsize=tick_fontsize - 4,
                fontweight="bold",
                va="bottom",
                ha="right",
                bbox=dict(
                    facecolor="white",
                    alpha=0.8,
                    edgecolor="red",
                    boxstyle="round,pad=0.3",
                ),
            )

        stats_list.append(stats_entry)

        ax.set_xlim(plot_lim_min, plot_lim_max)
        ax.set_ylim(plot_lim_min, plot_lim_max)
        ax.set_aspect("equal", adjustable="box")

        ax.set_xlabel(
            xlabel,
            fontsize=label_fontsize,
            fontweight=font_weight,
            fontfamily="Times New Roman",
        )

        ax.set_ylabel(
            ylabel,
            fontsize=label_fontsize,
            fontweight=font_weight,
            fontfamily="Times New Roman",
        )

        # label moved to title
        ax.set_title(
            f"{label_text} {title}",
            fontsize=title_fontsize,
            fontweight=font_weight,
            fontfamily="Times New Roman",
            pad=10,
        )

        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(5))

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=tick_fontsize,
            width=2,
            length=6,
        )

        for t in ax.get_xticklabels() + ax.get_yticklabels():
            t.set_fontweight(font_weight)
            t.set_fontfamily("Times New Roman")

        ax.grid(True, linestyle="--", alpha=0.3)

        return im

    # ======================================================
    # 5. Execution
    # ======================================================
    print("Plotting Figure 1...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))

    im1 = plot_density_scatter_fit(
        axes[0, 0],
        data_polar["lts"],
        data_polar["eis"],
        "LTS (K)",
        "EIS (K)",
        "EIS-LTS",
        "(a)",
        max_polar,
        "Polar",
        "EIS-LTS",
    )

    im2 = plot_density_scatter_fit(
        axes[0, 1],
        data_midlow["lts"],
        data_midlow["eis"],
        "LTS (K)",
        "EIS (K)",
        "EIS-LTS",
        "(b)",
        max_midlow,
        "Mid-Low Lat",
        "EIS-LTS",
    )

    im3 = plot_density_scatter_fit(
        axes[1, 0],
        data_polar["nss"],
        data_polar["eis"],
        "NSS (K)",
        "EIS (K)",
        "EIS-NSS",
        "(c)",
        max_polar,
        "Polar",
        "EIS-NSS",
    )

    im4 = plot_density_scatter_fit(
        axes[1, 1],
        data_midlow["nss"],
        data_midlow["eis"],
        "NSS (K)",
        "EIS (K)",
        "EIS-NSS",
        "(d)",
        max_midlow,
        "Mid-Low Lat",
        "EIS-NSS",
    )

    axes[0, 0].annotate(
        "Polar",
        xy=(0.5, 1.12),
        xycoords="axes fraction",
        fontsize=26,
        fontweight="bold",
        ha="center",
    )

    axes[0, 1].annotate(
        "Mid-Low Lat",
        xy=(0.5, 1.12),
        xycoords="axes fraction",
        fontsize=26,
        fontweight="bold",
        ha="center",
    )

    # unified colorbar thickness
    for ax, im in zip(axes.flat, [im1, im2, im3, im4]):
        cbar = plt.colorbar(
            im,
            ax=ax,
            fraction=cbar_fraction,
            pad=cbar_pad,
        )

        cbar.ax.tick_params(labelsize=cbar_tick_fontsize)

        for t in cbar.ax.get_yticklabels():
            t.set_fontweight(font_weight)
            t.set_fontfamily("Times New Roman")

    plt.subplots_adjust(
        left=0.08,
        right=0.98,
        top=0.88,
        bottom=0.08,
        hspace=0.25,
        wspace=0.25,
    )

    # Save figure
    save_path = os.path.join(output_dir, output_filename)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"[Success] Figure 1 saved to: {save_path}")

    # Save CSV
    pd.DataFrame(stats_list).to_csv(
        os.path.join(data_save_dir, summary_csv_filename),
        index=False,
    )

    pd.DataFrame(detailed_list).to_csv(
        os.path.join(data_save_dir, detailed_csv_filename),
        index=False,
    )

    plt.show()

except Exception as e:
    print(f"[Error] {e}")
    import traceback
    traceback.print_exc()
