# ==========================================================
# This script plots two line-figure groups fig4 and fig5:
#
# Figure A: LCC as a function of EIS under different Ts regimes.
# Figure B: LCC as a function of Ts under different EIS regimes.
#
# Outputs:
# 1. PNG figures
# 2. CSV files containing plotted point information
#
# Colors:
# CALIPSO: #3C5488
# MODIS:   #E64B35
# ==========================================================

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import traceback
from pathlib import Path

# ==========================================================
# 1. File Paths
# ==========================================================
project_dir = Path(__file__).resolve().parents[2]
file_path = project_dir / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc"
output_dir = project_dir / "outputs"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ==========================================================
# 2. Global Settings
# ==========================================================
min_data_points = 300

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["axes.linewidth"] = 1.5
plt.rcParams["xtick.major.width"] = 1.5
plt.rcParams["ytick.major.width"] = 1.5
plt.rcParams["axes.unicode_minus"] = False

TITLE_SIZE = 24
LABEL_SIZE = 24
TICK_SIZE = 24
LEGEND_SIZE = 17

# ==========================================================
# 3. Visual Style
# ==========================================================
calipso_color = "#3C5488"
modis_color = "#E64B35"

style_calipso = {
    "color": calipso_color,
    "linestyle": "-",
    "marker": "o",
    "linewidth": 4,
    "markersize": 13,
    "markerfacecolor": calipso_color,
    "markeredgecolor": "black",
    "markeredgewidth": 0.6,
    "label": "CALIPSO",
}

style_modis = {
    "color": modis_color,
    "linestyle": "-",
    "marker": "s",
    "linewidth": 4,
    "markersize": 13,
    "markerfacecolor": modis_color,
    "markeredgecolor": "black",
    "markeredgewidth": 0.6,
    "label": "MODIS",
}

# ==========================================================
# 4. Binning Parameters
# ==========================================================

# Figure A: LCC vs EIS under different Ts regimes
fig_eis_bins = np.arange(-5, 26, 2)
fig_eis_centers = (fig_eis_bins[:-1] + fig_eis_bins[1:]) / 2

ts_regimes = [
    {"label": "Ts < 4°C", "range": (-40, 4)},
    {"label": "4°C ≤ Ts < 15°C", "range": (4, 15)},
    {"label": "15°C ≤ Ts < 24°C", "range": (15, 24)},
    {"label": "Ts ≥ 24°C", "range": (24, 35)},
]

# Figure B: LCC vs Ts under different EIS regimes
fig_ts_bins = np.arange(-30, 31, 3)
fig_ts_centers = (fig_ts_bins[:-1] + fig_ts_bins[1:]) / 2

eis_regimes = [
    {"label": "EIS < 0 K", "range": (-10, 0)},
    {"label": "0 K ≤ EIS < 10 K", "range": (0, 10)},
    {"label": "EIS ≥ 10 K", "range": (10, 30)},
]

# ==========================================================
# 5. Helper Functions
# ==========================================================
def get_mean_std_count(stats, var_name, min_count):
    """
    Calculate bin mean, standard deviation, and sample count.

    Input LCC is assumed to be fraction.
    Output LCC is converted to percentage.
    """
    mean = stats[var_name]["mean"].values.astype(float) * 100.0
    std = stats[var_name]["std"].values.astype(float) * 100.0
    count = stats[var_name]["count"].values.astype(float)

    mask = count >= min_count

    mean[~mask] = np.nan
    std[~mask] = np.nan
    count[~mask] = np.nan

    return mean, std, count


def make_point_csv_rows(
    figure_name,
    panel_label,
    regime_label,
    x_name,
    x_centers,
    sensor_name,
    mean,
    std,
    count,
):
    rows = []

    for i in range(len(x_centers)):
        rows.append(
            {
                "figure": figure_name,
                "panel": panel_label,
                "regime": regime_label,
                "x_variable": x_name,
                "x_center": x_centers[i],
                "sensor": sensor_name,
                "lcc_mean_percent": mean[i],
                "lcc_std_percent": std[i],
                "sample_count": count[i],
            }
        )

    return rows


# ==========================================================
# 6. Main Program
# ==========================================================
def main():
    print(f"Loading file: {file_path}")

    try:
        # --------------------------------------------------
        # 6.1 Load and preprocess data
        # --------------------------------------------------
        with xr.open_dataset(file_path) as ds:

            if "sic" in ds:
                ocean_mask = ds["sic"].mean(dim="time").notnull()
            else:
                ocean_mask = ds["calipso_lcc"].mean(dim="time").notnull()

            required_vars = ["calipso_lcc", "modis_lcc", "ts", "eis"]

            for v in required_vars:
                if v not in ds:
                    raise ValueError(f"Missing variable in dataset: {v}")

            ds_analysis = ds[required_vars]
            ds_ocean = ds_analysis.where(ocean_mask)

            # Convert Ts from K to °C if needed
            if float(ds_ocean["ts"].max()) > 200:
                ds_ocean["ts"] = ds_ocean["ts"] - 273.15

            print("Converting to DataFrame...")
            stacked_ds = ds_ocean.stack(point=("time", "latitude", "longitude"))
            df_main = stacked_ds.to_dataframe()
            df_main.dropna(inplace=True)

            print(f"Data ready. Valid samples: {len(df_main):,}")

        # ==================================================
        # 7. Figure A: LCC vs EIS under different Ts regimes
        # ==================================================
        print("\n" + "=" * 55)
        print("Generating Figure A: LCC vs EIS under different Ts regimes")
        print("=" * 55)

        fig_a, axes_a = plt.subplots(
            nrows=2,
            ncols=2,
            figsize=(18, 14),
            sharex=True,
            sharey=True,
        )

        axes_a_flat = axes_a.flatten()
        subplot_labels_a = ["(a)", "(b)", "(c)", "(d)"]

        rows_a = []

        for idx, regime in enumerate(ts_regimes):
            ax = axes_a_flat[idx]

            df_sub = df_main[
                (df_main["ts"] >= regime["range"][0])
                & (df_main["ts"] < regime["range"][1])
            ]

            if df_sub.empty:
                continue

            eis_cut = pd.cut(
                df_sub["eis"],
                bins=fig_eis_bins,
                right=False,
            )

            stats = df_sub.groupby(eis_cut, observed=False)[
                ["calipso_lcc", "modis_lcc"]
            ].agg(["mean", "std", "count"])

            y_calipso, std_calipso, n_calipso = get_mean_std_count(
                stats, "calipso_lcc", min_data_points
            )

            y_modis, std_modis, n_modis = get_mean_std_count(
                stats, "modis_lcc", min_data_points
            )

            rows_a.extend(
                make_point_csv_rows(
                    figure_name="Figure_A_LCC_vs_EIS_by_Ts_Regimes",
                    panel_label=subplot_labels_a[idx],
                    regime_label=regime["label"],
                    x_name="EIS",
                    x_centers=fig_eis_centers,
                    sensor_name="CALIPSO",
                    mean=y_calipso,
                    std=std_calipso,
                    count=n_calipso,
                )
            )

            rows_a.extend(
                make_point_csv_rows(
                    figure_name="Figure_A_LCC_vs_EIS_by_Ts_Regimes",
                    panel_label=subplot_labels_a[idx],
                    regime_label=regime["label"],
                    x_name="EIS",
                    x_centers=fig_eis_centers,
                    sensor_name="MODIS",
                    mean=y_modis,
                    std=std_modis,
                    count=n_modis,
                )
            )

            ax.plot(fig_eis_centers, y_calipso, **style_calipso)
            ax.plot(fig_eis_centers, y_modis, **style_modis)

            ax.set_title(
                f"{subplot_labels_a[idx]} {regime['label']}",
                fontsize=TITLE_SIZE,
                fontweight="bold",
                fontfamily="Times New Roman",
                pad=10,
            )

            ax.set_ylim(0, 100)
            ax.set_xlim(-5, 25)
            ax.grid(True, linestyle="--", alpha=0.4)

            ax.tick_params(
                axis="both",
                which="major",
                labelsize=TICK_SIZE,
                width=1.8,
                length=6,
            )

            for t in ax.get_xticklabels() + ax.get_yticklabels():
                t.set_fontweight("bold")
                t.set_fontfamily("Times New Roman")

            if idx in [2, 3]:
                ax.set_xlabel(
                    "EIS (K)",
                    fontsize=LABEL_SIZE,
                    fontweight="bold",
                    fontfamily="Times New Roman",
                )

            if idx in [0, 2]:
                ax.set_ylabel(
                    "LCC (%)",
                    fontsize=LABEL_SIZE,
                    fontweight="bold",
                    fontfamily="Times New Roman",
                )

            ax.legend(
                fontsize=LEGEND_SIZE,
                loc="best",
                frameon=False,
            )

        save_path_a_png = os.path.join(
            output_dir,
            "Figure_LCC_vs_EIS_by_Ts_Regimes_ColorOnly.png",
        )

        save_path_a_csv = os.path.join(
            output_dir,
            "Figure_LCC_vs_EIS_by_Ts_Regimes_ColorOnly_points.csv",
        )

        plt.tight_layout()
        plt.savefig(save_path_a_png, dpi=300, bbox_inches="tight")

        pd.DataFrame(rows_a).to_csv(save_path_a_csv, index=False, encoding="utf-8-sig")

        print(f"Saved figure: {save_path_a_png}")
        print(f"Saved CSV:    {save_path_a_csv}")

        plt.show()

        # ==================================================
        # 8. Figure B: LCC vs Ts under different EIS regimes
        # ==================================================
        print("\n" + "=" * 55)
        print("Generating Figure B: LCC vs Ts under different EIS regimes")
        print("=" * 55)

        fig_b, axes_b = plt.subplots(
            nrows=1,
            ncols=3,
            figsize=(22, 7),
            sharey=True,
        )

        subplot_labels_b = ["(a)", "(b)", "(c)"]

        rows_b = []

        for idx, regime in enumerate(eis_regimes):
            ax = axes_b[idx]

            df_sub = df_main[
                (df_main["eis"] >= regime["range"][0])
                & (df_main["eis"] < regime["range"][1])
            ]

            if df_sub.empty:
                continue

            ts_cut = pd.cut(
                df_sub["ts"],
                bins=fig_ts_bins,
                right=False,
            )

            stats = df_sub.groupby(ts_cut, observed=False)[
                ["calipso_lcc", "modis_lcc"]
            ].agg(["mean", "std", "count"])

            y_calipso, std_calipso, n_calipso = get_mean_std_count(
                stats, "calipso_lcc", min_data_points
            )

            y_modis, std_modis, n_modis = get_mean_std_count(
                stats, "modis_lcc", min_data_points
            )

            rows_b.extend(
                make_point_csv_rows(
                    figure_name="Figure_B_LCC_vs_Ts_by_EIS_Regimes",
                    panel_label=subplot_labels_b[idx],
                    regime_label=regime["label"],
                    x_name="Ts",
                    x_centers=fig_ts_centers,
                    sensor_name="CALIPSO",
                    mean=y_calipso,
                    std=std_calipso,
                    count=n_calipso,
                )
            )

            rows_b.extend(
                make_point_csv_rows(
                    figure_name="Figure_B_LCC_vs_Ts_by_EIS_Regimes",
                    panel_label=subplot_labels_b[idx],
                    regime_label=regime["label"],
                    x_name="Ts",
                    x_centers=fig_ts_centers,
                    sensor_name="MODIS",
                    mean=y_modis,
                    std=std_modis,
                    count=n_modis,
                )
            )

            ax.plot(fig_ts_centers, y_calipso, **style_calipso)
            ax.plot(fig_ts_centers, y_modis, **style_modis)

            ax.set_title(
                f"{subplot_labels_b[idx]} {regime['label']}",
                fontsize=TITLE_SIZE,
                fontweight="bold",
                fontfamily="Times New Roman",
                pad=10,
            )

            ax.set_xlabel(
                "Ts (°C)",
                fontsize=LABEL_SIZE,
                fontweight="bold",
                fontfamily="Times New Roman",
            )

            ax.set_xlim(-30, 30)
            ax.set_ylim(0, 100)
            ax.grid(True, linestyle="--", alpha=0.4)

            ax.tick_params(
                axis="both",
                which="major",
                labelsize=TICK_SIZE,
                width=1.8,
                length=6,
            )

            for t in ax.get_xticklabels() + ax.get_yticklabels():
                t.set_fontweight("bold")
                t.set_fontfamily("Times New Roman")

            if idx == 0:
                ax.set_ylabel(
                    "LCC (%)",
                    fontsize=LABEL_SIZE,
                    fontweight="bold",
                    fontfamily="Times New Roman",
                )

            ax.legend(
                fontsize=LEGEND_SIZE,
                loc="best",
                frameon=False,
            )

        save_path_b_png = os.path.join(
            output_dir,
            "Figure_LCC_vs_Ts_by_EIS_Regimes_ColorOnly.png",
        )

        save_path_b_csv = os.path.join(
            output_dir,
            "Figure_LCC_vs_Ts_by_EIS_Regimes_ColorOnly_points.csv",
        )

        plt.tight_layout(w_pad=4.0)
        plt.savefig(save_path_b_png, dpi=300, bbox_inches="tight")

        pd.DataFrame(rows_b).to_csv(save_path_b_csv, index=False, encoding="utf-8-sig")

        print(f"Saved figure: {save_path_b_png}")
        print(f"Saved CSV:    {save_path_b_csv}")

        plt.show()

    except Exception as e:
        print(f"Error: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
