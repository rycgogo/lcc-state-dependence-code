# ==========================================================
# plot_fig2_lcc_sic_frequency_ts_eis.py
#
# Purpose
# -------
# This script generates Fig. 2 of the manuscript.
# The figure shows the distributions of low-cloud cover (LCC),
# sea ice concentration (SIC), and sample frequency in the
# surface temperature–estimated inversion strength (Ts–EIS)
# state space.
#
# Figure layout
# -------------
# The figure contains 2 rows and 2 columns:
#
#   (a) CALIPSO LCC
#   (b) MODIS LCC
#   (c) Sea ice concentration (SIC)
#   (d) Frequency of occurrence
#
# Processing notes
# ----------------
# - CALIPSO and MODIS LCC are assumed to be stored as fractions
#   and are converted to percentages.
# - Ts is converted from K to °C if necessary.
# - SIC is converted from fraction to percent if necessary.
# - Samples are binned in the Ts–EIS state space.
# - Mean fields are masked where the sample count is lower than
#   min_data_points.
# - The frequency panel uses LogNorm so that rare and frequent
#   states can be shown together.
#
# Outputs
# -------
# 1. Figure2_LCC_SIC_Frequency_TsEIS.png
# 2. Figure2_LCC_SIC_Frequency_TsEIS_Data.csv
#
# Recommended code path
# ---------------------
# E:\01_Master_Research\05_Code\03_Plotting\
# ==========================================================

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import warnings
from pathlib import Path
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator

# ==========================================================
# 0. Ignore warnings
# ==========================================================
warnings.simplefilter(action="ignore", category=FutureWarning)

# ==========================================================
# 1. Path settings
# ==========================================================
project_dir = Path(__file__).resolve().parents[2]
file_path = project_dir / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc"
output_dir = project_dir / "outputs"

output_filename = "Figure2_LCC_SIC_Frequency_TsEIS.png"
csv_filename = "Figure2_LCC_SIC_Frequency_TsEIS_Data.csv"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ==========================================================
# 2. Plot style settings
# ==========================================================
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

title_fontsize = 24
label_fontsize = 24
tick_fontsize = 22
cbar_tick_fontsize = 20
font_weight = "bold"

min_data_points = 100

# All panels use jet for consistency with Figure 1
cmap_all = "jet"

# Unified colorbar style
cbar_fraction = 0.050
cbar_pad = 0.040

# ==========================================================
# 3. Data processing
# ==========================================================
print(f"Loading file: {file_path}")

try:
    with xr.open_dataset(file_path) as ds:

        # Automatically identify SIC variable name
        if "siconc" in ds:
            sic_var = "siconc"
        elif "sic" in ds:
            sic_var = "sic"
        else:
            raise ValueError("No SIC variable found in dataset. Expected 'sic' or 'siconc'.")

        required_vars = [
            "calipso_lcc",
            "modis_lcc",
            "ts",
            "eis",
            sic_var,
        ]

        for v in required_vars:
            if v not in ds:
                raise ValueError(f"Missing variable in dataset: {v}")

        ds_subset = ds[required_vars]

        # Ocean / valid-data mask
        ocean_mask = ds[sic_var].mean(dim="time").notnull()
        ds_ocean = ds_subset.where(ocean_mask)

        # Convert Ts from K to °C if necessary
        if float(ds_ocean["ts"].max()) > 200:
            ds_ocean["ts"] = ds_ocean["ts"] - 273.15

        # Convert SIC to percent if stored as 0-1 fraction
        if float(ds_ocean[sic_var].max()) <= 1.1:
            ds_ocean[sic_var] = ds_ocean[sic_var] * 100.0

        print("Converting data to DataFrame...")

        df_main = ds_ocean.stack(
            point=("time", "latitude", "longitude")
        ).to_dataframe()

        df_main.dropna(
            subset=[
                "ts",
                "eis",
                "calipso_lcc",
                "modis_lcc",
                sic_var,
            ],
            inplace=True,
        )

        print(f"Valid sample size: {len(df_main)}")

        # ==================================================
        # 4. Binning in Ts-EIS state space
        # ==================================================
        ts_bins = np.arange(-31, 32, 2)
        eis_bins = np.arange(-5, 26, 2)

        ts_centers = (ts_bins[:-1] + ts_bins[1:]) / 2
        eis_centers = (eis_bins[:-1] + eis_bins[1:]) / 2

        df_main["ts_bin"] = pd.cut(
            df_main["ts"],
            bins=ts_bins,
            labels=ts_centers,
            right=False,
        )

        df_main["eis_bin"] = pd.cut(
            df_main["eis"],
            bins=eis_bins,
            labels=eis_centers,
            right=False,
        )

        def get_stat_matrix(df, col, stat="mean", multiply_100=False):
            """
            Calculate mean or count in each Ts-EIS bin.
            """
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

                if stat == "mean":
                    val = df.groupby(["eis_bin", "ts_bin"], observed=False)[col].mean()
                elif stat == "count":
                    val = df.groupby(["eis_bin", "ts_bin"], observed=False)[col].count()
                else:
                    raise ValueError("stat must be either 'mean' or 'count'.")

            if multiply_100:
                val = val * 100.0

            return val.unstack()

        print("Calculating statistics...")

        # LCC: convert from fraction to percent
        mat_cal = get_stat_matrix(
            df_main,
            "calipso_lcc",
            stat="mean",
            multiply_100=True,
        )

        mat_mod = get_stat_matrix(
            df_main,
            "modis_lcc",
            stat="mean",
            multiply_100=True,
        )

        # SIC is already in percent after the unit check above
        mat_sic = get_stat_matrix(
            df_main,
            sic_var,
            stat="mean",
            multiply_100=False,
        )

        # Frequency in percent
        mat_count = get_stat_matrix(
            df_main,
            "ts",
            stat="count",
        )

        total_count = mat_count.sum().sum()
        mat_freq = (mat_count / total_count) * 100.0

        # Mask empty cells
        mat_freq = mat_freq.where(mat_freq > 0)
        mat_sic = mat_sic.where(mat_sic >= 0)

        # Apply minimum-sample threshold to mean fields
        mask = mat_count < min_data_points
        mat_cal = mat_cal.where(~mask)
        mat_mod = mat_mod.where(~mask)
        mat_sic = mat_sic.where(~mask)

        # ==================================================
        # 5. Export plotting data
        # ==================================================
        print("Exporting plotting data...")

        df_csv = pd.DataFrame({
            "CALIPSO_LCC_percent": mat_cal.stack(dropna=False),
            "MODIS_LCC_percent": mat_mod.stack(dropna=False),
            "SIC_percent": mat_sic.stack(dropna=False),
            "Frequency_percent": mat_freq.stack(dropna=False),
        }).reset_index()

        df_csv.rename(
            columns={
                "eis_bin": "EIS_Center",
                "ts_bin": "Ts_Center",
            },
            inplace=True,
        )

        csv_save_path = os.path.join(output_dir, csv_filename)
        df_csv.to_csv(csv_save_path, index=False, float_format="%.4f")

        print(f"[Success] Plotting data saved to: {csv_save_path}")

        # ==================================================
        # 6. Plotting
        # ==================================================
        print("Plotting Figure 2...")

        fig, axes = plt.subplots(
            nrows=2,
            ncols=2,
            figsize=(18, 16),
            constrained_layout=False,
        )

        plot_configs = [
            {
                "ax": axes[0, 0],
                "data": mat_cal,
                "label": "(a)",
                "title": "CALIPSO LCC",
                "cmap": cmap_all,
                "vmin": 0,
                "vmax": 100,
                "ticks": np.arange(0, 101, 20),
                "cbar_label": "LCC (%)",
                "log": False,
            },
            {
                "ax": axes[0, 1],
                "data": mat_mod,
                "label": "(b)",
                "title": "MODIS LCC",
                "cmap": cmap_all,
                "vmin": 0,
                "vmax": 100,
                "ticks": np.arange(0, 101, 20),
                "cbar_label": "LCC (%)",
                "log": False,
            },
            {
                "ax": axes[1, 0],
                "data": mat_sic,
                "label": "(c)",
                "title": "Sea Ice Concentration (SIC)",
                "cmap": cmap_all,
                "vmin": 0,
                "vmax": 100,
                "ticks": np.arange(0, 101, 20),
                "cbar_label": "SIC (%)",
                "log": False,
            },
            {
                "ax": axes[1, 1],
                "data": mat_freq,
                "label": "(d)",
                "title": "Frequency of Occurrence",
                "cmap": cmap_all,
                "vmin": None,
                "vmax": None,
                "ticks": None,
                "cbar_label": "Frequency (%)",
                "log": True,
            },
        ]

        for cfg in plot_configs:
            ax = cfg["ax"]
            data = cfg["data"]

            if cfg["log"]:
                valid_values = data.values[np.isfinite(data.values)]
                valid_values = valid_values[valid_values > 0]

                if len(valid_values) == 0:
                    raise ValueError("Frequency matrix has no positive values.")

                vmin_log = max(np.nanmin(valid_values), 1e-4)
                vmax_log = np.nanmax(valid_values)

                mesh = ax.pcolormesh(
                    ts_bins,
                    eis_bins,
                    data.values,
                    cmap=cfg["cmap"],
                    norm=LogNorm(vmin=vmin_log, vmax=vmax_log),
                    shading="auto",
                )

            else:
                mesh = ax.pcolormesh(
                    ts_bins,
                    eis_bins,
                    data.values,
                    cmap=cfg["cmap"],
                    vmin=cfg["vmin"],
                    vmax=cfg["vmax"],
                    shading="auto",
                )

            # Panel labels are placed before the subplot title
            ax.set_title(
                f'{cfg["label"]} {cfg["title"]}',
                fontsize=title_fontsize,
                fontweight=font_weight,
                fontfamily="Times New Roman",
                pad=12,
            )

            ax.set_xlabel(
                "Ts (°C)",
                fontsize=label_fontsize,
                fontweight=font_weight,
                fontfamily="Times New Roman",
            )

            ax.set_ylabel(
                "EIS (K)",
                fontsize=label_fontsize,
                fontweight=font_weight,
                fontfamily="Times New Roman",
            )

            ax.set_xlim(ts_bins[0], ts_bins[-1])
            ax.set_ylim(eis_bins[0], eis_bins[-1])

            ax.xaxis.set_major_locator(MultipleLocator(10))
            ax.yaxis.set_major_locator(MultipleLocator(5))

            ax.grid(
                True,
                linestyle="--",
                alpha=0.4,
                color="gray",
            )

            ax.tick_params(
                labelsize=tick_fontsize,
                width=2,
                length=6,
            )

            for t in ax.get_xticklabels() + ax.get_yticklabels():
                t.set_fontweight(font_weight)
                t.set_fontfamily("Times New Roman")

            cbar = fig.colorbar(
                mesh,
                ax=ax,
                orientation="vertical",
                fraction=cbar_fraction,
                pad=cbar_pad,
                ticks=cfg["ticks"],
                extend="max" if cfg["log"] else "neither",
            )

            cbar.set_label(
                cfg["cbar_label"],
                fontsize=label_fontsize,
                fontweight=font_weight,
                fontfamily="Times New Roman",
            )

            cbar.ax.tick_params(
                labelsize=cbar_tick_fontsize,
                width=2,
            )

            for t in cbar.ax.get_yticklabels():
                t.set_fontweight(font_weight)
                t.set_fontfamily("Times New Roman")

        # Subplot spacing
        plt.subplots_adjust(
            left=0.08,
            right=0.95,
            bottom=0.08,
            top=0.93,
            wspace=0.30,
            hspace=0.30,
        )

        # ==================================================
        # 7. Save figure
        # ==================================================
        save_path = os.path.join(output_dir, output_filename)

        plt.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight",
        )

        print(f"[Success] Figure 2 saved to: {save_path}")

        plt.show()

except Exception as e:
    print(f"[Error] {e}")
    import traceback
    traceback.print_exc()
