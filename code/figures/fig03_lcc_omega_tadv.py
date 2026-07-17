# ==========================================================
# plot_fig3_lcc_omega_tadv_state_space.py
#
# Purpose
# -------
# This script generates Fig. 3 of the manuscript.
# The figure shows low-cloud cover (LCC) distributions in the
# Ts-EIS state space under different large-scale vertical-motion
# and horizontal temperature-advection regimes.
#
# Figure layout
# -------------
# The figure contains 3 rows and 4 columns:
#
#   Columns 1-2:
#       LCC distributions separated by 700 hPa vertical velocity.
#       Column 1 is based on CALIPSO LCC.
#       Column 2 is based on MODIS LCC.
#
#   Columns 3-4:
#       LCC distributions separated by horizontal temperature advection.
#       Column 3 is based on CALIPSO LCC.
#       Column 4 is based on MODIS LCC.
#
# Rows:
#   Row 1:
#       Ascending regime (omega700 < 0) and cold-advection regime
#       (T_adv < 0).
#
#   Row 2:
#       Subsiding regime (omega700 > 0) and warm-advection regime
#       (T_adv > 0).
#
#   Row 3:
#       Difference fields:
#       - omega decomposition: subsiding minus ascending
#       - T_adv decomposition: warm advection minus cold advection
#
# Variables
# ---------
# LCC:
#   calipso_lcc and modis_lcc are assumed to be fractions and are
#   converted to percentages for plotting.
#
# Ts:
#   Surface temperature. If stored in K, it is converted to °C.
#
# EIS:
#   Estimated inversion strength, in K.
#
# omega700:
#   700 hPa vertical velocity. Negative values represent ascending
#   motion, and positive values represent subsiding motion.
#
# T_adv:
#   Horizontal temperature advection. Negative values represent cold
#   advection, and positive values represent warm advection.
#
# Outputs
# -------
# 1. Figure3_LCC_Omega_Tadv_StateSpace.png
# 2. Figure3_LCC_Omega_Tadv_StateSpace_Data.csv
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
from matplotlib.ticker import MultipleLocator

# ==========================================================
# 0. Ignore warnings
# ==========================================================
warnings.simplefilter(action="ignore", category=FutureWarning)

# ==========================================================
# 1. Path settings
# ==========================================================
project_dir = Path(__file__).resolve().parents[2]
file_path_tadv = project_dir / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc"
file_path_omega = file_path_tadv
output_dir = project_dir / "outputs"

output_filename = "Figure3_LCC_Omega_Tadv_StateSpace.png"
csv_filename = "Figure3_LCC_Omega_Tadv_StateSpace_Data.csv"

omega_var_name = "omega700"
tadv_var_name = "tadv"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ==========================================================
# 2. Font and plotting settings
# ==========================================================
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Times New Roman"
plt.rcParams["mathtext.it"] = "Times New Roman:italic"
plt.rcParams["mathtext.bf"] = "Times New Roman:bold"
plt.rcParams["axes.unicode_minus"] = False

title_fontsize = 19
label_fontsize = 20
tick_fontsize = 18
font_weight = "bold"

min_data_points = 100

# ==========================================================
# 3. Load data
# ==========================================================
print(f"[Log] Loading Tadv file: {file_path_tadv}")

try:
    with xr.open_dataset(file_path_tadv) as ds_tadv:

        required_vars = [
            "calipso_lcc",
            "modis_lcc",
            "ts",
            "eis",
            tadv_var_name,
        ]

        for v in required_vars:
            if v not in ds_tadv:
                raise ValueError(f"Missing variable in Tadv file: {v}")

        # If the Tadv file already contains omega700, use it directly.
        if omega_var_name in ds_tadv:
            print("[Log] omega700 found in the Tadv file. Using this file directly.")

            ds = ds_tadv[
                [
                    "calipso_lcc",
                    "modis_lcc",
                    "ts",
                    "eis",
                    omega_var_name,
                    tadv_var_name,
                ]
            ].load()

        # The archived dataset labels 700-hPa vertical velocity as "omega".
        elif "omega" in ds_tadv:
            print("[Log] omega found in the archived dataset; using it as omega700.")
            ds = ds_tadv[
                [
                    "calipso_lcc",
                    "modis_lcc",
                    "ts",
                    "eis",
                    "omega",
                    tadv_var_name,
                ]
            ].rename({"omega": omega_var_name}).load()

        # If no omega variable is included, supplement it from another file.
        else:
            print("[Log] omega700 not found in the Tadv file. Loading omega700 from merged file.")

            with xr.open_dataset(file_path_omega) as ds_omega:
                if omega_var_name not in ds_omega:
                    raise ValueError(f"Missing variable in omega file: {omega_var_name}")

                ds_main = ds_tadv[
                    [
                        "calipso_lcc",
                        "modis_lcc",
                        "ts",
                        "eis",
                        tadv_var_name,
                    ]
                ]

                ds_omg = ds_omega[[omega_var_name]]

                # Align by common coordinates.
                ds_main, ds_omg = xr.align(ds_main, ds_omg, join="inner")
                ds = xr.merge([ds_main, ds_omg]).load()

        # ==================================================
        # 4. Ocean / valid-area mask
        # ==================================================
        if "sic" in ds_tadv:
            ocean_mask = ds_tadv["sic"].mean(dim="time").notnull()
            ocean_mask, ds = xr.align(ocean_mask, ds, join="inner")
            ds_ocean = ds.where(ocean_mask)
        else:
            ds_ocean = ds.where(ds["calipso_lcc"].mean(dim="time").notnull())

        # Convert Ts from K to °C if needed.
        if float(ds_ocean["ts"].max()) > 200:
            ds_ocean["ts"] = ds_ocean["ts"] - 273.15

        print("[Log] Converting to DataFrame...")

        df_main = ds_ocean.stack(
            point=("time", "latitude", "longitude")
        ).to_dataframe()

        df_main.dropna(
            subset=[
                "ts",
                "eis",
                omega_var_name,
                tadv_var_name,
                "calipso_lcc",
                "modis_lcc",
            ],
            inplace=True,
        )

        print(f"[Log] Valid samples: {len(df_main)}")

        # ==================================================
        # 5. Ts-EIS binning
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

        def get_matrix(df, col):
            """
            Calculate mean LCC in each Ts-EIS bin.

            LCC is assumed to be stored as a fraction and is converted
            to percentage. Bins with fewer than min_data_points samples
            are masked as NaN.
            """
            mean = df.groupby(["eis_bin", "ts_bin"], observed=False)[col].mean() * 100
            count = df.groupby(["eis_bin", "ts_bin"], observed=False)[col].count()

            mat_mean = mean.unstack()
            mat_count = count.unstack()

            return mat_mean.where(mat_count >= min_data_points)

        # ==================================================
        # 6. Split samples by omega700 and T_adv regimes
        # ==================================================
        print("[Log] Splitting samples by omega700 and T_adv regimes...")

        df_asc = df_main[df_main[omega_var_name] < 0]
        df_sub = df_main[df_main[omega_var_name] > 0]

        df_cold = df_main[df_main[tadv_var_name] < 0]
        df_warm = df_main[df_main[tadv_var_name] > 0]

        print("[Log] Calculating CALIPSO / MODIS Ts-EIS matrices...")

        # -------- omega decomposition --------
        mat_cal_asc = get_matrix(df_asc, "calipso_lcc")
        mat_mod_asc = get_matrix(df_asc, "modis_lcc")

        mat_cal_sub = get_matrix(df_sub, "calipso_lcc")
        mat_mod_sub = get_matrix(df_sub, "modis_lcc")

        mat_cal_omega_diff = mat_cal_sub - mat_cal_asc
        mat_mod_omega_diff = mat_mod_sub - mat_mod_asc

        # -------- T_adv decomposition --------
        mat_cal_cold = get_matrix(df_cold, "calipso_lcc")
        mat_mod_cold = get_matrix(df_cold, "modis_lcc")

        mat_cal_warm = get_matrix(df_warm, "calipso_lcc")
        mat_mod_warm = get_matrix(df_warm, "modis_lcc")

        mat_cal_tadv_diff = mat_cal_warm - mat_cal_cold
        mat_mod_tadv_diff = mat_mod_warm - mat_mod_cold

        # ==================================================
        # 7. Export plotting data to CSV
        # ==================================================
        print("[Log] Exporting plotting data to CSV...")

        df_csv = pd.DataFrame({
            "CALIPSO_Ascending_omega_lt_0": mat_cal_asc.stack(dropna=False),
            "MODIS_Ascending_omega_lt_0": mat_mod_asc.stack(dropna=False),
            "CALIPSO_Subsiding_omega_gt_0": mat_cal_sub.stack(dropna=False),
            "MODIS_Subsiding_omega_gt_0": mat_mod_sub.stack(dropna=False),
            "CALIPSO_Omega_Diff_Sub_minus_Asc": mat_cal_omega_diff.stack(dropna=False),
            "MODIS_Omega_Diff_Sub_minus_Asc": mat_mod_omega_diff.stack(dropna=False),

            "CALIPSO_ColdAdv_Tadv_lt_0": mat_cal_cold.stack(dropna=False),
            "MODIS_ColdAdv_Tadv_lt_0": mat_mod_cold.stack(dropna=False),
            "CALIPSO_WarmAdv_Tadv_gt_0": mat_cal_warm.stack(dropna=False),
            "MODIS_WarmAdv_Tadv_gt_0": mat_mod_warm.stack(dropna=False),
            "CALIPSO_Tadv_Diff_Warm_minus_Cold": mat_cal_tadv_diff.stack(dropna=False),
            "MODIS_Tadv_Diff_Warm_minus_Cold": mat_mod_tadv_diff.stack(dropna=False),
        }).reset_index()

        df_csv.rename(
            columns={
                "eis_bin": "EIS_Center",
                "ts_bin": "Ts_Center",
            },
            inplace=True,
        )

        csv_save_path = os.path.join(output_dir, csv_filename)
        df_csv.to_csv(csv_save_path, index=False, float_format="%.2f")

        print(f"[Success] Plotting data saved to: {csv_save_path}")

        # ==================================================
        # 8. Plotting: 3 rows × 4 columns
        # ==================================================
        print("[Log] Plotting combined Fig. 3...")

        fig, axes = plt.subplots(
            nrows=3,
            ncols=4,
            figsize=(28, 18),
            constrained_layout=False,
        )

        plots = [
            # Row 1: negative regimes
            {
                "ax": axes[0, 0],
                "data": mat_cal_asc,
                "label": "(a)",
                "title": "CALIPSO - Ascending (ω < 0)",
                "type": "abs",
            },
            {
                "ax": axes[0, 1],
                "data": mat_mod_asc,
                "label": "(b)",
                "title": "MODIS - Ascending (ω < 0)",
                "type": "abs",
            },
            {
                "ax": axes[0, 2],
                "data": mat_cal_cold,
                "label": "(c)",
                "title": r"CALIPSO - Cold Adv. ($\mathbf{T}_{\mathbf{adv}}$ < 0)",
                "type": "abs",
            },
            {
                "ax": axes[0, 3],
                "data": mat_mod_cold,
                "label": "(d)",
                "title": r"MODIS - Cold Adv. ($\mathbf{T}_{\mathbf{adv}}$ < 0)",
                "type": "abs",
            },

            # Row 2: positive regimes
            {
                "ax": axes[1, 0],
                "data": mat_cal_sub,
                "label": "(e)",
                "title": "CALIPSO - Subsiding (ω > 0)",
                "type": "abs",
            },
            {
                "ax": axes[1, 1],
                "data": mat_mod_sub,
                "label": "(f)",
                "title": "MODIS - Subsiding (ω > 0)",
                "type": "abs",
            },
            {
                "ax": axes[1, 2],
                "data": mat_cal_warm,
                "label": "(g)",
                "title": r"CALIPSO - Warm Adv. ($\mathbf{T}_{\mathbf{adv}}$ > 0)",
                "type": "abs",
            },
            {
                "ax": axes[1, 3],
                "data": mat_mod_warm,
                "label": "(h)",
                "title": r"MODIS - Warm Adv. ($\mathbf{T}_{\mathbf{adv}}$ > 0)",
                "type": "abs",
            },

            # Row 3: differences
            {
                "ax": axes[2, 0],
                "data": mat_cal_omega_diff,
                "label": "(i)",
                "title": "CALIPSO - Diff. (Sub - Asc)",
                "type": "diff",
            },
            {
                "ax": axes[2, 1],
                "data": mat_mod_omega_diff,
                "label": "(j)",
                "title": "MODIS - Diff. (Sub - Asc)",
                "type": "diff",
            },
            {
                "ax": axes[2, 2],
                "data": mat_cal_tadv_diff,
                "label": "(k)",
                "title": "CALIPSO - Diff. (Warm - Cold)",
                "type": "diff",
            },
            {
                "ax": axes[2, 3],
                "data": mat_mod_tadv_diff,
                "label": "(l)",
                "title": "MODIS - Diff. (Warm - Cold)",
                "type": "diff",
            },
        ]

        mesh_abs = None
        mesh_diff = None

        for p in plots:
            ax = p["ax"]

            if p["type"] == "abs":
                vmin, vmax = 10, 90
            else:
                vmin, vmax = -20, 20

            mesh = ax.pcolormesh(
                ts_bins,
                eis_bins,
                p["data"].values,
                cmap="RdBu_r",
                vmin=vmin,
                vmax=vmax,
                shading="auto",
            )

            if p["type"] == "abs":
                mesh_abs = mesh
            else:
                mesh_diff = mesh

            # Panel label is placed before the title.
            ax.set_title(
                f'{p["label"]} {p["title"]}',
                fontsize=title_fontsize,
                fontweight=font_weight,
                fontfamily="Times New Roman",
                pad=9,
            )

            ax.set_xlim(ts_bins[0], ts_bins[-1])
            ax.set_ylim(eis_bins[0], eis_bins[-1])

            ax.xaxis.set_major_locator(MultipleLocator(10))
            ax.yaxis.set_major_locator(MultipleLocator(5))

            ax.grid(True, linestyle="--", alpha=0.35, color="black")

            ax.tick_params(labelsize=tick_fontsize, width=1.8, length=5)

            for t in ax.get_xticklabels() + ax.get_yticklabels():
                t.set_fontweight(font_weight)
                t.set_fontfamily("Times New Roman")

            # Show y-axis labels only in the leftmost column.
            if ax in axes[:, 0]:
                ax.set_ylabel(
                    "EIS (K)",
                    fontsize=label_fontsize,
                    fontweight=font_weight,
                    fontfamily="Times New Roman",
                )
            else:
                ax.set_ylabel("")

            # Show x-axis labels only in the bottom row.
            if ax in axes[2, :]:
                ax.set_xlabel(
                    "Ts (°C)",
                    fontsize=label_fontsize,
                    fontweight=font_weight,
                    fontfamily="Times New Roman",
                )
            else:
                ax.set_xlabel("")

        # ==================================================
        # 9. Layout and colorbars
        # ==================================================
        plt.subplots_adjust(
            left=0.06,
            right=0.88,
            bottom=0.08,
            top=0.93,
            hspace=0.32,
            wspace=0.20,
        )

        # Shared LCC colorbar for absolute LCC panels in rows 1-2.
        cbar_ax_abs = fig.add_axes([0.90, 0.43, 0.018, 0.46])
        cb_abs = fig.colorbar(
            mesh_abs,
            cax=cbar_ax_abs,
            extend="both",
            ticks=np.arange(10, 91, 10),
        )

        cb_abs.set_label(
            "LCC (%)",
            fontsize=label_fontsize,
            fontweight=font_weight,
            fontfamily="Times New Roman",
        )

        cb_abs.ax.tick_params(labelsize=tick_fontsize)

        for t in cb_abs.ax.get_yticklabels():
            t.set_fontweight(font_weight)
            t.set_fontfamily("Times New Roman")

        # Shared difference colorbar for difference panels in row 3.
        cbar_ax_diff = fig.add_axes([0.90, 0.08, 0.018, 0.23])
        cb_diff = fig.colorbar(
            mesh_diff,
            cax=cbar_ax_diff,
            extend="both",
            ticks=np.arange(-20, 21, 10),
        )

        cb_diff.set_label(
            "Diff. in LCC (%)",
            fontsize=label_fontsize,
            fontweight=font_weight,
            fontfamily="Times New Roman",
        )

        cb_diff.ax.tick_params(labelsize=tick_fontsize)

        for t in cb_diff.ax.get_yticklabels():
            t.set_fontweight(font_weight)
            t.set_fontfamily("Times New Roman")

        # ==================================================
        # 10. Save
        # ==================================================
        save_path = os.path.join(output_dir, output_filename)

        plt.savefig(save_path, dpi=300, bbox_inches="tight")

        print(f"\n[Success] Combined Fig. 3 saved to: {save_path}")

        plt.show()

except Exception as e:
    print(f"\n[Error] {e}")
    import traceback
    traceback.print_exc()
