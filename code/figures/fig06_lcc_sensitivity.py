# ==========================================================
# plot_fig6_lcc_sensitivity_global_state_space.py
#
# Purpose
# -------
# This script generates Fig. 6 of the manuscript.
# The figure combines two diagnostics of low-cloud cover (LCC)
# sensitivity to five cloud-controlling factors (CCFs).
#
# Figure layout
# -------------
# The figure contains 5 rows and 4 columns:
#
#   Columns 1-2:
#       Global maps of standardized regression coefficients.
#       Column 1 is based on CALIPSO LCC.
#       Column 2 is based on MODIS LCC.
#
#   Columns 3-4:
#       Ts-EIS state-space maps of physical sensitivities.
#       Column 3 is based on CALIPSO LCC.
#       Column 4 is based on MODIS LCC.
#
# Variables shown from top to bottom:
#   1. surface temperature (Ts)
#   2. estimated inversion strength (EIS)
#   3. 700 hPa vertical velocity (omega)
#   4. horizontal temperature advection (T_adv)
#   5. sea ice concentration (SIC)
#
# Regression design
# -----------------
# For the global maps, monthly anomalies of LCC and CCFs are
# deseasonalized and standardized at each grid box before multiple
# linear regression. The resulting coefficients are standardized
# sensitivities.
#
# For the Ts-EIS state-space panels, monthly anomalies are grouped
# into Ts-EIS bins. Within each bin, non-standardized multiple linear
# regression is performed to estimate physical sensitivities:
#
#   dLCC/dTs, dLCC/dEIS, dLCC/domega, dLCC/dT_adv, and dLCC/dSIC.
#
# Unit conversions
# ----------------
# - Ts is converted from K to °C if necessary.
# - SIC is converted from fraction to percent if necessary.
# - omega is converted from Pa/s to hPa/day for state-space regression.
# - T_adv is converted from K/s to K/day for state-space regression.
#
# Outputs
# -------
# 1. Figure6_LCC_Sensitivity_Global_StateSpace.png
# 2. Figure6_LCC_Sensitivity_StateSpace_Data.csv
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
import statsmodels.api as sm
import warnings

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

from matplotlib.ticker import MultipleLocator
from pathlib import Path

# ==========================================================
# 0. Ignore warnings
# ==========================================================
warnings.simplefilter(action="ignore", category=FutureWarning)

# ==========================================================
# 1. Path settings
# ==========================================================
project_dir = Path(__file__).resolve().parents[2]
file_path = project_dir / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc"
file_path_omega = file_path
output_dir = project_dir / "outputs"

output_filename = "Figure6_LCC_Sensitivity_Global_StateSpace.png"
csv_filename = "Figure6_LCC_Sensitivity_StateSpace_Data.csv"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ==========================================================
# ==========================================================
# 2. Font and plotting parameters
# ==========================================================
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Times New Roman"
plt.rcParams["mathtext.it"] = "Times New Roman:italic"
plt.rcParams["mathtext.bf"] = "Times New Roman:bold"
plt.rcParams["axes.unicode_minus"] = False

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

title_fontsize = 25
label_fontsize = 26
tick_fontsize = 23
cbar_label_fontsize = 25
cbar_tick_fontsize = 25

font_weight = "bold"

cmap = "RdBu_r"

# Left two columns: global standardized regression coefficient color scale
global_vmin = -1.0
global_vmax = 1.0

# Right two columns: state-space physical sensitivity color scale
# omega has been converted to hPa/day, so the range is set to ±0.20
state_vlims = {
    "ts": 5.0,       # % / K
    "eis": 5.0,      # % / K
    "omega": 0.20,   # % / (hPa/day)
    "tadv": 5.0,     # % / (K/day)
    "sic": 1.0,      # % / %
}

state_units = {
    "ts": "% / K",
    "eis": "% / K",
    "omega": "% / (hPa/day)",
    "tadv": "% / (K/day)",
    "sic": "% / %",
}

min_data_points = 300

# ==========================================================
# 3. Helper functions
# ==========================================================
def compute_edges(centers):
    """
    Calculate pcolormesh edges from grid-center coordinates.
    """
    centers = np.asarray(centers)

    if len(centers) < 2:
        return np.array([centers[0] - 0.5, centers[0] + 0.5])

    diff = centers[1] - centers[0]
    start = centers[0] - diff / 2
    edges = np.arange(len(centers) + 1) * diff + start

    return edges


def deseasonalize_robust(da):
    """
    Remove the monthly climatological seasonal cycle.

    For each grid box and calendar month, the multi-year monthly mean
    is subtracted from the original monthly time series.
    """
    climatology = da.groupby("time.month").mean(dim="time")
    aligned = climatology.sel(month=da["time.month"])
    aligned["time"] = da["time"]

    return da - aligned


def standardize_robust_xarray(da):
    """
    Standardize an xarray DataArray along the time dimension.

    This is used for the global grid-box regression, so that the
    regression coefficients represent standardized sensitivities.
    """
    std = da.std(dim="time")
    mean = da.mean(dim="time")

    return ((da - mean) / std).where(std > 1e-6)


def ols_regression_global(y, ts, eis, omega, sic, tadv):
    """
    Perform standardized multiple linear regression at one grid box.

    Inputs are standardized monthly anomalies.

    Returned coefficient order:
        ts, eis, omega, sic, tadv
    """
    mask_core = (
        ~np.isnan(y)
        & ~np.isnan(ts)
        & ~np.isnan(eis)
        & ~np.isnan(omega)
        & ~np.isnan(tadv)
    )

    if np.sum(mask_core) < 20:
        return np.array([np.nan] * 5)

    y_c = y[mask_core]
    ts_c = ts[mask_core]
    eis_c = eis[mask_core]
    omega_c = omega[mask_core]
    sic_c = sic[mask_core]
    tadv_c = tadv[mask_core]

    sic_valid_mask = ~np.isnan(sic_c)
    use_sic = False

    if np.sum(sic_valid_mask) == len(y_c):
        if np.std(sic_c) > 1e-4:
            use_sic = True

    try:
        if use_sic:
            X = sm.add_constant(
                np.column_stack((ts_c, eis_c, omega_c, sic_c, tadv_c)),
                prepend=True,
            )

            params = sm.OLS(y_c, X).fit().params[1:]

            return params

        else:
            X = sm.add_constant(
                np.column_stack((ts_c, eis_c, omega_c, tadv_c)),
                prepend=True,
            )

            params = sm.OLS(y_c, X).fit().params[1:]

            return np.array([
                params[0],     # ts
                params[1],     # eis
                params[2],     # omega
                np.nan,        # sic
                params[3],     # tadv
            ])

    except Exception:
        return np.array([np.nan] * 5)


def run_global_regression(y, ts, eis, omega, sic, tadv):
    """
    Apply standardized multiple linear regression to all grid boxes.
    """
    return xr.apply_ufunc(
        ols_regression_global,
        y,
        ts,
        eis,
        omega,
        sic,
        tadv,
        input_core_dims=[
            ["time"],
            ["time"],
            ["time"],
            ["time"],
            ["time"],
            ["time"],
        ],
        output_core_dims=[["coefficient"]],
        exclude_dims=set(("time",)),
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float],
        output_sizes={"coefficient": 5},
    )


def get_state_space_physical_sensitivities(group):
    """
    Perform non-standardized physical-unit regression within each Ts-EIS bin.

    Returned sensitivities are for CALIPSO and MODIS LCC:
        ts, eis, omega, tadv, sic

    Units:
        dLCC/dTs      : % / K
        dLCC/dEIS     : % / K
        dLCC/domega   : % / (hPa/day)
        dLCC/dT_adv   : % / (K/day)
        dLCC/dSIC     : % / %
    """
    empty_result = {
        k: np.nan
        for k in [
            "cal_beta_ts",
            "cal_beta_eis",
            "cal_beta_omega",
            "cal_beta_tadv",
            "cal_beta_sic",
            "mod_beta_ts",
            "mod_beta_eis",
            "mod_beta_omega",
            "mod_beta_tadv",
            "mod_beta_sic",
        ]
    }

    if len(group) < min_data_points:
        return pd.Series(empty_result)

    group = group.copy()

    # Determine whether SIC should be included in the regression.
    # This avoids using SIC in bins where sea ice is almost always absent
    # or almost always saturated.
    sic_water_threshold = 1.0
    sic_ice_threshold = 99.0
    saturation_limit = 0.5

    n_saturated = (
        (group["sic_total"] <= sic_water_threshold)
        | (group["sic_total"] >= sic_ice_threshold)
    ).sum()

    if (n_saturated / len(group) > saturation_limit) or (group["sic_anom"].std() < 0.1):
        calc_sic = False
    else:
        calc_sic = True

    # Unit conversion
    # omega: Pa/s -> hPa/day
    # 1 Pa/s = 86400 Pa/day = 864 hPa/day
    group["omega_anom_hPa_day"] = group["omega_anom"] * 864.0

    # Tadv: K/s -> K/day
    group["tadv_anom_daily"] = group["tadv_anom"] * 86400.0

    x_cols = [
        "ts_anom",
        "eis_anom",
        "omega_anom_hPa_day",
        "tadv_anom_daily",
    ]

    if calc_sic:
        x_cols.append("sic_anom")

    results = {}

    def fit_model(y_col, prefix):
        try:
            reg_df = group[[y_col] + x_cols].dropna()

            if len(reg_df) < min_data_points:
                for var in ["ts", "eis", "omega", "tadv", "sic"]:
                    results[f"{prefix}_beta_{var}"] = np.nan
                return

            X = sm.add_constant(reg_df[x_cols], prepend=True)
            y = reg_df[y_col]

            model = sm.OLS(y, X).fit()

            results[f"{prefix}_beta_ts"] = model.params.get("ts_anom", np.nan)
            results[f"{prefix}_beta_eis"] = model.params.get("eis_anom", np.nan)
            results[f"{prefix}_beta_omega"] = model.params.get("omega_anom_hPa_day", np.nan)
            results[f"{prefix}_beta_tadv"] = model.params.get("tadv_anom_daily", np.nan)

            if calc_sic:
                results[f"{prefix}_beta_sic"] = model.params.get("sic_anom", np.nan)
            else:
                results[f"{prefix}_beta_sic"] = np.nan

        except Exception:
            for var in ["ts", "eis", "omega", "tadv", "sic"]:
                results[f"{prefix}_beta_{var}"] = np.nan

    fit_model("cal_lcc_anom", "cal")
    fit_model("mod_lcc_anom", "mod")

    return pd.Series(results)


def plot_global_map(ax, data, title, panel_label):
    """
    Plot a global map of standardized regression coefficients.
    Panel labels are placed before titles.
    """
    ax.set_global()

    data_values = data.values

    if "longitude" in data.coords:
        lons = data["longitude"].values
    else:
        lons = data["lon"].values

    if "latitude" in data.coords:
        lats = data["latitude"].values
    else:
        lats = data["lat"].values

    data_cyclic, lons_cyclic = add_cyclic_point(data_values, coord=lons)

    mappable = ax.pcolormesh(
        compute_edges(lons_cyclic),
        compute_edges(lats),
        data_cyclic,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        vmin=global_vmin,
        vmax=global_vmax,
        shading="flat",
        zorder=1,
    )

    ax.add_feature(cfeature.LAND, facecolor="silver", edgecolor="none", zorder=2)
    ax.coastlines(linewidth=0.75, color="black", zorder=3)

    ax.gridlines(
        draw_labels=False,
        linewidth=0.45,
        color="gray",
        alpha=0.55,
        linestyle="--",
        zorder=3,
    )

    ax.set_title(
        f"{panel_label} {title}",
        fontsize=title_fontsize,
        fontweight=font_weight,
        fontfamily="Times New Roman",
        pad=5,
    )

    ax.set_aspect(1.5)

    return mappable


def plot_state_space(ax, data, title, panel_label, ts_bins, eis_bins, vlim):
    """
    Plot a Ts-EIS state-space map of physical sensitivities.
    Panel labels are placed before titles.
    """
    mesh = ax.pcolormesh(
        ts_bins,
        eis_bins,
        data.values,
        cmap=cmap,
        vmin=-vlim,
        vmax=vlim,
        shading="auto",
    )

    ax.set_title(
        f"{panel_label} {title}",
        fontsize=title_fontsize,
        fontweight=font_weight,
        fontfamily="Times New Roman",
        pad=5,
    )

    ax.set_xlim(-30, 30)
    ax.set_ylim(-5, 25)

    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(5))

    ax.grid(True, linestyle="--", alpha=0.35, color="black")

    ax.tick_params(labelsize=tick_fontsize, width=1.8, length=6)

    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_fontweight(font_weight)
        t.set_fontfamily("Times New Roman")

    ax.set_box_aspect(0.82)

    return mesh


# ==========================================================
# 4. Data loading and preprocessing
# ==========================================================
print(f"[Log] Loading file: {file_path}")

try:
    with xr.open_dataset(file_path) as ds_input:

        required_vars = [
            "calipso_lcc",
            "modis_lcc",
            "ts",
            "eis",
            "sic",
            "tadv",
        ]

        for v in required_vars:
            if v not in ds_input:
                raise ValueError(f"Missing variable in input file: {v}")

        # Identify omega variable
        if "omega" in ds_input:
            omega_name = "omega"
            ds = ds_input[
                [
                    "calipso_lcc",
                    "modis_lcc",
                    "ts",
                    "eis",
                    "sic",
                    "tadv",
                    omega_name,
                ]
            ].load()

        elif "omega700" in ds_input:
            omega_name = "omega700"
            ds = ds_input[
                [
                    "calipso_lcc",
                    "modis_lcc",
                    "ts",
                    "eis",
                    "sic",
                    "tadv",
                    omega_name,
                ]
            ].load()

        else:
            print("[Log] No omega/omega700 in unified file. Loading omega700 from merged file.")

            with xr.open_dataset(file_path_omega) as ds_omega:
                if "omega700" not in ds_omega:
                    raise ValueError("omega700 is also missing in the merged file.")

                ds_main = ds_input[
                    [
                        "calipso_lcc",
                        "modis_lcc",
                        "ts",
                        "eis",
                        "sic",
                        "tadv",
                    ]
                ]

                ds_omg = ds_omega[["omega700"]]

                ds_main, ds_omg = xr.align(ds_main, ds_omg, join="inner")
                ds = xr.merge([ds_main, ds_omg]).load()

                omega_name = "omega700"

        # Ocean / valid-area mask
        ocean_mask = ds["calipso_lcc"].mean(dim="time").notnull()
        ds_ocean = ds.where(ocean_mask)

        # Unit harmonization
        if float(ds_ocean["ts"].max()) > 200:
            ds_ocean["ts"] = ds_ocean["ts"] - 273.15

        if float(ds_ocean["sic"].max()) <= 1.1:
            ds_ocean["sic"] = ds_ocean["sic"] * 100.0

        # ==================================================
        # 5. Global standardized regression
        # ==================================================
        print("[Log] Calculating global standardized regression coefficients...")

        ds_stand = xr.Dataset()

        for var in [
            "calipso_lcc",
            "modis_lcc",
            "ts",
            "eis",
            omega_name,
            "sic",
            "tadv",
        ]:
            ds_stand[var] = standardize_robust_xarray(
                deseasonalize_robust(ds_ocean[var])
            )

        coeff_names = [
            "ts",
            "eis",
            "omega",
            "sic",
            "tadv",
        ]

        coeffs_cal = run_global_regression(
            ds_stand["calipso_lcc"],
            ds_stand["ts"],
            ds_stand["eis"],
            ds_stand[omega_name],
            ds_stand["sic"],
            ds_stand["tadv"],
        ).assign_coords(coefficient=coeff_names)

        coeffs_mod = run_global_regression(
            ds_stand["modis_lcc"],
            ds_stand["ts"],
            ds_stand["eis"],
            ds_stand[omega_name],
            ds_stand["sic"],
            ds_stand["tadv"],
        ).assign_coords(coefficient=coeff_names)

        # ==================================================
        # 6. Ts-EIS state-space physical-unit regression
        # ==================================================
        print("[Log] Calculating Ts-EIS state-space physical sensitivities...")

        ds_anom = xr.Dataset()

        for var in [
            "calipso_lcc",
            "modis_lcc",
            "ts",
            "eis",
            omega_name,
            "sic",
            "tadv",
        ]:
            ds_anom[var] = deseasonalize_robust(ds_ocean[var])

        ds_final = xr.Dataset()

        ds_final["ts_total"] = ds_ocean["ts"]
        ds_final["eis_total"] = ds_ocean["eis"]
        ds_final["sic_total"] = ds_ocean["sic"]

        ds_final["cal_lcc_anom"] = ds_anom["calipso_lcc"] * 100.0
        ds_final["mod_lcc_anom"] = ds_anom["modis_lcc"] * 100.0

        ds_final["ts_anom"] = ds_anom["ts"]
        ds_final["eis_anom"] = ds_anom["eis"]
        ds_final["omega_anom"] = ds_anom[omega_name]
        ds_final["sic_anom"] = ds_anom["sic"]
        ds_final["tadv_anom"] = ds_anom["tadv"]

        print("[Log] Converting state-space data to DataFrame...")

        df_main = ds_final.stack(
            point=("time", "latitude", "longitude")
        ).to_dataframe()

        df_main.dropna(
            subset=[
                "cal_lcc_anom",
                "mod_lcc_anom",
                "ts_total",
                "eis_total",
                "sic_total",
                "ts_anom",
                "eis_anom",
                "omega_anom",
                "tadv_anom",
            ],
            inplace=True,
        )

        print(f"[Log] Valid samples for state-space regression: {len(df_main)}")

        ts_bins = np.arange(-31, 32, 2)
        eis_bins = np.arange(-5, 26, 2)

        ts_centers = (ts_bins[:-1] + ts_bins[1:]) / 2
        eis_centers = (eis_bins[:-1] + eis_bins[1:]) / 2

        df_main["ts_bin"] = pd.cut(
            df_main["ts_total"],
            bins=ts_bins,
            labels=ts_centers,
            right=False,
        )

        df_main["eis_bin"] = pd.cut(
            df_main["eis_total"],
            bins=eis_bins,
            labels=eis_centers,
            right=False,
        )

        print("[Log] Running grouped physical-unit regression...")

        sensitivity = (
            df_main.groupby(["eis_bin", "ts_bin"], observed=False)
            .apply(get_state_space_physical_sensitivities, include_groups=False)
        )

        grids = {}

        for metric in [
            "ts",
            "eis",
            "omega",
            "tadv",
            "sic",
        ]:
            grids[f"cal_{metric}"] = sensitivity[f"cal_beta_{metric}"].unstack()
            grids[f"mod_{metric}"] = sensitivity[f"mod_beta_{metric}"].unstack()

        # ==================================================
        # 7. Export state-space physical-sensitivity CSV
        # ==================================================
        print("[Log] Exporting state-space physical-sensitivity CSV...")

        df_csv = pd.DataFrame({
            "CALIPSO_partial_LCC_partial_Ts_percent_per_K": grids["cal_ts"].stack(dropna=False),
            "MODIS_partial_LCC_partial_Ts_percent_per_K": grids["mod_ts"].stack(dropna=False),

            "CALIPSO_partial_LCC_partial_EIS_percent_per_K": grids["cal_eis"].stack(dropna=False),
            "MODIS_partial_LCC_partial_EIS_percent_per_K": grids["mod_eis"].stack(dropna=False),

            "CALIPSO_partial_LCC_partial_omega_percent_per_hPa_day": grids["cal_omega"].stack(dropna=False),
            "MODIS_partial_LCC_partial_omega_percent_per_hPa_day": grids["mod_omega"].stack(dropna=False),

            "CALIPSO_partial_LCC_partial_Tadv_percent_per_K_day": grids["cal_tadv"].stack(dropna=False),
            "MODIS_partial_LCC_partial_Tadv_percent_per_K_day": grids["mod_tadv"].stack(dropna=False),

            "CALIPSO_partial_LCC_partial_SIC_percent_per_percent": grids["cal_sic"].stack(dropna=False),
            "MODIS_partial_LCC_partial_SIC_percent_per_percent": grids["mod_sic"].stack(dropna=False),
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

        print(f"[Success] CSV saved to: {csv_save_path}")

    # ======================================================
    # 8. Plotting: 5 rows × 4 columns
    # ======================================================
    print("[Log] Plotting combined 5-variable figure...")

    proj = ccrs.Robinson(central_longitude=180)

    fig = plt.figure(figsize=(31.5, 30.8))

    outer_gs = fig.add_gridspec(
        nrows=1,
        ncols=2,
        left=0.035,
        right=0.905,
        bottom=0.070,
        top=0.965,
        width_ratios=[1.06, 1.00],
        wspace=0.16,
    )

    left_gs = outer_gs[0, 0].subgridspec(
        nrows=5,
        ncols=2,
        wspace=0.10,
        hspace=0.19,
    )

    right_gs = outer_gs[0, 1].subgridspec(
        nrows=5,
        ncols=2,
        wspace=0.12,
        hspace=0.19,
    )

    axes = np.empty((5, 4), dtype=object)

    for r in range(5):
        axes[r, 0] = fig.add_subplot(left_gs[r, 0], projection=proj)
        axes[r, 1] = fig.add_subplot(left_gs[r, 1], projection=proj)

        axes[r, 2] = fig.add_subplot(right_gs[r, 0])
        axes[r, 3] = fig.add_subplot(right_gs[r, 1])

    variables = [
        {"key": "ts", "title": "∂LCC/∂Ts"},
        {"key": "eis", "title": "∂LCC/∂EIS"},
        {"key": "omega", "title": "∂LCC/∂ω"},
        {"key": "tadv", "title": r"$\mathrm{\partial LCC/\partial T_{adv}}$"},
        {"key": "sic", "title": "∂LCC/∂SIC"},
    ]

    panel_labels = [
        "(a)", "(b)", "(c)", "(d)",
        "(e)", "(f)", "(g)", "(h)",
        "(i)", "(j)", "(k)", "(l)",
        "(m)", "(n)", "(o)", "(p)",
        "(q)", "(r)", "(s)", "(t)",
    ]

    panel_idx = 0

    global_mesh = None
    state_meshes = {}

    for r, var_info in enumerate(variables):
        key = var_info["key"]
        title = var_info["title"]

        global_mesh = plot_global_map(
            axes[r, 0],
            coeffs_cal.sel(coefficient=key),
            f"CALIPSO Global {title}",
            panel_labels[panel_idx],
        )
        panel_idx += 1

        global_mesh = plot_global_map(
            axes[r, 1],
            coeffs_mod.sel(coefficient=key),
            f"MODIS Global {title}",
            panel_labels[panel_idx],
        )
        panel_idx += 1

        state_mesh_cal = plot_state_space(
            axes[r, 2],
            grids[f"cal_{key}"],
            f"CALIPSO State-space {title}",
            panel_labels[panel_idx],
            ts_bins,
            eis_bins,
            state_vlims[key],
        )
        panel_idx += 1

        state_mesh_mod = plot_state_space(
            axes[r, 3],
            grids[f"mod_{key}"],
            f"MODIS State-space {title}",
            panel_labels[panel_idx],
            ts_bins,
            eis_bins,
            state_vlims[key],
        )
        panel_idx += 1

        state_meshes[key] = state_mesh_mod

    # ======================================================
    # 9. Axis labels
    # ======================================================
    for r in range(5):
        for c in range(4):

            if c < 2:
                continue

            if c == 2:
                axes[r, c].set_ylabel(
                    "EIS (K)",
                    fontsize=label_fontsize,
                    fontweight=font_weight,
                    fontfamily="Times New Roman",
                    labelpad=5,
                )
            else:
                axes[r, c].set_ylabel("")

            if r == 4:
                axes[r, c].set_xlabel(
                    "Ts (°C)",
                    fontsize=label_fontsize,
                    fontweight=font_weight,
                    fontfamily="Times New Roman",
                    labelpad=5,
                )
            else:
                axes[r, c].set_xlabel("")

    # ======================================================
    # 10. Colorbars
    # ======================================================

    # ------------------------------------------------------
    # 10.1 Shared colorbar for global maps, columns 1-2
    # ------------------------------------------------------
    pos_left0 = axes[4, 0].get_position()
    pos_left1 = axes[4, 1].get_position()

    cbar_ax_global = fig.add_axes([
        pos_left0.x0 + 0.015,
        0.032,
        pos_left1.x1 - pos_left0.x0 - 0.030,
        0.018,
    ])

    cb_global = fig.colorbar(
        global_mesh,
        cax=cbar_ax_global,
        orientation="horizontal",
        extend="both",
        ticks=np.arange(-1.0, 1.01, 0.5),
    )

    cb_global.set_label(
        "Global standardized coefficient",
        fontsize=cbar_label_fontsize,
        fontweight=font_weight,
        fontfamily="Times New Roman",
        labelpad=7,
    )

    cb_global.ax.tick_params(labelsize=cbar_tick_fontsize, width=1.8, length=6)

    for t in cb_global.ax.get_xticklabels():
        t.set_fontweight(font_weight)
        t.set_fontfamily("Times New Roman")

    # ------------------------------------------------------
    # 10.2 One state-space colorbar for each row, columns 3-4
    # ------------------------------------------------------
    for r, key in enumerate(["ts", "eis", "omega", "tadv", "sic"]):

        pos = axes[r, 3].get_position()

        cbar_height = pos.height * 0.82
        cbar_bottom = pos.y0 + (pos.height - cbar_height) / 2

        cbar_ax = fig.add_axes([
            pos.x1 + 0.008,
            cbar_bottom,
            0.012,
            cbar_height,
        ])

        ticks = [-state_vlims[key], 0, state_vlims[key]]

        cb = fig.colorbar(
            state_meshes[key],
            cax=cbar_ax,
            orientation="vertical",
            extend="both",
            ticks=ticks,
        )

        cb.set_label(
            state_units[key],
            fontsize=cbar_label_fontsize - 2,
            fontweight=font_weight,
            fontfamily="Times New Roman",
            labelpad=6,
        )

        cb.ax.tick_params(labelsize=cbar_tick_fontsize, width=1.6, length=5)

        for t in cb.ax.get_yticklabels():
            t.set_fontweight(font_weight)
            t.set_fontfamily("Times New Roman")

    # ======================================================
    # 11. Save figure
    # ======================================================
    save_path = os.path.join(output_dir, output_filename)

    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"[Success] Combined 5-variable figure saved to: {save_path}")

    plt.show()

except Exception as e:
    print(f"[Error] {e}")
    traceback.print_exc()
