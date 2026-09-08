"""Generate manuscript Figure 3 from monthly fields aggregated to 10° × 10°.

The native 2° latitude × 5° longitude monthly fields are first averaged to
10° × 10° (5 latitude cells × 2 longitude cells). Global standardized
grid-box regressions and unstandardized Ts–EIS-bin regressions are then
applied to the aggregated monthly fields.
"""

from __future__ import annotations

import os
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr


PROJECT_DIR = Path(__file__).resolve().parents[2]
INPUT_FILE = PROJECT_DIR / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2degx5deg_monthly_2007-2021.nc"
OUTPUT_DIR = PROJECT_DIR / "outputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

VARIABLES = ["ts", "eis", "omega", "tadv", "sic"]
GLOBAL_COEFFICIENT_ORDER = ["ts", "eis", "omega", "sic", "tadv"]
COEFFICIENT_LABELS = {
    "ts": r"$\partial\mathrm{LCC}/\partial \mathrm{Ts}$",
    "eis": r"$\partial\mathrm{LCC}/\partial\mathrm{EIS}$",
    "omega": r"$\partial\mathrm{LCC}/\partial\omega$",
    "tadv": r"$\partial\mathrm{LCC}/\partial T_{adv}$",
    "sic": r"$\partial\mathrm{LCC}/\partial\mathrm{SIC}$",
}
STATE_LIMITS = {"ts": 5.0, "eis": 5.0, "omega": 0.20, "tadv": 5.0, "sic": 1.0}
STATE_UNITS = {
    "ts": "% / K",
    "eis": "% / K",
    "omega": "% / (hPa day$^{-1}$)",
    "tadv": "% / (K day$^{-1}$)",
    "sic": "% / %",
}
# Coarsening the monthly fields to 10° × 10° reduces the number of
# spatiotemporal samples in each Ts–EIS bin.  A threshold of 150 retains a
# reasonably sampled fit while avoiding the unnecessarily sparse coverage
# produced by the native-grid threshold of 300.
MIN_STATE_SAMPLES = 150


def harmonize(dataset: xr.Dataset) -> xr.Dataset:
    names = ["calipso_lcc", "modis_lcc", "ts", "eis", "omega", "sic", "tadv"]
    ds = dataset[names].copy()
    for name in names:
        ds[name] = ds[name].transpose("time", "latitude", "longitude")
    if float(ds["ts"].max()) > 200:
        ds["ts"] = ds["ts"] - 273.15
    if float(ds["sic"].max()) <= 1.1:
        ds["sic"] = ds["sic"] * 100.0
    ocean_mask = ds["calipso_lcc"].mean("time", skipna=True).notnull()
    return ds.where(ocean_mask)


def deseasonalize(data: xr.DataArray) -> xr.DataArray:
    return data.groupby("time.month") - data.groupby("time.month").mean("time", skipna=True)


def standardize(data: xr.DataArray) -> xr.DataArray:
    mean = data.mean("time", skipna=True)
    std = data.std("time", skipna=True)
    return ((data - mean) / std).where(std > 1e-6)


def regress_one_grid(y, predictors):
    core = np.all(np.isfinite(predictors[:, [0, 1, 2, 4]]), axis=1) & np.isfinite(y)
    if core.sum() < 20:
        return np.full(5, np.nan)
    yy = y[core]
    xx = predictors[core]
    use_sic = np.all(np.isfinite(xx[:, 3])) and np.nanstd(xx[:, 3]) > 1e-4
    columns = [0, 1, 2, 3, 4] if use_sic else [0, 1, 2, 4]
    try:
        design = sm.add_constant(xx[:, columns], prepend=True)
        params = sm.OLS(yy, design).fit().params[1:]
    except Exception:
        return np.full(5, np.nan)
    result = np.full(5, np.nan)
    result[columns] = params
    return result


def global_coefficients(ds: xr.Dataset) -> tuple[xr.DataArray, xr.DataArray]:
    anomalies = {name: standardize(deseasonalize(ds[name])) for name in ds.data_vars}
    predictors = np.stack(
        [anomalies[name].values for name in ["ts", "eis", "omega", "sic", "tadv"]],
        axis=-1,
    )
    nlat, nlon = ds.sizes["latitude"], ds.sizes["longitude"]

    outputs = []
    for sensor in ["calipso_lcc", "modis_lcc"]:
        y = anomalies[sensor].values
        beta = np.full((nlat, nlon, 5), np.nan)
        for ilat in range(nlat):
            for ilon in range(nlon):
                beta[ilat, ilon, :] = regress_one_grid(y[:, ilat, ilon], predictors[:, ilat, ilon, :])
        outputs.append(
            xr.DataArray(
                beta,
                dims=("latitude", "longitude", "coefficient"),
                coords={
                    "latitude": ds.latitude,
                    "longitude": ds.longitude,
                    "coefficient": GLOBAL_COEFFICIENT_ORDER,
                },
            )
        )
    return outputs[0], outputs[1]


def aggregate_to_10deg(ds: xr.Dataset) -> xr.Dataset:
    # Native spacing is 2° latitude x 5° longitude; 5 x 2 cells = 10° x 10°.
    return ds.coarsen(latitude=5, longitude=2, boundary="trim").mean(skipna=True)


def fit_state_bin(frame: pd.DataFrame, sensor: str, use_sic: bool):
    x_names = ["ts_anom", "eis_anom", "omega_hpa_day", "tadv_k_day"]
    if use_sic:
        x_names.append("sic_anom")
    regression = frame[[sensor] + x_names].dropna()
    if len(regression) < MIN_STATE_SAMPLES:
        return np.full(5, np.nan), np.full(5, np.nan)
    model = sm.OLS(regression[sensor], sm.add_constant(regression[x_names], prepend=True)).fit()
    result = np.full(5, np.nan)
    p_values = np.full(5, np.nan)
    result[0] = model.params.get("ts_anom", np.nan)
    result[1] = model.params.get("eis_anom", np.nan)
    result[2] = model.params.get("omega_hpa_day", np.nan)
    result[3] = model.params.get("tadv_k_day", np.nan)
    p_values[0] = model.pvalues.get("ts_anom", np.nan)
    p_values[1] = model.pvalues.get("eis_anom", np.nan)
    p_values[2] = model.pvalues.get("omega_hpa_day", np.nan)
    p_values[3] = model.pvalues.get("tadv_k_day", np.nan)
    if use_sic:
        result[4] = model.params.get("sic_anom", np.nan)
        p_values[4] = model.pvalues.get("sic_anom", np.nan)
    return result, p_values


def state_space_coefficients(ds: xr.Dataset):
    anomalies = {name: deseasonalize(ds[name]) for name in ds.data_vars}
    frame = pd.DataFrame(
        {
            "ts_total": ds["ts"].values.ravel(),
            "eis_total": ds["eis"].values.ravel(),
            "sic_total": ds["sic"].values.ravel(),
            "cal": (anomalies["calipso_lcc"].values * 100.0).ravel(),
            "mod": (anomalies["modis_lcc"].values * 100.0).ravel(),
            "ts_anom": anomalies["ts"].values.ravel(),
            "eis_anom": anomalies["eis"].values.ravel(),
            "omega_hpa_day": (anomalies["omega"].values * 864.0).ravel(),
            "sic_anom": anomalies["sic"].values.ravel(),
            "tadv_k_day": (anomalies["tadv"].values * 86400.0).ravel(),
        }
    )
    required = ["ts_total", "eis_total", "sic_total", "cal", "mod", "ts_anom", "eis_anom", "omega_hpa_day", "tadv_k_day"]
    frame = frame.dropna(subset=required)

    ts_edges = np.arange(-31, 32, 2)
    eis_edges = np.arange(-5, 26, 2)
    ts_centers = (ts_edges[:-1] + ts_edges[1:]) / 2
    eis_centers = (eis_edges[:-1] + eis_edges[1:]) / 2
    frame["ts_bin"] = pd.cut(frame.ts_total, ts_edges, labels=False, right=False)
    frame["eis_bin"] = pd.cut(frame.eis_total, eis_edges, labels=False, right=False)
    frame = frame.dropna(subset=["ts_bin", "eis_bin"])
    frame[["ts_bin", "eis_bin"]] = frame[["ts_bin", "eis_bin"]].astype(int)

    cal = np.full((5, len(eis_centers), len(ts_centers)), np.nan)
    mod = np.full_like(cal, np.nan)
    cal_p = np.full_like(cal, np.nan)
    mod_p = np.full_like(cal, np.nan)
    for (ieis, its), group in frame.groupby(["eis_bin", "ts_bin"], sort=False):
        if len(group) < MIN_STATE_SAMPLES:
            continue
        saturated = ((group.sic_total <= 1.0) | (group.sic_total >= 99.0)).sum() / len(group)
        use_sic = saturated <= 0.5 and group.sic_anom.std() >= 0.1
        cal[:, ieis, its], cal_p[:, ieis, its] = fit_state_bin(group, "cal", use_sic)
        mod[:, ieis, its], mod_p[:, ieis, its] = fit_state_bin(group, "mod", use_sic)
    return cal, mod, cal_p, mod_p, ts_edges, eis_edges


def coordinate_edges(centers):
    centers = np.asarray(centers, dtype=float)
    delta = np.diff(centers).mean()
    return np.concatenate(([centers[0] - delta / 2], centers + delta / 2))


def map_panel(ax, data: xr.DataArray, title: str, panel: str, vmin=-1.0, vmax=1.0, cmap="RdBu_r"):
    cyclic, cyclic_lon = add_cyclic_point(data.values, coord=data.longitude.values)
    mesh = ax.pcolormesh(
        coordinate_edges(cyclic_lon),
        coordinate_edges(data.latitude.values),
        cyclic,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
    )
    ax.add_feature(cfeature.LAND, facecolor="0.72", edgecolor="none", zorder=3)
    ax.coastlines(linewidth=0.55, zorder=4)
    ax.gridlines(draw_labels=False, linewidth=0.45, color="gray", alpha=0.55, linestyle="--", zorder=4)
    ax.set_global()
    # Preserve the taller map-panel geometry used by the original Figure 3.
    ax.set_aspect(1.5)
    ax.set_title(f"{panel} {title}", fontsize=10.5, fontweight="bold", pad=3)
    return mesh


def state_panel(ax, values, ts_edges, eis_edges, title, panel, limit, p_values=None):
    mesh = ax.pcolormesh(ts_edges, eis_edges, values, cmap="RdBu_r", vmin=-limit, vmax=limit, shading="flat")
    ax.set_xlim(-30, 30)
    ax.set_ylim(-5, 25)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.tick_params(labelsize=23, width=1.8, length=6)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
        label.set_fontfamily("Times New Roman")
    ax.grid(linestyle="--", linewidth=0.4, alpha=0.35)
    if p_values is not None:
        significant = np.where((p_values < 0.05) & np.isfinite(values))
        ts_centers = (ts_edges[:-1] + ts_edges[1:]) / 2
        eis_centers = (eis_edges[:-1] + eis_edges[1:]) / 2
        ax.scatter(ts_centers[significant[1]], eis_centers[significant[0]], s=10, c="black", marker=".", zorder=3)
    ax.set_title(f"{panel} {title}", fontsize=10.5, fontweight="bold", pad=3)
    return mesh


def plot_full_figure(cal10, mod10, state_cal, state_mod, state_cal_p, state_mod_p, ts_edges, eis_edges, output_stem):
    projection = ccrs.Robinson(central_longitude=180)
    plt.rcParams.update({"font.family": "Times New Roman", "font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold"})
    fig = plt.figure(figsize=(31.5, 30.8))
    outer = fig.add_gridspec(1, 2, left=0.035, right=0.905, bottom=0.070, top=0.965, width_ratios=[1.06, 1.00], wspace=0.16)
    left = outer[0, 0].subgridspec(5, 2, wspace=0.10, hspace=0.19)
    right = outer[0, 1].subgridspec(5, 2, wspace=0.12, hspace=0.19)
    labels = [f"({chr(97 + i)})" for i in range(20)]
    global_mesh = None
    axes = np.empty((5, 4), dtype=object)
    state_meshes = []
    for row, key in enumerate(VARIABLES):
        ax0 = fig.add_subplot(left[row, 0], projection=projection)
        ax1 = fig.add_subplot(left[row, 1], projection=projection)
        ax2 = fig.add_subplot(right[row, 0])
        ax3 = fig.add_subplot(right[row, 1])
        axes[row] = [ax0, ax1, ax2, ax3]
        title = COEFFICIENT_LABELS[key]
        global_mesh = map_panel(ax0, cal10.sel(coefficient=key), f"CALIPSO Global {title}", labels[row * 4])
        global_mesh = map_panel(ax1, mod10.sel(coefficient=key), f"MODIS Global {title}", labels[row * 4 + 1])
        state_mesh = state_panel(ax2, state_cal[row], ts_edges, eis_edges, f"CALIPSO State-space {title}", labels[row * 4 + 2], STATE_LIMITS[key], p_values=state_cal_p[row])
        state_mesh = state_panel(ax3, state_mod[row], ts_edges, eis_edges, f"MODIS State-space {title}", labels[row * 4 + 3], STATE_LIMITS[key], p_values=state_mod_p[row])
        for axis in (ax0, ax1, ax2, ax3):
            axis.title.set_fontsize(25)
            axis.title.set_fontweight("bold")
            axis.title.set_fontfamily("Times New Roman")
        ax2.set_ylabel("EIS (K)", fontsize=26, fontweight="bold", fontfamily="Times New Roman", labelpad=5)
        if row == 4:
            ax2.set_xlabel(r"$T_s$ (°C)", fontsize=9, fontweight="bold")
            ax3.set_xlabel(r"$T_s$ (°C)", fontsize=9, fontweight="bold")
        pos = ax3.get_position()
        cbar_height = pos.height * 0.82
        cax = fig.add_axes([pos.x1 + 0.008, pos.y0 + (pos.height - cbar_height) / 2, 0.012, cbar_height])
        cbar = fig.colorbar(
            state_mesh,
            cax=cax,
            extend="both",
            ticks=np.linspace(-STATE_LIMITS[key], STATE_LIMITS[key], 5),
        )
        cbar.set_label(STATE_UNITS[key], fontsize=23, fontweight="bold", fontfamily="Times New Roman", labelpad=6)
        cbar.ax.tick_params(labelsize=25, width=1.6, length=5)
        for label in cbar.ax.get_yticklabels():
            label.set_fontweight("bold")
            label.set_fontfamily("Times New Roman")
        state_meshes.append(state_mesh)
        if row == 4:
            ax2.set_xlabel("Ts (°C)", fontsize=26, fontweight="bold", fontfamily="Times New Roman", labelpad=5)
            ax3.set_xlabel("Ts (°C)", fontsize=26, fontweight="bold", fontfamily="Times New Roman", labelpad=5)
    pos0, pos1 = axes[4, 0].get_position(), axes[4, 1].get_position()
    cax = fig.add_axes([pos0.x0 + 0.015, 0.032, pos1.x1 - pos0.x0 - 0.030, 0.018])
    cbar = fig.colorbar(global_mesh, cax=cax, orientation="horizontal", extend="both", ticks=np.arange(-1.0, 1.01, 0.5))
    cbar.set_label("Global standardized coefficient", fontsize=25, fontweight="bold", fontfamily="Times New Roman")
    cbar.ax.tick_params(labelsize=25, width=1.8, length=6)
    for label in cbar.ax.get_xticklabels():
        label.set_fontweight("bold")
        label.set_fontfamily("Times New Roman")
    fig.savefig(OUTPUT_DIR / f"{output_stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT_DIR / f"{output_stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_comparison(cal2, mod2, cal10, mod10):
    projection = ccrs.Robinson(central_longitude=180)
    fig = plt.figure(figsize=(16, 13))
    grid = fig.add_gridspec(5, 4, left=0.06, right=0.985, bottom=0.08, top=0.94, wspace=0.08, hspace=0.17)
    coefficient_axes = []
    coefficient_mesh = None
    for row, key in enumerate(VARIABLES):
        datasets = [
            cal2.sel(coefficient=key),
            cal10.sel(coefficient=key),
            mod2.sel(coefficient=key),
            mod10.sel(coefficient=key),
        ]
        titles = [
            "CALIPSO 2°×5° regression", "CALIPSO 10°×10° regression",
            "MODIS 2°×5° regression", "MODIS 10°×10° regression",
        ]
        for col, (data, title) in enumerate(zip(datasets, titles)):
            ax = fig.add_subplot(grid[row, col], projection=projection)
            mesh = map_panel(ax, data, "", "", -1.0, 1.0)
            ax.set_title(title if row == 0 else "", fontsize=9.5, fontweight="bold", pad=4)
            if col == 0:
                ax.text(-0.08, 0.5, COEFFICIENT_LABELS[key], transform=ax.transAxes, rotation=90,
                        ha="center", va="center", fontsize=10, fontweight="bold")
            coefficient_axes.append(ax)
            coefficient_mesh = mesh
    cbar = fig.colorbar(coefficient_mesh, ax=coefficient_axes, orientation="horizontal", fraction=0.025, pad=0.025, extend="both")
    cbar.set_label("Standardized partial regression coefficient", fontweight="bold")
    fig.suptitle("Global regressions: native 2°×5° data vs. monthly fields aggregated to 10°×10° before regression", fontsize=13, fontweight="bold", y=0.985)
    fig.savefig(OUTPUT_DIR / "Figure3_2x5_vs_10x10_global_comparison.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_global_comparison_with_difference(cal2, mod2, cal10, mod10):
    """Add common-grid differences without relabeling them as native results."""
    projection = ccrs.Robinson(central_longitude=180)
    fig = plt.figure(figsize=(22, 13))
    grid = fig.add_gridspec(5, 6, left=0.045, right=0.985, bottom=0.10, top=0.93, wspace=0.08, hspace=0.17)
    coefficient_axes = []
    difference_axes = []
    coefficient_mesh = None
    difference_mesh = None
    cal2_block = cal2.coarsen(latitude=5, longitude=2, boundary="trim").mean(skipna=True)
    mod2_block = mod2.coarsen(latitude=5, longitude=2, boundary="trim").mean(skipna=True)
    for row, key in enumerate(VARIABLES):
        datasets = [
            cal2.sel(coefficient=key),
            cal10.sel(coefficient=key),
            cal10.sel(coefficient=key) - cal2_block.sel(coefficient=key),
            mod2.sel(coefficient=key),
            mod10.sel(coefficient=key),
            mod10.sel(coefficient=key) - mod2_block.sel(coefficient=key),
        ]
        titles = [
            "CALIPSO 2°×5° regression", "CALIPSO 10°×10° regression", "CALIPSO difference",
            "MODIS 2°×5° regression", "MODIS 10°×10° regression", "MODIS difference",
        ]
        for col, (data, title) in enumerate(zip(datasets, titles)):
            ax = fig.add_subplot(grid[row, col], projection=projection)
            is_difference = col in (2, 5)
            limit = 0.5 if is_difference else 1.0
            mesh = map_panel(ax, data, "", "", -limit, limit)
            ax.set_title(title if row == 0 else "", fontsize=9.2, fontweight="bold", pad=4)
            if col == 0:
                ax.text(-0.08, 0.5, COEFFICIENT_LABELS[key], transform=ax.transAxes, rotation=90,
                        ha="center", va="center", fontsize=10, fontweight="bold")
            if is_difference:
                difference_axes.append(ax)
                difference_mesh = mesh
            else:
                coefficient_axes.append(ax)
                coefficient_mesh = mesh
    cbar1 = fig.colorbar(coefficient_mesh, ax=coefficient_axes, orientation="horizontal", fraction=0.022, pad=0.025, extend="both")
    cbar1.set_label("Standardized partial regression coefficient", fontweight="bold")
    cbar2 = fig.colorbar(difference_mesh, ax=difference_axes, orientation="horizontal", fraction=0.022, pad=0.025, extend="both")
    cbar2.set_label("Difference: 10°×10° regression minus mean 2°×5° coefficient in the same block", fontweight="bold")
    fig.suptitle("Resolution sensitivity of global regressions", fontsize=14, fontweight="bold", y=0.985)
    fig.savefig(OUTPUT_DIR / "Figure3_2x5_vs_10x10_global_comparison_with_difference.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_state_comparison(state_cal2, state_mod2, state_cal10, state_mod10, ts_edges, eis_edges):
    fig, axes = plt.subplots(5, 4, figsize=(15, 16), constrained_layout=True)
    titles = [
        "CALIPSO 2°×5° samples", "CALIPSO 10°×10° samples",
        "MODIS 2°×5° samples", "MODIS 10°×10° samples",
    ]
    for row, key in enumerate(VARIABLES):
        fields = [state_cal2[row], state_cal10[row], state_mod2[row], state_mod10[row]]
        row_mesh = None
        for col, (field, title) in enumerate(zip(fields, titles)):
            row_mesh = state_panel(
                axes[row, col], field, ts_edges, eis_edges,
                title if row == 0 else "", "", STATE_LIMITS[key]
            )
            if col == 0:
                axes[row, col].set_ylabel(f"{COEFFICIENT_LABELS[key]}\nEIS (K)", fontsize=9, fontweight="bold")
            if row == 4:
                axes[row, col].set_xlabel(r"$T_s$ (°C)", fontsize=9, fontweight="bold")
        cbar = fig.colorbar(row_mesh, ax=axes[row, :], fraction=0.018, pad=0.01, extend="both")
        cbar.set_label(STATE_UNITS[key], fontsize=8.5, fontweight="bold")
    fig.suptitle("Ts–EIS-bin regressions after changing only the spatial grid", fontsize=14, fontweight="bold")
    fig.savefig(OUTPUT_DIR / "Figure3_2x5_vs_10x10_state_space_comparison.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_state_differences(state_cal2, state_mod2, state_cal10, state_mod10, ts_edges, eis_edges):
    fig, axes = plt.subplots(5, 2, figsize=(8.5, 16), constrained_layout=True)
    for row, key in enumerate(VARIABLES):
        fields = [state_cal10[row] - state_cal2[row], state_mod10[row] - state_mod2[row]]
        row_mesh = None
        for col, (field, sensor) in enumerate(zip(fields, ["CALIPSO", "MODIS"])):
            row_mesh = state_panel(
                axes[row, col], field, ts_edges, eis_edges,
                f"{sensor} difference" if row == 0 else "", "", STATE_LIMITS[key]
            )
            if col == 0:
                axes[row, col].set_ylabel(f"{COEFFICIENT_LABELS[key]}\nEIS (K)", fontsize=9, fontweight="bold")
            if row == 4:
                axes[row, col].set_xlabel(r"$T_s$ (°C)", fontsize=9, fontweight="bold")
        cbar = fig.colorbar(row_mesh, ax=axes[row, :], fraction=0.025, pad=0.015, extend="both")
        cbar.set_label(STATE_UNITS[key], fontsize=8.5, fontweight="bold")
    fig.suptitle("State-space difference: 10°×10°-sample regression minus 2°×5°-sample regression", fontsize=12, fontweight="bold")
    fig.savefig(OUTPUT_DIR / "Figure3_2x5_vs_10x10_state_space_difference.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    print("Loading native monthly fields")
    with xr.open_dataset(INPUT_FILE) as source:
        native = harmonize(source).load()
    print("Aggregating monthly fields to 10° × 10°")
    coarse = aggregate_to_10deg(native)
    print("Computing 10° × 10° global coefficients")
    cal10, mod10 = global_coefficients(coarse)
    print("Computing 10° × 10° Ts–EIS state-space sensitivities")
    state_cal10, state_mod10, state_cal10_p, state_mod10_p, ts_edges, eis_edges = state_space_coefficients(coarse)
    plot_full_figure(
        cal10, mod10, state_cal10, state_mod10, state_cal10_p, state_mod10_p, ts_edges, eis_edges,
        "Figure3_LCC_Sensitivity_Global_StateSpace_10degx10deg",
    )
    print(f"Outputs written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
