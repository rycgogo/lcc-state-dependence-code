import os
import numpy as np
import pandas as pd
import xarray as xr
from pyhdf.SD import SD, SDC
from pathlib import Path


# ==============================================================================
# 最终任务: 计算MODIS C_low并以正确的纬度顺序保存为CALIPSO格式
# ==============================================================================

def find_all_hdf_files(root_folder):
    """遍历所有子目录，找到所有.hdf文件。"""
    hdf_files = []
    for dirpath, _, filenames in os.walk(root_folder):
        for f in filenames:
            if f.lower().endswith('.hdf'):
                hdf_files.append(os.path.join(dirpath, f))
    return sorted(hdf_files)


def get_low_cloud_indices(sample_file_path, standard_pressure_intervals_var, threshold):
    """根据标准压力区间变量，计算出哪些索引属于低云。"""
    try:
        hdf = SD(sample_file_path, SDC.READ)
        intervals = hdf.select(standard_pressure_intervals_var).get()
        hdf.end()
        bin_centers = (intervals[:-1] + intervals[1:]) / 2.0
        low_cloud_indices = np.where(bin_centers > threshold)[0]
        print("--- 低云分箱定义 ---")
        print(f"使用 '{standard_pressure_intervals_var}' 作为压力基准。")
        print(f"使用科学阈值: CTP > {threshold} hPa")
        print(f"计算出的标准低云分箱索引为: {low_cloud_indices}")
        return low_cloud_indices, True
    except Exception as e:
        print(f"错误: 无法读取标准压力区间变量 -> {e}")
        return None, False


def calculate_and_save_modis_clow_correct_order():
    """主函数，以正确的纬度顺序计算C_low并将其重采样至CALIPSO网格后保存。"""
    print("=" * 50)
    print("开始计算MODIS C_low (V12.2 - 随机重叠算法)...")
    print("=" * 50)

    # --- 用户配置区 ---
    project_dir = Path(__file__).resolve().parents[2]
    modis_data_root = project_dir / "raw_inputs" / "modis"
    # <-- [修改] 已更新为正确的文件名
    calipso_template_file = project_dir / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc"
    start_year = 2007
    end_year = 2021
    output_filename = project_dir / "derived_intermediate" / "modis_lcc_2.5deg_monthly_2007-2021.nc"
    output_filename.parent.mkdir(parents=True, exist_ok=True)

    # --- 变量定义 ---
    total_cloud_fraction_var = 'Cloud_Fraction_Mean_Mean'
    pressure_histogram_var = 'Cloud_Top_Pressure_Histogram_Counts'
    standard_pressure_intervals_var = 'Cloud_Top_Pressure_Histo_Intervals'
    low_cloud_threshold = 680
    CF_SCALE_FACTOR = 9.999999747378752e-05

    # (文件筛选和低云索引获取部分与之前相同)
    all_files = find_all_hdf_files(modis_data_root)
    file_list_filtered = [fp for fp in all_files if
                          start_year <= int(os.path.basename(fp).split('.')[1][1:5]) <= end_year]
    if not file_list_filtered:
        print(f"错误: 未找到 {start_year}-{end_year} 年期间的文件。")
        return
    print(f"找到 {len(file_list_filtered)} 个文件进行处理。")
    low_cloud_indices, success = get_low_cloud_indices(file_list_filtered[0], standard_pressure_intervals_var,
                                                       low_cloud_threshold)
    if not success: return

    regridded_clow_list = []
    print("\n--- 开始计算并重采样每月数据 ---")
    for file_path in file_list_filtered:
        try:
            # --- 计算各层云量 ---
            hdf_file = SD(file_path, SDC.READ)
            raw_total_cf = hdf_file.select(total_cloud_fraction_var).get().astype(np.float32)
            fill_val_cf = hdf_file.select(total_cloud_fraction_var).attributes().get('_FillValue')
            if fill_val_cf is not None: raw_total_cf[raw_total_cf == fill_val_cf] = np.nan
            total_cf = raw_total_cf * CF_SCALE_FACTOR
            p_histo_counts = hdf_file.select(pressure_histogram_var).get().astype(np.float32)
            fill_val_h = hdf_file.select(pressure_histogram_var).attributes().get('_FillValue')
            if fill_val_h is not None: p_histo_counts[p_histo_counts == fill_val_h] = np.nan
            total_histo_counts = np.nansum(p_histo_counts, axis=0)
            with np.errstate(divide='ignore', invalid='ignore'):
                p_histo_fractions = p_histo_counts / total_histo_counts
            p_histo_fractions = np.nan_to_num(p_histo_fractions, nan=0.0)
            cf_at_each_level = total_cf * p_histo_fractions
            hdf_file.end()

            # --- 使用随机重叠假设计算 C_low ---
            num_pressure_bins = cf_at_each_level.shape[0]
            valid_low_cloud_indices = [idx for idx in low_cloud_indices if idx < num_pressure_bins]

            # 1. 计算分子: 低云的总云量 (CF_CTP>680hPa)
            c_low_part = np.nansum(cf_at_each_level[valid_low_cloud_indices, :, :], axis=0)

            # 2. 计算分母中的中高云部分 (CF_CTP<680hPa)
            all_indices = np.arange(num_pressure_bins)
            high_mid_cloud_indices = np.setdiff1d(all_indices, valid_low_cloud_indices)
            c_high_mid_part = np.nansum(cf_at_each_level[high_mid_cloud_indices, :, :], axis=0)

            # 3. 应用随机重叠公式，并进行安全处理
            denominator = 1.0 - c_high_mid_part
            with np.errstate(divide='ignore', invalid='ignore'):
                c_low = c_low_part / denominator

            # 4. 将无效值(如除以0)和超出物理范围的值处理掉
            c_low[np.isinf(c_low)] = np.nan
            c_low = np.clip(c_low, 0, 1)

            # --- 重采样 ---
            lats_highres = np.linspace(89.5, -89.5, 180)
            lons_highres = np.linspace(-179.5, 179.5, 360)
            year_doy = os.path.basename(file_path).split('.')[1][1:]
            time = pd.to_datetime(year_doy, format='%Y%j')

            da_highres = xr.DataArray(c_low, coords={'latitude': lats_highres, 'longitude': lons_highres},
                                      dims=['latitude', 'longitude'])
            da_highres = da_highres.assign_coords(longitude=(((da_highres.longitude + 360) % 360))).sortby('longitude')
            da_highres_sorted = da_highres.sortby('latitude')
            da_regridded = da_highres_sorted.coarsen(latitude=2, longitude=5, boundary='trim').mean()
            da_regridded['time'] = time

            regridded_clow_list.append(da_regridded)
            print(f"  √ 已处理并重采样: {os.path.basename(file_path)}")

        except Exception as e:
            print(f"  × 处理文件时发生错误: {os.path.basename(file_path)} -> {e}")
            continue

    if regridded_clow_list:
        print("\n--- 正在合并所有月份数据并保存到最终文件 ---")
        ds_template = xr.open_dataset(calipso_template_file)
        final_da = xr.concat(regridded_clow_list, dim='time')
        final_da = final_da.rename('modis_lcc')

        final_ds = xr.Dataset(
            {'modis_lcc': (('time', 'latitude', 'longitude'), final_da.values)},
            coords=ds_template.coords
        )
        final_ds.attrs[
            'description'] = 'MODIS low cloud fraction (LCC) using random overlap assumption, regridded to CALIPSO 2x5 degree grid. Latitude sorted.'

        final_ds.to_netcdf(output_filename)
        print(f"\n🎉 修正后的文件已成功保存到: {output_filename}")
    else:
        print("\n没有数据被成功处理。")


if __name__ == '__main__':
    calculate_and_save_modis_clow_correct_order()
