# # -*- coding: utf-8 -*-
# from pyhdf.SD import SD, SDC
# import os
# import numpy as np
# from netCDF4 import Dataset
# from tqdm import tqdm  # 导入tqdm库以使用进度条
#
# # ================= 配置参数 =================
# # <--- 修改点 1: 更新数据根目录 ---
# base_dir = r"G:\\"  # 数据根目录
#
# lat_res = 1.0  # 1°纬度分辨率
# lon_res = 1.0  # 1°经度分辨率
#
# # <--- 修改点 2: 更新处理的年份范围 ---
# start_year = 2021
# end_year = 2021
#
#
# # ================= 路径验证与创建 =================
# def verify_and_create_dir(path):
#     """验证并创建目录（如果不存在）"""
#     try:
#         os.makedirs(path, exist_ok=True)
#         return True
#     except Exception as e:
#         print(f"❌ 无法创建目录 {path}: {str(e)}")
#         return False
#
#
# def verify_output_file(file_path):
#     """验证输出文件路径是否可写"""
#     dir_path = os.path.dirname(file_path)
#     return verify_and_create_dir(dir_path)
#
#
# # ================= 初始化网格 =================
# lat_bins = np.arange(-90, 90 + lat_res, lat_res)
# lon_bins = np.arange(-180, 180 + lon_res, lon_res)
# n_lat = len(lat_bins) - 1
# n_lon = len(lon_bins) - 1
# n_months = 12
# lat_centers = (lat_bins[:-1] + lat_bins[1:]) / 2
# lon_centers = (lon_bins[:-1] + lon_bins[1:]) / 2
#
#
# # ================= 核心处理函数 =================
# def process_file(file_path, lcc_sum_month, counts_month):
#     """处理单个HDF文件，更新当月的统计数组"""
#     try:
#         hdf = SD(file_path, SDC.READ)
#     except Exception:
#         return 0
#
#     try:
#         feature_flags = hdf.select("Feature_Classification_Flags")[:]
#         layer_top_pressure = hdf.select("Layer_Top_Pressure")[:]
#         opacity_flag = hdf.select("Opacity_Flag")[:]
#         num_layers = hdf.select("Number_Layers_Found")[:].flatten()
#         latitude = hdf.select("Latitude")[:].flatten()
#         longitude = hdf.select("Longitude")[:].flatten()
#         hdf.end()
#     except Exception:
#         hdf.end()
#         return 0
#
#     if len(latitude) == 0:
#         return 0
#
#     longitude = np.where(longitude > 180, longitude - 360, longitude)
#     longitude = np.where(longitude < -180, longitude + 360, longitude)
#
#     lat_bin = np.digitize(latitude, lat_bins) - 1
#     lon_bin = np.digitize(longitude, lon_bins) - 1
#
#     valid_mask = (lat_bin >= 0) & (lat_bin < n_lat) & (lon_bin >= 0) & (lon_bin < n_lon)
#     if not np.any(valid_mask):
#         return 0
#
#     valid_profiles = np.where(valid_mask)[0]
#
#     for i in valid_profiles:
#         layers = int(num_layers[i])
#         if layers <= 0:
#             lcc = 0
#         else:
#             top_pressures = layer_top_pressure[i, :layers]
#             flags = feature_flags[i, :layers]
#             opacity = opacity_flag[i, :layers]
#             is_cloud = (flags & 7) == 2
#             has_opaque_mid_high = np.any(is_cloud & (top_pressures <= 680) & (opacity == 1))
#             lcc = 0
#             if not has_opaque_mid_high and np.any(is_cloud & (top_pressures > 680)):
#                 lcc = 1
#         lcc_sum_month[lat_bin[i], lon_bin[i]] += lcc
#         counts_month[lat_bin[i], lon_bin[i]] += 1
#
#     return len(valid_profiles)
#
#
# # ================= 年份处理函数 =================
# def process_year(year):
#     """处理单个年份的所有月份数据"""
#     year_dir = os.path.join(base_dir, str(year))
#     if not os.path.exists(year_dir):
#         print(f"⚠️ 年份目录不存在：{year_dir}，跳过该年份")
#         return False
#
#     output_file = os.path.join(year_dir, f"calipso_lcc_monthly_1deg_{year}.nc")
#
#     if os.path.exists(output_file):
#         print(f"🔄 发现旧结果文件，将强制重新计算并覆盖: {output_file}")
#
#     lcc_sum_year = np.full((n_months, n_lat, n_lon), 0, dtype=np.float32)
#     counts_year = np.full((n_months, n_lat, n_lon), 0, dtype=np.int32)
#
#     # <--- 修改点 3: 修正临时文件的保存路径，使其与base_dir保持一致 ---
#     temp_save = os.path.join(year_dir, f"temp_lcc_stats_{year}.npz")
#     if os.path.exists(temp_save):
#         os.remove(temp_save)
#         print(f"🧹 已删除旧的中间结果文件，将从头开始计算: {temp_save}")
#
#     for month in range(1, 13):
#         month_idx = month - 1
#         month_dir = os.path.join(year_dir, f"{month:02d}")
#
#         if not os.path.exists(month_dir):
#             print(f"  - {year}年{month}月目录不存在，跳过。")
#             continue
#
#         print(f"\n📅 处理 {year}年{month}月")
#         file_list = [f for f in os.listdir(month_dir) if f.endswith(".hdf")]
#         if not file_list:
#             print("  📂 该月无HDF文件，跳过")
#             continue
#
#         total_profiles = 0
#         with tqdm(total=len(file_list), desc="处理文件", unit="个", leave=False) as pbar:
#             for filename in file_list:
#                 file_path = os.path.join(month_dir, filename)
#                 profiles = process_file(file_path, lcc_sum_year[month_idx], counts_year[month_idx])
#                 total_profiles += profiles
#                 pbar.update(1)
#
#         print(f"  ✅ 完成该月处理，有效剖面数：{total_profiles}")
#
#         try:
#             np.savez_compressed(temp_save, lcc_sum=lcc_sum_year, counts=counts_year)
#         except Exception as e:
#             print(f"  ⚠️ 中间结果保存失败：{e}")
#
#     print(f"\n📊 计算 {year}年逐月平均LCC...")
#     with np.errstate(divide='ignore', invalid='ignore'):
#         lcc_monthly_year = np.where(counts_year > 0, lcc_sum_year / counts_year, np.nan)
#
#     try:
#         with Dataset(output_file, "w", format="NETCDF4") as nc:
#             nc.createDimension("month", n_months)
#             nc.createDimension("lat", n_lat)
#             nc.createDimension("lon", n_lon)
#             year_var = nc.createVariable("year", "i4", ())
#             month_var = nc.createVariable("month", "i4", ("month",))
#             lat_var = nc.createVariable("lat", "f4", ("lat",))
#             lon_var = nc.createVariable("lon", "f4", ("lon",))
#             lcc_var = nc.createVariable("lcc", "f4", ("month", "lat", "lon"), fill_value=np.nan)
#             year_var[:] = year
#             month_var[:] = np.arange(1, 13)
#             lat_var[:] = lat_centers
#             lon_var[:] = lon_centers
#             lcc_var[:] = lcc_monthly_year
#             nc.description = f"CALIPSO Monthly Low Cloud Cover (LCC) on 1x1 degree grid for {year}"
#             lat_var.units = "degrees_north"
#             lon_var.units = "degrees_east"
#             lcc_var.units = "fraction (0-1)"
#             lcc_var.long_name = "Low Cloud Cover"
#         print(f"💾 该年份结果已保存：{output_file}")
#         if os.path.exists(temp_save):
#             os.remove(temp_save)
#             print(f"🧹 清理该年份中间文件：{temp_save}")
#         return True
#     except Exception as e:
#         print(f"❌ 保存该年份结果失败：{e}")
#         return False
#
#
# def main():
#     print(f"===== 开始处理 {start_year}-{end_year}年CALIPSO数据 (单线程稳定模式) =====")
#     for year in range(start_year, end_year + 1):
#         print(f"\n========== 处理年份：{year} ==========")
#         success = process_year(year)
#         status = "成功" if success else "失败"
#         print(f"========== {year}年处理{status} ==========")
#     print("\n🎉 所有年份处理完成！")
#
#
# if __name__ == "__main__":
#     main()
#


# -*- coding: utf-8 -*-
# FILENAME: calculate_calipso_monthly_structured_output.py
from pyhdf.SD import SD, SDC
import os
import numpy as np
from netCDF4 import Dataset
from tqdm import tqdm
import time
from pathlib import Path

# ===================================================================
# --- 1. 配置 (已按您的要求更新) ---
# ===================================================================
# CALIPSO HDF 数据所在的根目录
PROJECT_DIR = Path(__file__).resolve().parents[2]
base_dir = PROJECT_DIR / "raw_inputs" / "calipso"

# 要处理的年份范围 (设置为2014年)
start_year = 2007
end_year = 2021

# 输出文件的【主】目录
main_output_dir = PROJECT_DIR / "derived_intermediate" / "calipso_lcc_1deg"

# 网格分辨率
lat_res = 1.0
lon_res = 1.0

# ===================================================================
# --- 2. 初始化网格 (无需修改) ---
# ===================================================================
lat_bins = np.arange(-90, 90 + lat_res, lat_res)
lon_bins = np.arange(-180, 180 + lon_res, lon_res)
n_lat = len(lat_bins) - 1
n_lon = len(lon_bins) - 1
lat_centers = (lat_bins[:-1] + lat_bins[1:]) / 2
lon_centers = (lon_bins[:-1] + lon_bins[1:]) / 2


# ===================================================================
# --- 3. 核心处理函数 (无需修改) ---
# ===================================================================
def process_file(file_path, numerator_sum, denominator_unobscured_sum):
    """
    处理单个HDF文件，严格遵循 LCC = C(low) / (1 - C(opaque_high)) 的逻辑
    """
    try:
        hdf = SD(file_path, SDC.READ)
        feature_flags = hdf.select("Feature_Classification_Flags")[:]
        layer_top_pressure = hdf.select("Layer_Top_Pressure")[:]
        opacity_flag = hdf.select("Opacity_Flag")[:]
        num_layers = hdf.select("Number_Layers_Found")[:].flatten()
        latitude = hdf.select("Latitude")[:].flatten()
        longitude = hdf.select("Longitude")[:].flatten()
        hdf.end()
    except Exception:
        if 'hdf' in locals() and hdf: hdf.end()
        return

    if len(latitude) == 0: return

    longitude = np.where(longitude > 180, longitude - 360, longitude)

    lat_bin_indices = np.digitize(latitude, lat_bins) - 1
    lon_bin_indices = np.digitize(longitude, lon_bins) - 1

    valid_mask = (lat_bin_indices >= 0) & (lat_bin_indices < n_lat) & \
                 (lon_bin_indices >= 0) & (lon_bin_indices < n_lon)

    for i in np.where(valid_mask)[0]:
        lat_idx, lon_idx = lat_bin_indices[i], lon_bin_indices[i]
        layers = int(num_layers[i])

        has_opaque_mid_high = False

        if layers > 0:
            top_pressures = layer_top_pressure[i, :layers]
            flags = feature_flags[i, :layers]
            opacity = opacity_flag[i, :layers]
            is_cloud = (flags & 7) == 2
            has_opaque_mid_high = np.any(is_cloud & (top_pressures <= 680) & (opacity == 1))

        if not has_opaque_mid_high:
            denominator_unobscured_sum[lat_idx, lon_idx] += 1
            if layers > 0 and np.any(is_cloud & (top_pressures > 680)):
                numerator_sum[lat_idx, lon_idx] += 1


# ===================================================================
# --- 4. 按月份处理和保存的函数 (无需修改) ---
# ===================================================================
def process_and_save_month(year, month, month_dir, output_path):
    """处理并保存单个月份的数据"""
    if not os.path.exists(month_dir):
        return

    file_list = [f for f in os.listdir(month_dir) if f.endswith(".hdf")]
    if not file_list:
        return

    numerator_sum_month = np.zeros((n_lat, n_lon), dtype=np.float32)
    denominator_unobscured_month = np.zeros((n_lat, n_lon), dtype=np.int32)

    with tqdm(total=len(file_list), desc=f"处理 {year}-{month:02d}", unit="个", leave=False) as pbar:
        for filename in file_list:
            file_path = os.path.join(month_dir, filename)
            process_file(file_path, numerator_sum_month, denominator_unobscured_month)
            pbar.update(1)

    with np.errstate(divide='ignore', invalid='ignore'):
        lcc_month = np.where(denominator_unobscured_month > 0, numerator_sum_month / denominator_unobscured_month,
                             np.nan)

    try:
        with Dataset(output_path, "w", format="NETCDF4") as nc:
            nc.createDimension("latitude", n_lat)
            nc.createDimension("longitude", n_lon)
            lat_var = nc.createVariable("latitude", "f4", ("latitude",))
            lon_var = nc.createVariable("longitude", "f4", ("longitude",))
            lcc_var = nc.createVariable("lcc", "f4", ("latitude", "longitude"), fill_value=np.nan)

            lat_var[:] = lat_centers
            lon_var[:] = lon_centers
            lcc_var[:] = lcc_month

            nc.description = f"CALIPSO Conditional LCC for {year}-{month:02d}, using formula LCC = C(low) / (1 - C(opaque_high))."
            lat_var.units = "degrees_north"
            lon_var.units = "degrees_east"

        print(f"  -> ✅ 已保存到: {os.path.relpath(output_path, main_output_dir)}")
    except Exception as e:
        print(f"  -> ❌ 保存失败: {os.path.basename(output_path)} -> {e}")


# ===================================================================
# --- 5. 程序入口 (修改为处理2014年5月) ---
# ===================================================================
def main():
    print("开始处理CALIPSO LCC（按月/文件夹保存，用户指定公式）...")
    print(f"数据源: {base_dir}")
    print(f"主输出目录: {main_output_dir}")

    start_time = time.time()
    for year in range(start_year, end_year + 1):
        print(f"\n{'=' * 20} 处理年份: {year} {'=' * 20}")
        year_dir = os.path.join(base_dir, str(year))

        # --- 【修改点】设置为只处理5月 ---
        # 只处理5月
        start_month = 5
        end_month = 5
        # --------------------

        for month in range(start_month, end_month + 1):
            month_dir = os.path.join(year_dir, f"{month:02d}")

            output_month_dir = os.path.join(main_output_dir, str(year), f"{month:02d}")
            os.makedirs(output_month_dir, exist_ok=True)

            output_file = os.path.join(output_month_dir, f"CALIPSO_LCC_Conditional_{year}-{month:02d}.nc")

            process_and_save_month(year, month, month_dir, output_file)

    end_time = time.time()
    print(f"\n🎉 所有年份 ({start_year}-{end_year}) 处理完成！总耗时: {(end_time - start_time) / 60:.2f} 分钟。")


if __name__ == "__main__":
    main()
