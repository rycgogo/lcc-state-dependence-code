# import xarray as xr
# import numpy as np
# import os
#
# # 定义基础路径和文件名
# base_dir = r"E:/"
# nc_file = os.path.join(base_dir, "Unified_Data_Clean.nc")
# grib_file = os.path.join(base_dir, "10m wind.grib")
# output_file = os.path.join(base_dir, "Unified_Data_Clean_with_Tadv.nc")
#
# # 读取原有的干净 NC 数据和新的风场 GRIB 数据
# ds_nc = xr.open_dataset(nc_file)
# ds_grib = xr.open_dataset(grib_file, engine='cfgrib')
#
# # 提取经向和纬向风速变量
# # 注：在 cfgrib 引擎加载时，10u 和 10v 的短名通常会被标准化识别为 u10 和 v10
# u10 = ds_grib['u10']
# v10 = ds_grib['v10']
#
# # 将风场数据对齐到原有 NC 数据的经纬度和时间网格上
# # 这一步确保了后续风速和温度在同一个点上进行乘积运算
# u10_aligned = u10.interp_like(ds_nc)
# v10_aligned = v10.interp_like(ds_nc)
#
# # 设定地球平均半径常量 (单位：米)
# R_E = 6371000.0
#
# # 提取纬度数据并转换为弧度，用于计算纬度圈上的实际物理距离
# lat_rad = np.deg2rad(ds_nc['latitude'])
# cos_lat = np.cos(lat_rad)
#
# # 使用 xarray 的方法计算表面温度（ts）在经纬度上的中心差分
# # 此时求导结果的单位是 K/degree
# dTs_dlon_deg = ds_nc['ts'].differentiate('longitude')
# dTs_dlat_deg = ds_nc['ts'].differentiate('latitude')
#
# # 将偏导数中的度数转换为弧度，进而转化为实际物理距离的偏导数 (K/m)
# dTs_dlon_rad = dTs_dlon_deg * (180.0 / np.pi)
# dTs_dlat_rad = dTs_dlat_deg * (180.0 / np.pi)
#
# dT_dx = dTs_dlon_rad / (R_E * cos_lat)
# dT_dy = dTs_dlat_rad / R_E
#
# # 严格按照文献公式计算表面温度平流 Tadv
# tadv = - (u10_aligned * dT_dx + v10_aligned * dT_dy)
#
# # 将计算得出的平流结果合并回原始数据集，并增加必要的属性元数据
# ds_nc['tadv'] = tadv
# ds_nc['tadv'].attrs['long_name'] = 'Surface Temperature Advection'
# ds_nc['tadv'].attrs['units'] = 'K s**-1'
# ds_nc['tadv'].attrs['description'] = 'Calculated using 10m wind from ERA5 and surface temperature gradients'
#
# # 导出为一个带有全新变量的整合数据文件
# ds_nc.to_netcdf(output_file)
#
# # 释放内存并关闭文件读取
# ds_nc.close()
# ds_grib.close()
#
# print(f"数据处理已成功完成，包含 tadv 的新文件已存放至：{output_file}")

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs  # 引入 cartopy 用于绘制地图轮廓
from pathlib import Path

# Use the archived final analysis dataset. The original pre-Tadv merged input
# is not a separate publication artifact; its derived Tadv field is included
# in this NetCDF file.
project_dir = Path(__file__).resolve().parents[2]
file_path = project_dir / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc"

# 读取 NetCDF 文件
ds = xr.open_dataset(file_path)

# 提取 tadv 变量并计算时间平均
tadv_mean = ds['tadv'].mean(dim='time', keep_attrs=True)

# 单位转换：从 K/s 转换为 K/day (乘以 86400)
# 这样色轴范围 ±3 代表的是 ±3 K/day
tadv_day = tadv_mean * 86400

# 创建画布，指定投影
plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())

# 绘图
# xarray 的 plot 方法支持传入 ax 参数
tadv_day.plot(
    ax=ax,                # 将图画在带有投影的坐标轴上
    cmap='coolwarm',
    center=0,
    vmin=-3,              # 范围保持不变：-3
    vmax=3,               # 范围保持不变：3
    cbar_kwargs={'label': 'Surface Temperature Advection (K/day)'}
)

# 添加海岸线轮廓
# linewidth=0.8 设置线条粗细，color='black' 设置颜色
ax.coastlines(linewidth=0.8, color='black')

# 添加标题
plt.title('Global Distribution of Time-Mean Surface Temperature Advection (K/day)')

plt.show()

# 释放内存
ds.close()
