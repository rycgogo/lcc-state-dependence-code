# # # # import xarray as xr
# # # # import numpy as np
# # # #
# # # # # 定义输入和输出文件
# # # # input_file = r"E:\data2\eis\8d25c5beb882a1957d092c7c78d7a349\data_0.nc"
# # # # output_file = r"E:\data2\eis\lts_eis_2019-2021_1deg.nc"
# # # #
# # # # # 定义常量
# # # # Rd = 287.0      # 干空气气体常数, J/(kg·K)
# # # # cp = 1004.0     # 干空气定压比热, J/(kg·K)
# # # # g = 9.80665     # 重力加速度, m/s²
# # # # p0 = 1000.0     # 基准气压, hPa
# # # # Lv = 2.5e6      # 汽化潜热, J/kg
# # # #
# # # # # 定义计算 LCL 的函数
# # # # def calc_lcl_bolton(T, rh):
# # # #     """
# # # #     计算 LCL 高度 (m)
# # # #     T  : 温度 (K)
# # # #     rh : 相对湿度 (0~1)
# # # #     使用 Bolton (1980) 的近似公式：
# # # #         LCL ~ 125 * (T_c - Td_c)
# # # #     其中 T_c 为摄氏温度，Td_c 为露点温度 (°C)
# # # #     """
# # # #     T_c = T - 273.15
# # # #     e_s = 6.112 * np.exp(17.67 * T_c / (T_c + 243.5))
# # # #     e = rh * e_s
# # # #     Td_c = 243.5 * np.log(e/6.112) / (17.67 - np.log(e/6.112))
# # # #     lcl = 125.0 * (T_c - Td_c)
# # # #     lcl = xr.where(lcl < 0, 0, lcl)
# # # #     return lcl
# # # #
# # # # # 定义计算饱和比湿的函数
# # # # def saturation_specific_humidity(T, p_hpa):
# # # #     """
# # # #     计算饱和比湿 q_s (kg/kg)
# # # #     T    : 温度 (K)
# # # #     p_hpa: 压力 (hPa)
# # # #     """
# # # #     T_c = T - 273.15
# # # #     e_s_hpa = 6.112 * np.exp(17.67 * T_c / (T_c + 243.5))
# # # #     e_s_pa = e_s_hpa * 100.0
# # # #     p_pa = p_hpa * 100.0
# # # #     epsilon = 0.622
# # # #     q_s = epsilon * e_s_pa / (p_pa - e_s_pa)
# # # #     return q_s
# # # #
# # # # # 打开数据集（利用 dask 分块降低内存占用）
# # # # ds = xr.open_dataset(input_file, chunks={"valid_time": 10, "latitude": 100, "longitude": 100})
# # # #
# # # # # 如果存在 valid_time，将其重命名为 time
# # # # if 'valid_time' in ds.coords or 'valid_time' in ds.dims:
# # # #     ds = ds.rename({'valid_time': 'time'})
# # # #
# # # # # 删除 expver 和 number 维度或坐标（如果存在）
# # # # if 'expver' in ds.dims:
# # # #     ds = ds.isel(expver=0)
# # # # if 'expver' in ds.coords:
# # # #     ds = ds.drop('expver')
# # # # if 'number' in ds.dims:
# # # #     ds = ds.isel(number=0)
# # # # if 'number' in ds.coords:
# # # #     ds = ds.drop('number')
# # # #
# # # # # 计算36个月的平均值
# # # # ds_mean = ds.mean(dim='time')
# # # #
# # # # # 重新网格化为 1°×1°
# # # # # 定义新的 1°×1° 网格
# # # # lat_bins = np.arange(-90, 91, 1.0)  # -90 to 90, 1° resolution
# # # # lon_bins = np.arange(0, 360, 1.0)   # 0 to 360, 1° resolution
# # # #
# # # # # 使用 xarray 的 groupby_bins 进行重网格化（平均）
# # # # # 先对 latitude 进行分组
# # # # ds_regrid = ds_mean.groupby_bins('latitude', bins=lat_bins, labels=lat_bins[:-1]+0.5).mean()
# # # # # 再对 longitude 进行分组
# # # # ds_regrid = ds_regrid.groupby_bins('longitude', bins=lon_bins, labels=lon_bins[:-1]+0.5).mean()
# # # #
# # # # # 重命名 bins 后的维度
# # # # ds_regrid = ds_regrid.rename({'latitude_bins': 'latitude', 'longitude_bins': 'longitude'})
# # # #
# # # # # 确保维度顺序正确
# # # # ds_regrid = ds_regrid.transpose('pressure_level', 'latitude', 'longitude')
# # # #
# # # # # 分别选取 700 hPa 和 1000 hPa 层的数据
# # # # ds_700 = ds_regrid.sel(pressure_level=700)
# # # # ds_1000 = ds_regrid.sel(pressure_level=1000)
# # # #
# # # # # 计算位温（Potential Temperature）
# # # # theta_700 = ds_700['t'] * (p0 / 700.0)**(Rd/cp)
# # # # theta_1000 = ds_1000['t'] * (p0 / 1000.0)**(Rd/cp)
# # # #
# # # # # LTS 定义为：θ_700 - θ_1000
# # # # LTS = theta_700 - theta_1000
# # # #
# # # # # 计算几何高度（由 geopotential 转换，单位：m）
# # # # z_700 = ds_700['z'] / g
# # # # z_1000 = ds_1000['z'] / g
# # # #
# # # # # 计算 LCL（抬升凝结高度）
# # # # rh_1000 = ds_1000['r'] / 100.0
# # # # z_LCL = calc_lcl_bolton(ds_1000['t'], rh_1000)
# # # #
# # # # # 计算湿绝热递减率 Gamma_m（近似取850 hPa）
# # # # p_850 = 850.0
# # # # T_850 = 0.5 * (ds_700['t'] + ds_1000['t'])
# # # # q_s_850 = saturation_specific_humidity(T_850, p_850)
# # # #
# # # # Gamma_m = (1000.0 / p_850) * (g / cp) * (
# # # #     1.0 - (1.0 + (Lv * q_s_850) / (Rd * T_850)) /
# # # #     (1.0 + (Lv**2 * q_s_850) / (cp * Rd * T_850**2))
# # # # )
# # # #
# # # # # 计算 EIS
# # # # EIS = LTS - Gamma_m * (z_700 - z_LCL)
# # # #
# # # # # 确保 LTS 和 EIS 的维度正确
# # # # LTS = LTS.transpose('latitude', 'longitude')
# # # # EIS = EIS.transpose('latitude', 'longitude')
# # # #
# # # # # 打印维度以调试
# # # # print("LTS dimensions:", LTS.dims, LTS.shape)
# # # # print("EIS dimensions:", EIS.dims, EIS.shape)
# # # # print("Latitude length:", len(ds_regrid['latitude']))
# # # # print("Longitude length:", len(ds_regrid['longitude']))
# # # #
# # # # # 创建输出数据集
# # # # ds_out = xr.Dataset(
# # # #     {
# # # #         'LTS': (['latitude', 'longitude'], LTS.data),
# # # #         'EIS': (['latitude', 'longitude'], EIS.data)
# # # #     },
# # # #     coords={
# # # #         'latitude': ds_regrid['latitude'],
# # # #         'longitude': ds_regrid['longitude']
# # # #     }
# # # # )
# # # #
# # # # # 保存结果到 netCDF 文件
# # # # ds_out.to_netcdf(output_file)
# # # # print(f"处理完成，结果已保存到 {output_file}。")
# # # #
# # # # # 计算并输出统计信息
# # # # def calculate_stats(data_array):
# # # #     stats = {
# # # #         'min': float(data_array.min().compute()),
# # # #         'max': float(data_array.max().compute()),
# # # #         'mean': float(data_array.mean().compute()),
# # # #         'std': float(data_array.std().compute()),
# # # #         '25%': float(np.nanpercentile(data_array.compute(), 25)),
# # # #         '50%': float(np.nanpercentile(data_array.compute(), 50)),
# # # #         '75%': float(np.nanpercentile(data_array.compute(), 75))
# # # #     }
# # # #     return stats
# # # #
# # # # # 计算 EIS 的统计信息
# # # # eis_stats = calculate_stats(ds_out['EIS'])
# # # # print("EIS 完整统计信息：")
# # # # for stat, value in eis_stats.items():
# # # #     print(f"{stat}: {value}")
# # # #
# # # # # 计算 LTS 的统计信息
# # # # lts_stats = calculate_stats(ds_out['LTS'])
# # # # print("\nLTS 完整统计信息：")
# # # # for stat, value in lts_stats.items():
# # # #     print(f"{stat}: {value}")
# # # #
# # # # # 关闭数据集
# # # # ds.close()
# # # # ds_out.close()
# # #
# # # import xarray as xr
# # # import numpy as np
# # #
# # # # 定义输入和输出文件
# # # input_file = r"E:\data2\eis\data_0.nc"
# # # output_file_monthly = r"E:\data2\eis\lts_eis_2019-2021_1deg_monthly.nc"
# # # # output_file = r"E:\data2\eis\lts_eis_2019-2021_1deg.nc"  # 原输出文件注释
# # #
# # # # 定义常量
# # # Rd = 287.0      # 干空气气体常数, J/(kg·K)
# # # cp = 1004.0     # 干空气定压比热, J/(kg·K)
# # # g = 9.80665     # 重力加速度, m/s²
# # # p0 = 1000.0     # 基准气压, hPa
# # # Lv = 2.5e6      # 汽化潜热, J/kg
# # #
# # # # 定义计算 LCL 的函数
# # # def calc_lcl_bolton(T, rh):
# # #     """
# # #     计算 LCL 高度 (m)
# # #     T  : 温度 (K)
# # #     rh : 相对湿度 (0~1)
# # #     使用 Bolton (1980) 的近似公式：
# # #         LCL ~ 125 * (T_c - Td_c)
# # #     其中 T_c 为摄氏温度，Td_c 为露点温度 (°C)
# # #     """
# # #     T_c = T - 273.15
# # #     e_s = 6.112 * np.exp(17.67 * T_c / (T_c + 243.5))
# # #     e = rh * e_s
# # #     Td_c = 243.5 * np.log(e/6.112) / (17.67 - np.log(e/6.112))
# # #     lcl = 125.0 * (T_c - Td_c)
# # #     lcl = xr.where(lcl < 0, 0, lcl)
# # #     return lcl
# # #
# # # # 定义计算饱和比湿的函数
# # # def saturation_specific_humidity(T, p_hpa):
# # #     """
# # #     计算饱和比湿 q_s (kg/kg)
# # #     T    : 温度 (K)
# # #     p_hpa: 压力 (hPa)
# # #     """
# # #     T_c = T - 273.15
# # #     e_s_hpa = 6.112 * np.exp(17.67 * T_c / (T_c + 243.5))
# # #     e_s_pa = e_s_hpa * 100.0
# # #     p_pa = p_hpa * 100.0
# # #     epsilon = 0.622
# # #     q_s = epsilon * e_s_pa / (p_pa - e_s_pa)
# # #     return q_s
# # #
# # # # 打开数据集（利用 dask 分块降低内存占用）
# # # ds = xr.open_dataset(input_file, chunks={"valid_time": 10, "latitude": 100, "longitude": 100})
# # #
# # # # 如果存在 valid_time，将其重命名为 time
# # # if 'valid_time' in ds.coords or 'valid_time' in ds.dims:
# # #     ds = ds.rename({'valid_time': 'time'})
# # #
# # # # 删除 expver 和 number 维度或坐标（如果存在）
# # # if 'expver' in ds.dims:
# # #     ds = ds.isel(expver=0)
# # # if 'expver' in ds.coords:
# # #     ds = ds.drop('expver')
# # # if 'number' in ds.dims:
# # #     ds = ds.isel(number=0)
# # # if 'number' in ds.coords:
# # #     ds = ds.drop('number')
# # #
# # # # 定义新的 1°×1° 网格
# # # lat_bins = np.arange(-90, 91, 1.0)  # -90 to 90, 1° resolution
# # # lon_bins = np.arange(0, 360, 1.0)   # 0 to 360, 1° resolution
# # #
# # # # 按月重新网格化（保留时间维度）
# # # # 使用 groupby_bins 进行重网格化（平均）
# # # ds_regrid = ds.groupby_bins('latitude', bins=lat_bins, labels=lat_bins[:-1]+0.5).mean()
# # # ds_regrid = ds_regrid.groupby_bins('longitude', bins=lon_bins, labels=lon_bins[:-1]+0.5).mean()
# # #
# # # # 重命名 bins 后的维度
# # # ds_regrid = ds_regrid.rename({'latitude_bins': 'latitude', 'longitude_bins': 'longitude'})
# # #
# # # # 确保维度顺序正确
# # # ds_regrid = ds_regrid.transpose('time', 'pressure_level', 'latitude', 'longitude')
# # #
# # # # 分别选取 700 hPa 和 1000 hPa 层的数据
# # # ds_700 = ds_regrid.sel(pressure_level=700)
# # # ds_1000 = ds_regrid.sel(pressure_level=1000)
# # #
# # # # 计算位温（Potential Temperature）
# # # theta_700 = ds_700['t'] * (p0 / 700.0)**(Rd/cp)
# # # theta_1000 = ds_1000['t'] * (p0 / 1000.0)**(Rd/cp)
# # #
# # # # LTS 定义为：θ_700 - θ_1000
# # # LTS = theta_700 - theta_1000
# # #
# # # # 计算几何高度（由 geopotential 转换，单位：m）
# # # z_700 = ds_700['z'] / g
# # # z_1000 = ds_1000['z'] / g
# # #
# # # # 计算 LCL（抬升凝结高度）
# # # rh_1000 = ds_1000['r'] / 100.0
# # # z_LCL = calc_lcl_bolton(ds_1000['t'], rh_1000)
# # #
# # # # 计算湿绝热递减率 Gamma_m（近似取850 hPa）
# # # p_850 = 850.0
# # # T_850 = 0.5 * (ds_700['t'] + ds_1000['t'])
# # # q_s_850 = saturation_specific_humidity(T_850, p_850)
# # #
# # # Gamma_m = (1000.0 / p_850) * (g / cp) * (
# # #     1.0 - (1.0 + (Lv * q_s_850) / (Rd * T_850)) /
# # #     (1.0 + (Lv**2 * q_s_850) / (cp * Rd * T_850**2))
# # # )
# # #
# # # # 计算 EIS
# # # EIS = LTS - Gamma_m * (z_700 - z_LCL)
# # #
# # # # 确保 EIS 的维度正确
# # # EIS = EIS.transpose('time', 'latitude', 'longitude')
# # #
# # # # 打印维度以调试
# # # print("EIS dimensions:", EIS.dims, EIS.shape)
# # # print("Latitude length:", len(ds_regrid['latitude']))
# # # print("Longitude length:", len(ds_regrid['longitude']))
# # #
# # # # 创建输出数据集（仅保存 EIS 和时间、经纬度）
# # # ds_out_monthly = xr.Dataset(
# # #     {
# # #         'EIS': (['time', 'latitude', 'longitude'], EIS.data)
# # #     },
# # #     coords={
# # #         'time': ds_regrid['time'],
# # #         'latitude': ds_regrid['latitude'],
# # #         'longitude': ds_regrid['longitude']
# # #     }
# # # )
# # #
# # # # 保存月度结果到 netCDF 文件
# # # ds_out_monthly.to_netcdf(output_file_monthly)
# # # print(f"处理完成，月度结果已保存到 {output_file_monthly}。")
# # #
# # # # # 以下为原代码的注释部分
# # # # """
# # # # # 计算36个月的平均值
# # # # ds_mean = ds.mean(dim='time')
# # # #
# # # # # 重新网格化为 1°×1°
# # # # ds_regrid = ds_mean.groupby_bins('latitude', bins=lat_bins, labels=lat_bins[:-1]+0.5).mean()
# # # # ds_regrid = ds_regrid.groupby_bins('longitude', bins=lon_bins, labels=lon_bins[:-1]+0.5).mean()
# # # # ds_regrid = ds_regrid.rename({'latitude_bins': 'latitude', 'longitude_bins': 'longitude'})
# # # # ds_regrid = ds_regrid.transpose('pressure_level', 'latitude', 'longitude')
# # # #
# # # # ds_700 = ds_regrid.sel(pressure_level=700)
# # # # ds_1000 = ds_regrid.sel(pressure_level=1000)
# # # #
# # # # theta_700 = ds_700['t'] * (p0 / 700.0)**(Rd/cp)
# # # # theta_1000 = ds_1000['t'] * (p0 / 1000.0)**(Rd/cp)
# # # # LTS = theta_700 - theta_1000
# # # #
# # # # z_700 = ds_700['z'] / g
# # # # z_1000 = ds_1000['z'] / g
# # # # rh_1000 = ds_1000['r'] / 100.0
# # # # z_LCL = calc_lcl_bolton(ds_1000['t'], rh_1000)
# # # #
# # # # p_850 = 850.0
# # # # T_850 = 0.5 * (ds_700['t'] + ds_1000['t'])
# # # # q_s_850 = saturation_specific_humidity(T_850, p_850)
# # # # Gamma_m = (1000.0 / p_850) * (g / cp) * (
# # # #     1.0 - (1.0 + (Lv * q_s_850) / (Rd * T_850)) /
# # # #     (1.0 + (Lv**2 * q_s_850) / (cp * Rd * T_850**2))
# # # # )
# # # #
# # # # EIS = LTS - Gamma_m * (z_700 - z_LCL)
# # # # LTS = LTS.transpose('latitude', 'longitude')
# # # # EIS = EIS.transpose('latitude', 'longitude')
# # # #
# # # # print("LTS dimensions:", LTS.dims, LTS.shape)
# # # # print("EIS dimensions:", EIS.dims, EIS.shape)
# # # # print("Latitude length:", len(ds_regrid['latitude']))
# # # # print("Longitude length:", len(ds_regrid['longitude']))
# # # #
# # # # ds_out = xr.Dataset(
# # # #     {
# # # #         'LTS': (['latitude', 'longitude'], LTS.data),
# # # #         'EIS': (['latitude', 'longitude'], EIS.data)
# # # #     },
# # # #     coords={
# # # #         'latitude': ds_regrid['latitude'],
# # # #         'longitude': ds_regrid['longitude']
# # # #     }
# # # # )
# # # #
# # # # ds_out.to_netcdf(output_file)
# # # # print(f"处理完成，结果已保存到 {output_file}。")
# # # #
# # # # def calculate_stats(data_array):
# # # #     stats = {
# # # #         'min': float(data_array.min().compute()),
# # # #         'max': float(data_array.max().compute()),
# # # #         'mean': float(data_array.mean().compute()),
# # # #         'std': float(data_array.std().compute()),
# # # #         '25%': float(np.nanpercentile(data_array.compute(), 25)),
# # # #         '50%': float(np.nanpercentile(data_array.compute(), 50)),
# # # #         '75%': float(np.nanpercentile(data_array.compute(), 75))
# # # #     }
# # # #     return stats
# # # #
# # # # eis_stats = calculate_stats(ds_out['EIS'])
# # # # print("EIS 完整统计信息：")
# # # # for stat, value in eis_stats.items():
# # # #     print(f"{stat}: {value}")
# # # #
# # # # lts_stats = calculate_stats(ds_out['LTS'])
# # # # print("\nLTS 完整统计信息：")
# # # # for stat, value in lts_stats.items():
# # # #     print(f"{stat}: {value}")
# # # #
# # # # ds_out.close()
# # # # """
# # #
# # # # 关闭数据集
# # # ds.close()
# # # ds_out_monthly.close()
# #
# #
# # import xarray as xr
# # import numpy as np
# # import pandas as pd
# # from datetime import datetime
# # import os
# # import warnings  # For handling potential warnings
# #
# # # Suppress harmless warnings, e.g., about NaN comparisons
# # warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered')
# #
# # # --- 物理常数 (SI units) ---
# # g = 9.80665  # Standard acceleration of gravity (m/s^2)
# # R_d = 287.058  # Gas constant for dry air (J/(kg*K))
# # R_v = 461.5  # Gas constant for water vapor (J/(kg*K))
# # c_pa = 1004.0  # Specific heat at constant pressure for dry air (J/(kg*K))
# # epsilon = R_d / R_v  # Ratio of gas constants (approx 0.622)
# # kappa = R_d / c_pa  # (approx 0.286)
# #
# # # Lv (Latent heat of vaporization) - using a common constant value
# # # Note: For highly precise calculations, Lv is temperature-dependent.
# # # However, for monthly mean atmospheric profiles and this formula, a constant is often used.
# # L_v = 2.501e6  # Latent heat of vaporization (J/kg) at 0°C (approx)
# #
# # P0_Pa = 100000.0  # Reference pressure for potential temperature (Pa), 1000 hPa
# #
# #
# # # --- 辅助函数：根据提供的公式和物理原理 ---
# #
# # def calculate_potential_temperature(T_K: xr.DataArray, P_hPa: [xr.DataArray, float]) -> xr.DataArray:
# #     """
# #     计算位温 (Potential Temperature).
# #     Formula: theta = T * (P0 / P)^(R_d / c_p)
# #     Args:
# #         T_K (xarray.DataArray): Temperature in Kelvin.
# #         P_hPa (float or xarray.DataArray): Pressure in hPa.
# #     Returns:
# #         xarray.DataArray: Potential temperature in Kelvin.
# #     """
# #     P_Pa = P_hPa * 100  # Convert hPa to Pa
# #     return T_K * (P0_Pa / P_Pa) ** kappa
# #
# #
# # def calculate_saturation_vapor_pressure(T_K: xr.DataArray) -> xr.DataArray:
# #     """
# #     计算饱和水汽压 (Saturation Vapor Pressure) 采用 Magnus-Tetens 公式.
# #     Args:
# #         T_K (xarray.DataArray): Temperature in Kelvin.
# #     Returns:
# #         xarray.DataArray: Saturation vapor pressure in Pa.
# #     """
# #     T_C = T_K - 273.15  # Convert T from K to Celsius
# #     # Coefficients for Magnus-Tetens formula over water (Bolton, 1980)
# #     es = 611.2 * np.exp((17.67 * T_C) / (T_C + 243.5))  # Pa
# #     return es
# #
# #
# # def calculate_saturation_specific_humidity(es_Pa: xr.DataArray, P_Pa: xr.DataArray) -> xr.DataArray:
# #     """
# #     计算饱和比湿 (Saturation Specific Humidity).
# #     Formula: qs = epsilon * es / (P - (1 - epsilon) * es)
# #     Args:
# #         es_Pa (xarray.DataArray): Saturation vapor pressure in Pa.
# #         P_Pa (xarray.DataArray): Pressure in Pa.
# #     Returns:
# #         xarray.DataArray: Saturation specific humidity (kg/kg).
# #     """
# #     # Use .broadcast_like if dimensions differ but are compatible
# #     # In xarray, arithmetic operations handle broadcasting automatically if dimensions align or are compatible.
# #
# #     denominator = P_Pa - (1 - epsilon) * es_Pa
# #
# #     # Handle cases where denominator might be zero or negative
# #     # Set qs to 0 where denominator is non-positive or es_Pa is NaN/very small.
# #     # Using np.fmax to ensure denominator is at least a very small positive number to avoid div by zero,
# #     # or handle explicitly with xr.where.
# #
# #     # More robust: use xr.where to set to NaN or 0 where denominator is problematic
# #     qs = xr.where(denominator > 0, epsilon * es_Pa / denominator, np.nan)
# #     return qs
# #
# #
# # def calculate_LCL_height(T_K: xr.DataArray, RH_percent: xr.DataArray, P_hPa: [xr.DataArray, float]) -> xr.DataArray:
# #     """
# #     计算举升凝结高度 (Lifted Condensation Level, LCL) 的高度 z_LCL.
# #     使用 Bolton (1980) 经验公式计算露点温度，然后用简化公式计算LCL高度。
# #     z_LCL = 125 * (T_dry_bulb_C - T_dew_point_C) (approx. in meters)
# #     Args:
# #         T_K (xarray.DataArray): Temperature in Kelvin.
# #         RH_percent (xarray.DataArray): Relative humidity in percentage (0-100).
# #         P_hPa (float or xarray.DataArray): Pressure in hPa.
# #     Returns:
# #         xarray.DataArray: LCL height in meters.
# #     """
# #     T_C = T_K - 273.15  # Temperature in Celsius
# #
# #     # Calculate actual vapor pressure (e) from RH and saturation vapor pressure (es)
# #     # Using Bolton's coefficients for es_hPa (from es_Pa / 100)
# #     e_s_hPa = calculate_saturation_vapor_pressure(T_K) / 100  # Convert Pa to hPa
# #     e_actual_hPa = e_s_hPa * (RH_percent / 100.0)  # Actual vapor pressure in hPa
# #
# #     # Calculate dew point temperature (Td_C) in Celsius (Bolton, 1980)
# #     # Handle cases where e_actual_hPa might be non-positive or very small, which would cause log(0) or log(negative)
# #     # By setting Td_C to NaN if e_actual_hPa is problematic.
# #     Td_C = xr.where(e_actual_hPa > 0,
# #                     (243.5 * np.log(e_actual_hPa / 6.112)) / (17.67 - np.log(e_actual_hPa / 6.112)),
# #                     np.nan)
# #
# #     # Simplified formula for LCL height (Stull, 1988), approximate in meters
# #     # This formula assumes height above the initial level. Given the context, this is usually acceptable.
# #     z_LCL = 125.0 * (T_C - Td_C)  # meters
# #
# #     # Ensure z_LCL is not negative (LCL cannot be below the surface)
# #     z_LCL = z_LCL.where(z_LCL >= 0, 0)  # Set negative LCL heights to 0 (at or near surface)
# #
# #     return z_LCL
# #
# #
# # def calculate_moist_adiabatic_gradient(T_K: xr.DataArray, P_hPa: [xr.DataArray, float], Lv: float, cp_a: float,
# #                                        Rv: float) -> xr.DataArray:
# #     """
# #     计算湿绝热位温梯度 Gamma_m.
# #     Formula: Gamma_m(T,p) = (g / cp_a) * (1 + Lv * qs(T,p) / (Rv * T)) / (1 + Lv^2 * qs(T,p) / (cp_a * Rv * T^2))
# #     Args:
# #         T_K (xarray.DataArray): Temperature in Kelvin.
# #         P_hPa (xarray.DataArray or float): Pressure in hPa.
# #         Lv (float): Latent heat of vaporization (J/kg).
# #         cp_a (float): Specific heat at constant pressure for dry air (J/(kg*K)).
# #         Rv (float): Gas constant for water vapor (J/(kg*K)).
# #     Returns:
# #         xarray.DataArray: Moist adiabatic potential temperature gradient (K/m).
# #     """
# #     P_Pa = P_hPa * 100  # Convert hPa to Pa
# #     es_Pa = calculate_saturation_vapor_pressure(T_K)
# #     qs = calculate_saturation_specific_humidity(es_Pa, P_Pa)  # kg/kg
# #
# #     # Avoid division by zero when qs is zero, or T_K is zero (unlikely for atmospheric temp)
# #     # Using .fillna(0) for qs might be too aggressive if NaN means missing data rather than zero moisture.
# #     # We proceed with NaNs to propagate them if inputs are NaN.
# #
# #     numerator_term = 1 + (Lv * qs) / (Rv * T_K)
# #     denominator_term = 1 + (Lv ** 2 * qs) / (cp_a * Rv * T_K ** 2)
# #
# #     # Handle division by zero or problematic denominator
# #     gamma_m = xr.where(denominator_term != 0, (g / cp_a) * (numerator_term / denominator_term), np.nan)
# #     return gamma_m
# #
# #
# # # --- 配置参数 ---
# # input_era5_file_path = r"F:\data\ERA5\7184c7e96123f1bb2da0d9d5af1a900.nc"
# # output_eis_file_path = r"F:\data\ERA5\ERA5_EIS_Omega_Monthly_1deg_2006-2023.nc"
# #
# # target_lat_res = 1.0
# # target_lon_res = 1.0
# #
# # start_year = 2006
# # end_year = 2023
# #
# # print(f"Starting EIS and Omega calculation from: {input_era5_file_path}")
# #
# # try:
# #     ds = xr.open_dataset(input_era5_file_path)
# #     print("ERA5 NetCDF file loaded successfully. Original dataset info:")
# #     # print(ds) # Uncomment for detailed debug info
# #
# #     # --- Pre-processing: Handle number/expver dimensions and longitude conversion ---
# #     # Squeeze or average 'number' dimension
# #     if 'number' in ds.dims:
# #         if ds.dims['number'] > 1:
# #             print("Warning: 'number' dimension detected with multiple members. Averaging across 'number'.")
# #             ds = ds.mean(dim='number', keep_attrs=True)  # keep_attrs to retain metadata
# #         else:
# #             ds = ds.squeeze('number', drop=True)
# #     elif 'number' in ds.coords:  # If 'number' is a coordinate but not a dimension
# #         ds = ds.drop_vars('number')
# #
# #     # Squeeze or average 'expver' dimension
# #     if 'expver' in ds.dims:
# #         if ds.dims['expver'] > 1:
# #             print("Warning: 'expver' dimension detected with multiple members. Averaging across 'expver'.")
# #             ds = ds.mean(dim='expver', keep_attrs=True)  # keep_attrs to retain metadata
# #         else:
# #             ds = ds.squeeze('expver', drop=True)
# #     elif 'expver' in ds.coords:  # If 'expver' is a coordinate but not a dimension
# #         ds = ds.drop_vars('expver')
# #
# #     # Convert longitudes from 0-360 to -180-180 range
# #     lon_min_orig = ds['longitude'].min().item()
# #     lon_max_orig = ds['longitude'].max().item()
# #
# #     if lon_min_orig >= 0 and lon_max_orig <= 360 and lon_max_orig > 180:
# #         print("Converting dataset longitudes from 0-360 to -180-180 range.")
# #         ds = ds.assign_coords(longitude=(((ds['longitude'] + 180) % 360) - 180))
# #         ds = ds.sortby('longitude')  # Must sort after conversion
# #     else:
# #         print("Longitudes appear to be already in -180 to 180 range or similar. No conversion needed.")
# #
# #     # Invert latitude if it's decreasing for consistency in interpolation
# #     if ds['latitude'].values[0] > ds['latitude'].values[-1]:
# #         print("Inverting latitude dimension to be increasing.")
# #         ds = ds.isel(latitude=slice(None, None, -1))  # Slice with step -1 to reverse order
# #         # Update attributes to reflect new order if necessary, or xarray handles it.
# #
# #     print(
# #         f"Dataset Longitude Range AFTER conversion: Min={ds['longitude'].min().item()}, Max={ds['longitude'].max().item()}")
# #     print(
# #         f"Dataset Latitude Range AFTER inversion: Min={ds['latitude'].min().item()}, Max={ds['latitude'].max().item()}")
# #
# #     # --- Extract variables at specific pressure levels ---
# #     # Using .sel with method='nearest' to handle exact level matching
# #     P_1000_hPa = 1000.0
# #     P_700_hPa = 700.0
# #     P_500_hPa = 500.0
# #
# #     # Ensure we select by the correct pressure level dimension name 'pressure_level'
# #     t_1000 = ds['t'].sel(pressure_level=P_1000_hPa, method='nearest')
# #     r_1000 = ds['r'].sel(pressure_level=P_1000_hPa, method='nearest')
# #     z_1000 = ds['z'].sel(pressure_level=P_1000_hPa, method='nearest')
# #
# #     t_700 = ds['t'].sel(pressure_level=P_700_hPa, method='nearest')
# #     z_700 = ds['z'].sel(pressure_level=P_700_hPa, method='nearest')
# #
# #     w_500 = ds['w'].sel(pressure_level=P_500_hPa, method='nearest')
# #
# #     # Convert geopotential (z, m^2/s^2) to geometric height (m)
# #     height_1000 = z_1000 / g
# #     height_700 = z_700 / g
# #
# #     # --- Calculate Potential Temperature ---
# #     theta_1000 = calculate_potential_temperature(t_1000, P_1000_hPa)
# #     theta_700 = calculate_potential_temperature(t_700, P_700_hPa)
# #
# #     # --- Calculate LCL Height (z_LCL) at 1000hPa ---
# #     # z_LCL uses T_K, RH_percent, P_hPa from the lowest level (1000hPa)
# #     z_LCL = calculate_LCL_height(t_1000, r_1000, P_1000_hPa)
# #
# #     # --- Calculate Moist Adiabatic Gradient (Gamma_m) ---
# #     # As per image: Gamma_m is at 850hPa using mean temperature of 1000 and 700 hPa levels.
# #     # We approximate the mean temperature of the layer.
# #     T_mean_layer = (t_1000 + t_700) / 2.0
# #     # The pressure for Gamma_m calculation will be 850hPa as specified in the formula image.
# #     P_gamma_m_hPa = 850.0
# #
# #     gamma_m = calculate_moist_adiabatic_gradient(T_mean_layer, P_gamma_m_hPa, L_v, c_pa, R_v)  # K/m
# #
# #     # --- Calculate LTS and EIS ---
# #     LTS = theta_700 - theta_1000  # K
# #
# #     # Ensure all terms have compatible dimensions for subtraction
# #     # EIS = LTS - Gamma_m * (z_700hPa - z_LCL)
# #     # Here, z_700hPa and z_LCL are in meters. Gamma_m is K/m. Result will be K.
# #     EIS = LTS - gamma_m * (height_700 - z_LCL)  # K
# #
# #     # --- Prepare output DataArray structure ---
# #     years_unique = np.arange(start_year, end_year + 1)
# #     months_unique = np.arange(1, 13)
# #
# #     target_latitudes = np.arange(-89.5, 90.5, target_lat_res)
# #     target_longitudes = np.arange(-179.5, 180.5, target_lon_res)
# #
# #     n_years = len(years_unique)
# #     n_months = len(months_unique)
# #     n_target_lat = len(target_latitudes)
# #     n_target_lon = len(target_longitudes)
# #
# #     eis_monthly_1deg = xr.DataArray(
# #         np.full((n_years, n_months, n_target_lat, n_target_lon), np.nan, dtype=np.float32),
# #         coords={
# #             'year': years_unique,
# #             'month': months_unique,
# #             'latitude': target_latitudes,
# #             'longitude': target_longitudes
# #         },
# #         dims=['year', 'month', 'latitude', 'longitude'],
# #         name='eis'
# #     )
# #
# #     omega_500hPa_monthly_1deg = xr.DataArray(
# #         np.full((n_years, n_months, n_target_lat, n_target_lon), np.nan, dtype=np.float32),
# #         coords={
# #             'year': years_unique,
# #             'month': months_unique,
# #             'latitude': target_latitudes,
# #             'longitude': target_longitudes
# #         },
# #         dims=['year', 'month', 'latitude', 'longitude'],
# #         name='vertical_velocity_500hPa'
# #     )
# #
# #     print("\nProcessing data for each month and interpolating to 1 degree grid...")
# #     # Loop through each time step to get monthly data and interpolate
# #     for time_idx in range(len(ds['valid_time'])):
# #         current_time_val = ds['valid_time'].isel(valid_time=time_idx)
# #         current_year = pd.to_datetime(current_time_val.item()).year
# #         current_month = pd.to_datetime(current_time_val.item()).month
# #
# #         if not (start_year <= current_year <= end_year):
# #             continue
# #
# #         year_idx = current_year - start_year
# #         month_idx = current_month - 1
# #
# #         # Select current month's calculated EIS and 500hPa omega
# #         # Use .copy() to ensure we're working with a separate copy if needed, though xarray usually handles views well.
# #         eis_current_month = EIS.isel(valid_time=time_idx)
# #         omega_current_month = w_500.isel(valid_time=time_idx)
# #
# #         # Interpolate to the target 1x1 degree grid
# #         # .interp handles coordinates based on their values, not just indices, so lat/lon order is important.
# #         eis_interpolated = eis_current_month.interp(
# #             latitude=target_latitudes, longitude=target_longitudes, method='linear'
# #         )
# #         omega_interpolated = omega_current_month.interp(
# #             latitude=target_latitudes, longitude=target_longitudes, method='linear'
# #         )
# #
# #         # Assign interpolated values to the pre-allocated DataArrays
# #         eis_monthly_1deg[year_idx, month_idx, :, :] = eis_interpolated.values
# #         omega_500hPa_monthly_1deg[year_idx, month_idx, :, :] = omega_interpolated.values
# #
# #         if (time_idx + 1) % 12 == 0 or (time_idx + 1) == len(ds['valid_time']):  # Print yearly progress
# #             print(f"  Processed data up to {current_year}/{current_month:02d}")
# #
# #     print("\nAll data processed and re-gridded.")
# #
# #     # --- Save to NetCDF file ---
# #     print(f"\nSaving results to {output_eis_file_path}")
# #
# #     output_ds = xr.Dataset(
# #         {
# #             eis_monthly_1deg.name: eis_monthly_1deg,
# #             omega_500hPa_monthly_1deg.name: omega_500hPa_monthly_1deg
# #         }
# #     )
# #
# #     # Add global attributes
# #     output_ds.attrs[
# #         'description'] = f"ERA5 Monthly Averaged EIS and 500hPa Vertical Velocity, re-gridded to {target_lat_res}°x{target_lon_res}°."
# #     output_ds.attrs['original_source'] = "ERA5 Reanalysis - Copernicus Climate Change Service (C3S)"
# #     output_ds.attrs['data_period'] = f"{start_year}-{end_year}"
# #     output_ds.attrs['processing_date'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
# #     output_ds.attrs[
# #         'eis_calculation_method'] = "Based on Romps (2017) formulas; LCL height is a common approximation (Stull, 1988)."
# #     output_ds.attrs[
# #         'notes'] = "Calculated across all grid points (ocean and land). Future work may include land/ocean masking."
# #
# #     # Add variable attributes
# #     output_ds['eis'].attrs['units'] = 'K'
# #     output_ds['eis'].attrs['long_name'] = 'Estimated Inversion Strength'
# #     output_ds['vertical_velocity_500hPa'].attrs['units'] = 'Pa/s'
# #     output_ds['vertical_velocity_500hPa'].attrs['long_name'] = 'Vertical Velocity at 500 hPa (Pressure Coordinate)'
# #
# #     # Specify compression for smaller file size
# #     comp = dict(zlib=True, complevel=5)
# #     encoding = {var: comp for var in output_ds.data_vars}
# #
# #     output_ds.to_netcdf(output_eis_file_path, format="NETCDF4", encoding=encoding)
# #
# #     print(f"Processing complete! Results saved to {output_eis_file_path}")
# #     print(
# #         f"Output data dimensions: Years({n_years}) x Months({n_months}) x Latitudes({n_target_lat}) x Longitudes({n_target_lon})")
# #
# # except FileNotFoundError:
# #     print(f"Error: Input file not found at {input_era5_file_path}. Please check the path.")
# # except KeyError as e:
# #     print(
# #         f"Error: Required variable or pressure level not found: {e}. Please check the ERA5 file content (e.g., 't', 'r', 'z', 'w' or pressure_level values).")
# # except Exception as e:
# #     print(f"An unexpected error occurred during processing: {e}")
#
# import xarray as xr
# import numpy as np
#
# # =============================================================================
# # --- 核心物理计算模块 (手动实现) ---
# # =============================================================================
# G = 9.80665
# CP = 1005.7
# RD = 287.058
# RV = 461.5
# LV = 2.501e6
# EPSILON = RD / RV
#
#
# def potential_temperature(pressure_pa, temperature_k):
#     p0 = 100000.0
#     return temperature_k * (p0 / pressure_pa) ** (RD / CP)
#
#
# def saturation_vapor_pressure(temperature_k):
#     temp_c = temperature_k - 273.15
#     return 611.2 * np.exp(17.67 * temp_c / (temp_c + 243.5))
#
#
# def saturation_mixing_ratio(pressure_pa, temperature_k):
#     es = saturation_vapor_pressure(temperature_k)
#     es = np.minimum(es, pressure_pa * 0.99)
#     return EPSILON * es / (pressure_pa - es)
#
#
# def lcl_height_approx(temperature_k, rh_fraction):
#     temp_c = temperature_k - 273.15
#     rh_percent = rh_fraction * 100
#     dew_point_c = temp_c - ((100 - rh_percent) / 5.0)
#     z_lcl = 125.0 * (temp_c - dew_point_c)
#     return z_lcl
#
#
# def moist_lapse_rate_manual(pressure_pa, temperature_k):
#     qs = saturation_mixing_ratio(pressure_pa, temperature_k)
#     numerator = 1 + (LV * qs) / (RD * temperature_k)
#     denominator = 1 + (LV ** 2 * qs) / (CP * RV * temperature_k ** 2)
#     gamma_m = (G / CP) * (1 - numerator / denominator)
#     return gamma_m
#
#
# # =============================================================================
# # --- 主处理函数 ---
# # =============================================================================
# def process_era5_data(input_file_path, output_file_path):
#     print(f"开始处理文件: {input_file_path}")
#     print("采用最终稳健计算方法，请耐心等待...")
#
#     with xr.open_dataset(input_file_path) as ds:
#         # --- 步骤 1 & 2: 重采样并提取变量 ---
#         new_lat = np.arange(90, -90.1, -1.0)
#         new_lon = np.arange(0, 360.1, 1.0)
#         print("第1/7步：将数据重采样至1x1度网格...")
#         ds_resampled = ds.interp(latitude=new_lat, longitude=new_lon, method='linear')
#         print("重采样完成。")
#
#         print("第2/7步：提取所需变量...")
#         # 根据检查结果，数据是3维的，没有'time'维度
#         t_sfc_np = ds_resampled['t'].sel(pressure_level=1000).values
#         t_700_np = ds_resampled['t'].sel(pressure_level=700).values
#         z_700_geopotential = ds_resampled['z'].sel(pressure_level=700).values
#         omega_np = ds_resampled['w'].sel(pressure_level=500).values
#         print("变量提取完成。")
#
#         # --- 步骤 3 - 7: 执行所有计算 ---
#         print("第3/7步：计算低层大气稳定度 (LTS)...")
#         p_sfc_pa = 100000.0
#         p_700_pa = 70000.0
#         theta_sfc = potential_temperature(p_sfc_pa, t_sfc_np)
#         theta_700 = potential_temperature(p_700_pa, t_700_np)
#         lts = theta_700 - theta_sfc
#         print("LTS计算完成。")
#
#         print("第4/7步：计算700hPa高度...")
#         z_700 = z_700_geopotential / G
#         print("700hPa高度计算完成。")
#
#         print("第5/7步：计算抬升凝结高度 (ZLCL)...")
#         z_lcl = lcl_height_approx(t_sfc_np, 0.8)
#         print("ZLCL计算完成。")
#
#         print("第6/7步：计算湿绝热递减率 (Gamma_m)...")
#         t_mean = (t_sfc_np + t_700_np) / 2
#         p_mid_pa = 85000.0
#         gamma_m_850 = np.full(t_mean.shape, np.nan)
#         valid_indices = np.isfinite(t_mean)
#         if np.any(valid_indices):
#             gamma_valid = moist_lapse_rate_manual(p_mid_pa, t_mean[valid_indices])
#             gamma_m_850[valid_indices] = gamma_valid
#         print("Gamma_m计算完成。")
#
#         print("第7/7步：计算最终EIS...")
#         height_diff = z_700 - z_lcl
#         eis = lts - gamma_m_850 * height_diff
#         print("EIS计算完成。")
#
#         # --- 【决定性最终修正】根据文件的真实维度结构创建NetCDF文件 ---
#         print("正在创建并保存最终NetCDF文件...")
#
#         # 1. 定义最终文件的坐标 (不再包含不存在的 'time' 坐标)
#         final_coords = {
#             'valid_time': ds_resampled.coords['valid_time'],
#             'latitude': ds_resampled.coords['latitude'],
#             'longitude': ds_resampled.coords['longitude']
#         }
#
#         # 2. 创建一个包含坐标和主要数据变量的 Dataset
#         #    维度 'dims' 也相应地移除了 'time'
#         final_ds = xr.Dataset(
#             data_vars={
#                 'eis': (('valid_time', 'latitude', 'longitude'), eis),
#                 'omega': (('valid_time', 'latitude', 'longitude'), omega_np)
#             },
#             coords=final_coords
#         )
#
#         # 3. 再将 year 和 month 作为新的数据变量添加进去
#         final_ds['year'] = ('valid_time', ds_resampled['valid_time'].dt.year.values)
#         final_ds['month'] = ('valid_time', ds_resampled['valid_time'].dt.month.values)
#
#         # 4. 为变量添加描述性属性
#         final_ds['eis'].attrs = {'units': 'K', 'long_name': 'Estimated Inversion Strength'}
#         final_ds['omega'].attrs = {'units': 'Pa/s', 'long_name': 'Vertical Velocity at 500hPa (Omega)'}
#         final_ds['year'].attrs = {'long_name': 'Year'}
#         final_ds['month'].attrs = {'long_name': 'Month'}
#
#         # 5. 为整个文件添加全局属性
#         final_ds.attrs['description'] = '月平均EIS和500hPa垂直速度(Omega)，1x1度网格'
#         final_ds.attrs['source_file'] = input_file_path
#
#         # 6. 保存到文件
#         final_ds.to_netcdf(output_file_path)
#         print(f"\n成功！整合后的数据已保存至: {output_file_path}")
#
#
# # --- 运行主程序 ---
# if __name__ == '__main__':
#     input_file = r"F:\data\ERA5\7184c7e96123f1bb2da0d9d5af1a900.nc"
#     output_file = r"F:\data\ERA5\ERA5_EIS_Omega_Monthly_1deg_2006-2023.nc"
#
#     try:
#         process_era5_data(input_file, output_file)
#     except Exception as e:
#         print(f"\n处理过程中发生错误: {e}")

import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path


def plot_global_eis(input_file_path):
    """
    读取计算好的EIS文件，绘制全球多年平均EIS分布图。
    """
    print(f"正在读取文件: {input_file_path}")

    # --- 1. 读取数据 ---
    try:
        ds = xr.open_dataset(input_file_path)
        eis_data = ds['eis']
    except FileNotFoundError:
        print(f"错误：文件未找到，请确认路径是否正确: {input_file_path}")
        return
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return

    # --- 2. 计算多年平均 ---
    print("正在计算多年平均值...")
    # 使用 .mean() 函数沿着 'valid_time' 维度计算平均
    # The archived final data uses ``time``; retained legacy files use
    # ``valid_time``. Support either convention without hard-coded paths.
    time_dim = 'time' if 'time' in eis_data.dims else 'valid_time'
    eis_mean = eis_data.mean(dim=time_dim)
    print("计算完成。")

    # --- 3. 开始绘图 ---
    print("正在生成全球地图...")

    # 创建一个画布和子图，指定地图投影为 PlateCarree (最常用的等经纬度投影)
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # 设置绘图范围为全球
    ax.set_global()

    # --- 4. 绘制数据 ---
    # 使用 xarray 自带的 .plot() 函数，它能很好地与 cartopy 配合
    # vmin, vmax: 设置颜色条的最小值和最大值，对应您的 -10 到 10 的要求
    # cmap: 选择一个色板，'RdBu_r' 是一个很好的红-白-蓝色板，适合展示正负偏差
    # transform: 告诉cartopy我们的数据坐标是标准的经纬度
    im = eis_mean.plot(ax=ax,
                       vmin=-10, vmax=10,
                       cmap='RdBu_r',
                       transform=ccrs.PlateCarree(),
                       add_colorbar=False)  # 我们将手动创建颜色条以更好地控制

    # --- 5. 添加地图元素 ---
    # 添加海岸线
    ax.coastlines()

    # 添加经纬度网格线，并标注在左边和底边
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # 添加颜色条 (colorbar)
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', pad=0.05, shrink=0.8)
    cbar.set_label('多年平均估算逆温强度 (EIS) [K]')

    # 添加标题
    ax.set_title('全球多年平均EIS分布图 (2006-2023)', fontsize=16)

    # --- 6. 显示图像 ---
    print("绘图完成，正在显示图像...")
    plt.show()


# --- 运行主程序 ---
if __name__ == '__main__':
    # 我们刚刚创建的文件的路径
    project_dir = Path(__file__).resolve().parents[2]
    final_file = project_dir / "data" / "LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc"
    plot_global_eis(final_file)
