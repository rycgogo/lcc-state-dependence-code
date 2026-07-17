# Upstream data sources and retrieval record

This file documents the third-party source products used to create the
processed dataset in `data/`. Do **not** upload copies of these source products
to the Zenodo record unless their licences and repository policies explicitly
allow redistribution. Cite them in the manuscript and retain this file in the
Zenodo software/data package.

## CALIPSO

- **Product:** CALIPSO Lidar Level 2 1 km Cloud Layer, Version 4-51
  (`CAL_LID_L2_01kmCLay-Standard-V4-51`), CALIOP instrument, HDF4.
- **Coverage used:** January 2007--December 2021.
- **Official landing page:**
  https://asdc.larc.nasa.gov/project/CALIPSO/CAL_LID_L2_01kmCLay-Standard-V4-51_V4-51
- **DOI:** https://doi.org/10.5067/CALIOP/CALIPSO/CAL_LID_L2_01kmCLay-Standard-V4-51
- **Citation:** NASA/LARC/SD/ASDC. (n.d.). *CALIPSO Lidar Level 2 1 km Cloud
  Layer, V4-51* [Dataset]. NASA Langley Atmospheric Science Data Center DAAC.
- **Processing:** `code/provenance/derive_calipso_lcc.py` calculates
  conditional low-cloud cover from the layer fields using
  `LCC = C(low cloud) / [1 - C(opaque high cloud)]` and then aggregates the
  result to the study grid.
- **Access date:** `[replace with the date on which the original granules were
  obtained; YYYY-MM-DD]`.

## MODIS Aqua

- **Product:** MODIS/Aqua Atmosphere L3 Monthly Product, Collection 6.1,
  1-degree global grid (`MYD08_M3.061`), HDF-EOS.
- **Coverage used:** January 2007--December 2021.
- **Official product page:** https://modis-atmos.gsfc.nasa.gov/products/monthly
- **DOI:** https://doi.org/10.5067/MODIS/MYD08_M3.061
- **Requested citation:** Platnick, S., King, M., & Hubanks, P. (2017).
  *MODIS Atmosphere L3 Monthly Product*. NASA MODIS Adaptive Processing
  System, Goddard Space Flight Center.
- **Processing:** `code/provenance/derive_modis_lcc.py` extracts the selected
  cloud field and regrids it to the CALIPSO analysis grid.
- **Access date:** `[replace with the original download date; YYYY-MM-DD]`.

## ERA5

All ERA5 input is monthly averaged reanalysis downloaded from the Copernicus
Climate Data Store (CDS), then subset to January 2007--December 2021 for the
final analysis dataset.

- **Pressure-level monthly means:** air temperature (`t`) at 1000, 925, and
  700 hPa; pressure velocity (`w`) at 700 hPa; and, for derivation of EIS,
  geopotential (`z`) and relative humidity (`r`) at the pressure levels used by
  the EIS formulation.
- **Single-level monthly means:** skin temperature (`skt`), sea-ice
  concentration (`siconc`), 10-m zonal and meridional wind (`u10`, `v10`), and
  land-sea mask (`lsm`).
- **Purpose:** `t`, `z`, and `r` support EIS/LTS/NSS; `w` is the manuscript's
  700-hPa omega; `skt` and `siconc` provide Ts and SIC; `u10` and `v10` are
  used for Tadv.
- **Required manuscript attribution:** “Contains modified Copernicus Climate
  Change Service information [year]. Neither the European Commission nor ECMWF
  is responsible for any use that may be made of the Copernicus information or
  data it contains.”
- **Citation/access date:** copy the current CDS-provided citation for each
  exact catalogue entry used and replace `[CDS access date]`. The ERA5 paper is
  Hersbach et al. (2020), *Quarterly Journal of the Royal Meteorological
  Society*, 146, 1999–2049, https://doi.org/10.1002/qj.3803.

## Scope and retained source files

The public package preserves the final 2.5-degree monthly analysis dataset,
not the multi-gigabyte source files. Local ERA5 source files may be retained
outside the Zenodo package for audit. The historical
`ERA5_EIS_Omega_Monthly_1deg_2006-2023.nc` file is excluded: its omega variable
is documented as 500 hPa and is not the 700-hPa omega used in this study.

The local raw files currently do not include the ERA5 geopotential and relative
humidity fields required to recompute EIS from first principles. This does not
prevent reproduction of the paper figures because the deposited final dataset
already includes EIS, LTS, and NSS. To make preprocessing fully reproducible,
save the CDS request or download record for `z` and `r` alongside this file.
