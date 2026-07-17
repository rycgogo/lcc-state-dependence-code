State-dependent marine low-cloud-cover analysis

Python code and processed data used to generate Figures 1–6 for the manuscript *Observational Evidence for the State Dependence of Marine Low-Cloud Feedback*.

## Overview

This repository documents the analysis of state-dependent relationships between marine low-cloud cover and cloud-controlling factors using CALIPSO and MODIS satellite observations together with ERA5 reanalysis data. The analysis covers January 2007 through December 2021.

The repository contains figure-generation scripts and provenance scripts for the processed analysis dataset. The final processed data and the archival release of the software are preserved on Zenodo.

## Archived versions and licenses

- **Processed analysis dataset:** Ren, Y. (2026). *Processed monthly global marine low-cloud-cover and environmental dataset, 2007–2021* (Version 1.0.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.21409565 (CC-BY-4.0).
- **Archival software release:** Ren, Y. (2026). *Python code for state-dependent marine low-cloud-cover analysis* (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.21407760 (MIT License).
- **Active-development repository:** https://github.com/rycgogo/lcc-state-dependence-code

For reproducibility, use the Zenodo software record corresponding to Version 1.0.0. If this GitHub repository is updated, the exact Git commit or release tag used for a manuscript version should be recorded here.

## Repository structure

```text
.
├── code/
│   ├── figures/
│   │   ├── fig01_stability_indices.py
│   │   ├── fig02_lcc_sic_frequency.py
│   │   ├── fig03_lcc_omega_tadv.py
│   │   ├── fig04_fig05_lcc_relationships.py
│   │   └── fig06_lcc_sensitivity.py
│   └── provenance/
│       ├── derive_calipso_lcc.py
│       ├── derive_modis_lcc.py
│       ├── derive_era5_eis.py
│       └── derive_temperature_advection.py
├── requirements.txt
├── SOURCE_DATA.md
├── LICENSE-CODE-MIT.txt
└── README.md
```

## Software environment

The scripts were developed in Python. Install the required packages in a clean environment:

```bash
python -m venv .venv
source .venv/bin/activate       # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

The required packages are listed in `requirements.txt`, including `xarray`, `numpy`, `pandas`, `matplotlib`, `scipy`, `statsmodels`, `cartopy`, `netCDF4`, `pyhdf`, `tqdm`, and `cfgrib`.

## Data access and setup

1. Download the processed dataset from the Zenodo dataset record above.
2. Create a local `data/` directory in the repository root.
3. Place the downloaded NetCDF file in `data/` using the filename expected by the figure scripts:

```text
data/LCC_CALIPSO_MODIS_ERA5_Tadv_2.5deg_monthly_2007-2021.nc
```

4. The scripts in `code/figures/` resolve this path automatically and write generated figures and CSV files to `outputs/`.

> **Important consistency check before release:** the current manuscript describes a common grid of 2° latitude × 5° longitude, while the archived package filename uses `2.5deg`. Confirm the actual grid in the deposited NetCDF metadata and make the manuscript, filename, and this README consistent before final submission.

## Reproducing the manuscript figures

Run the following scripts from the repository root after the environment and dataset have been prepared:

```bash
python code/figures/fig01_stability_indices.py
python code/figures/fig02_lcc_sic_frequency.py
python code/figures/fig03_lcc_omega_tadv.py
python code/figures/fig04_fig05_lcc_relationships.py
python code/figures/fig06_lcc_sensitivity.py
```

These scripts generate the analyses and graphics for Figures 1–6. Generated files are written to `outputs/`.

## Provenance and source products

The scripts in `code/provenance/` document the construction of the processed data from CALIPSO, MODIS, and ERA5 products. They are not required to reproduce the manuscript figures when the processed Zenodo dataset is used.

Detailed source-product versions, official landing pages, citations, and processing notes are provided in [`SOURCE_DATA.md`](SOURCE_DATA.md). The public archive contains the processed analysis dataset rather than redistributable copies of the multi-gigabyte upstream source products.

## Citation

If you use this code or dataset, please cite the accompanying manuscript and the relevant Zenodo records:

```text
Ren, Y. (2026). Python code for state-dependent marine low-cloud-cover analysis
(Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.21407760

Ren, Y. (2026). Processed monthly global marine low-cloud-cover and environmental
dataset, 2007–2021 (Version 1.0.0) [Dataset]. Zenodo.
https://doi.org/10.5281/zenodo.21409565
```

## License

The code is distributed under the MIT License. The processed dataset is released under CC-BY-4.0. See the Zenodo records and `LICENSE-CODE-MIT.txt` for details.
