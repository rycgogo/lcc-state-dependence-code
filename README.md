# State-dependent marine low-cloud-cover analysis

Python code used to generate Figures 1--3 and Supporting Figures S1--S3 for the accompanying *Geophysical Research Letters* manuscript, "Observational Evidence for the State Dependence of Marine Low-Cloud Sensitivity."

## Data access

The processed monthly analysis dataset on a 2° latitude × 5° longitude grid for January 2007--December
2021 is archived separately on Zenodo: **https://doi.org/10.5281/zenodo.21409565**.
Download the NetCDF file and place it at:

```text
data/LCC_CALIPSO_MODIS_ERA5_Tadv_2degx5deg_monthly_2007-2021.nc
```

The five scripts in `code/figures/` then resolve this path automatically.

## Run the figure scripts

```bash
pip install -r requirements.txt
python code/figures/fig01_lcc_sic_frequency.py
python code/figures/fig02_lcc_omega_tadv.py
python code/figures/fig03_lcc_sensitivity.py
python code/figures/fig_s01_stability_indices.py
python code/figures/fig_s02_s03_lcc_relationships.py
```

Generated figures and CSVs are written to `outputs/`.

## Manuscript figure mapping

| Script | Current manuscript output | Earlier manuscript numbering |
| --- | --- | --- |
| `fig01_lcc_sic_frequency.py` | Figure 1 | Figure 2 |
| `fig02_lcc_omega_tadv.py` | Figure 2 | Figure 3 |
| `fig03_lcc_sensitivity.py` | Figure 3 | Figure 6 |
| `fig_s01_stability_indices.py` | Supporting Figure S1 | Figure 1 |
| `fig_s02_s03_lcc_relationships.py` | Supporting Figures S2--S3 | Figures 4--5 |

## Provenance scripts

Scripts in `code/provenance/` document the processing of CALIPSO, MODIS, and
ERA5 data. They are not required to reproduce the manuscript figures. See
`SOURCE_DATA.md` for the official upstream products, citations, and expected
optional raw-input directories.

## Release and license

This archive corresponds to release `v1.0.0` and is distributed under the MIT
License; see `LICENSE-CODE-MIT.txt`.
