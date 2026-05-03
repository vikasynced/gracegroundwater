India Aquifer Specific Yield — Code Repository
Paper: Aquifer Specific Yield Is Associated with Decadal Groundwater Storage Trajectories Across India's Major Hydrogeological Zones
Author: Vikash Kumar, Independent Researcher, India
Contact: vikasynced@gmail.com
Journal: Hydrogeology Journal (Springer Nature)
Code DOI: https://doi.org/10.5281/zenodo.19999102
---
What This Repository Contains
All analysis scripts used to produce the results, figures, and supplementary tables reported in the paper. The analysis uses two decades (2003–2024) of GRACE/GRACE-FO satellite gravimetry, GLDAS Noah 2.1 land surface model output, and CGWB borehole observations to test whether aquifer specific yield is associated with decadal groundwater storage trajectories across India's six major hydrogeological zones.
All input data are freely available from their original sources listed below. This repository contains code only — no data files are redistributed.
---
Scripts — Run in This Exact Order
Step	Script	Task	Key Output
1	`01\_process\_gws.py`	Isolate groundwater storage from GRACE + GLDAS	`gws\_monthly.csv`
2	`02\_trend\_analysis.py`	Per-pixel Mann-Kendall + Sen's slope	`gws\_pixel\_trend.csv`, `gws\_pixel\_acceleration.csv`
3	`03\_zone\_summaries.py`	Zone-level monthly summaries and trends	`gws\_zone\_timeseries.csv`, `gws\_zone\_trend.csv`, `gws\_zone\_acceleration.csv`
4	`04\_validation\_b5\_igp.py`	Strict-criterion CGWB borehole validation (IGP only)	`validation\_matched\_pairs\_igp.csv`, `validation\_summary\_igp.csv`
5	`05\_validation\_diagnostics\_b5.py`	Sy sensitivity sweep, pixel-wise R², RMSE/σ	`validation\_sy\_sensitivity.csv`, `validation\_pixelwise\_r2.csv`, `validation\_signal\_noise.csv`
6	`06\_gap\_robustness\_b8.py`	GRACE/GRACE-FO transition gap robustness test	`gap\_sensitivity.csv`
7	`07\_FIGURES\_B7.py`	Publication figures (Fig 1a, 1b, 2a, 2b)	PNG figures at 300 DPI
8	`08\_zone\_level\_regression.py`	Zone-level Sy regression + exact permutation + bootstrap	`zone\_level\_regression\_results.csv`, Fig 3
Important: Steps 1–3 must complete before any downstream script is run. Steps 4–8 can be run independently of each other once Step 3 is complete.
---
Input Data — Download These First
All input datasets are freely available. No account is required except for GLDAS (NASA Earthdata, free registration).
1. GRACE/GRACE-FO CSR RL06.3 Mascon Solution
Source: University of Texas Center for Space Research
URL: https://www2.csr.utexas.edu/grace/RL06_mascons.html
File: `CSR\_GRACE\_GRACE-FO\_RL0603\_Mascons\_all-corrections.nc`
Version: RL06.3 (March 2026 release; includes corrected GRACE-FO accelerometer transplant data post-July 2023)
Size: ~105 MB
Place at: `D:\\INDIA\_AQUIFER\_STUDY\\ROLE\_B\_GROUNDWATER\\01\_raw\_data\\grace\_fo\\`
2. GLDAS Noah 2.1 Monthly
Source: NASA GES DISC
URL: https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1
Files: 264 monthly NetCDF4 files, January 2003 – December 2024
Naming: `GLDAS\_NOAH025\_M.A\[YYYYMM].021.nc4`
Size: ~6.3 GB total
Place at: `D:\\INDIA\_AQUIFER\_STUDY\\ROLE\_B\_GROUNDWATER\\01\_raw\_data\\gldas\\`
Note: Free NASA Earthdata account required — register at https://urs.earthdata.nasa.gov
3. CGWB National Groundwater Level Database
Source: Central Ground Water Board, Ministry of Jal Shakti, Government of India
URL: https://cgwb.gov.in (Data → Ground Water Level Data)
File: `1\_India\_GWLs\_2000\_2024\_wells\_within\_India.csv`
Size: ~250 MB
Place at: `D:\\INDIA\_AQUIFER\_STUDY\\ROLE\_B\_GROUNDWATER\\01\_raw\_data\\cgwb\\`
4. Aquifer Characteristics (CGWB 2020)
Source: CGWB Ground Water Resource Assessment 2020, Table 2.1 (Norms for Specific Yield)
URL: https://cgwb.gov.in (Publications → Dynamic Ground Water Resources of India 2020)
File: `aquifer\_characteristics.csv` — compile manually from Table 2.1 using the six zone Sy values in the paper's Methods section
Place at: `D:\\INDIA\_AQUIFER\_STUDY\\ROLE\_D\_INTEGRATION\\00\_inputs\_from\_other\_roles\\`
---
Required Folder Structure
```
D:\\INDIA\_AQUIFER\_STUDY\\
│
├── ROLE\_B\_GROUNDWATER\\
│   ├── 01\_raw\_data\\
│   │   ├── grace\_fo\\
│   │   │   └── CSR\_GRACE\_GRACE-FO\_RL0603\_Mascons\_all-corrections.nc
│   │   ├── gldas\\
│   │   │   └── GLDAS\_NOAH025\_M.A\[YYYYMM].021.nc4  (264 files)
│   │   └── cgwb\\
│   │       └── 1\_India\_GWLs\_2000\_2024\_wells\_within\_India.csv
│   ├── 02\_processed\\         ← created automatically by 01\_process\_gws.py
│   ├── 03\_outputs\\           ← created automatically by 03\_zone\_summaries.py and validation scripts
│   └── 05\_figures\\           ← created automatically by 07\_FIGURES\_B7.py
│
└── ROLE\_D\_INTEGRATION\\
    ├── 00\_inputs\_from\_other\_roles\\
    │   └── aquifer\_characteristics.csv
    ├── 02\_panel\_regression\\  ← created automatically by 08\_zone\_level\_regression.py
    └── 04\_paper\_figures\\     ← created automatically by 08\_zone\_level\_regression.py
```
All output folders are created automatically by the scripts. Only the raw input folders and files need to be set up manually before running.
---
Python Environment
Tested on Python 3.10 (Anaconda base environment), Windows 11.
```bash
pip install numpy pandas scipy statsmodels matplotlib xarray netCDF4 pymannkendall rasterio
```
Or with conda:
```bash
conda install numpy pandas scipy statsmodels matplotlib xarray netcdf4 rasterio
pip install pymannkendall
```
No GPU or cloud compute required. All analysis runs locally on a standard laptop.
Approximate run times (Intel i7, 8 GB RAM):
Script	Approx. run time
`01\_process\_gws.py`	15–30 minutes
`02\_trend\_analysis.py`	20–40 minutes
`03\_zone\_summaries.py`	2–5 minutes
`04\_validation\_b5\_igp.py`	5–10 minutes
`05\_validation\_diagnostics\_b5.py`	10–20 minutes
`06\_gap\_robustness\_b8.py`	< 1 minute
`07\_FIGURES\_B7.py`	5–10 minutes
`08\_zone\_level\_regression.py`	< 1 minute
The most memory-intensive step is `01\_process\_gws.py` (loads 264 GLDAS files; ~6 GB RAM peak). All other scripts run comfortably within 2 GB RAM.
---
The Six Hydrogeological Zones
Zone ID	Name	States	Aquifer Type	Sy (%)
Z1_IGP_NW	IGP North-West	Punjab, Haryana, W UP, Delhi	Unconfined alluvial	10
Z2_IGP_E	IGP East	Bihar, E UP, W Bengal	Alluvial	10
Z3_DECCAN	Deccan Basalt	Maharashtra, Telangana, N Karnataka	Fractured basalt	2
Z4_CRYST_S	Crystalline South	Tamil Nadu, S Karnataka, S Andhra	Fractured granite/gneiss	2
Z5_ARID_NW	Arid Alluvial NW	Rajasthan (W), Gujarat (N)	Confined/semi-confined alluvial	8
Z6_HIMALAYA	Himalayan Foothills	Uttarakhand, HP, Assam	Mountain-front recharge	8
Zone bounding boxes are defined inside each script and documented in Supplementary Table S2 of the paper.
---
Key Results (for verification)
After running all scripts in order, these statistics should be reproduced exactly:
Statistic	Expected value
Z1 mean GWS trend (2003–2024)	−4.76 ± 0.12 cm/yr
Z1 median pixel-wise validation R²	0.81 (n = 73 pixels)
Z2 median pixel-wise validation R²	0.06 (n = 14 pixels) — validation failed
Validation R² invariance across Sy 0.05–0.20	R² = 0.27–0.28
Zone-level OLS R² (Sy vs GWS trend, n = 6)	0.44
OLS slope	−0.34 cm/yr per % Sy
Spearman ρ	−0.96
Exact permutation p (all 720 orderings)	0.022
Bootstrap 95% CI for slope	[−1.51, −0.11]
Gap robustness: sign flips across cutoffs	0 (zero)
---
Citation
Paper:
Kumar, V. (2026). Aquifer specific yield is associated with decadal groundwater storage trajectories across India's major hydrogeological zones. Hydrogeology Journal. [DOI to be added upon acceptance]
Code:
Kumar, V. (2026). Analysis code — aquifer specific yield and groundwater storage trajectories across India (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.19999102
---
License
MIT License — see LICENSE file for details.
