# India Aquifer Specific Yield — Code Repository

**Paper:** Aquifer Specific Yield Is Associated with Decadal Groundwater Storage Trajectories Across India's Major Hydrogeological Zones

**Author:** Vikash Kumar, Independent Researcher, New Delhi, India
**Contact:** vikasynced@gmail.com
**Target Journal:** Hydrogeology Journal (Springer Nature)

---

## What This Repository Contains

All analysis scripts used to produce the results, figures, and supplementary tables reported in the paper. The analysis uses two decades (2003–2024) of GRACE/GRACE-FO satellite gravimetry, GLDAS Noah 2.1 land surface model output, and CGWB borehole observations to test whether aquifer specific yield is associated with decadal groundwater storage trajectories across India's six major hydrogeological zones.

---

## Scripts — Run in This Exact Order

| Step | Script | Task | Key Output |
|------|--------|------|------------|
| 1 | `process_gws.py` | Isolate groundwater storage from GRACE + GLDAS | `gws_monthly.csv` |
| 2 | `trend_analysis.py` | Per-pixel Mann-Kendall + Sen's slope | `gws_pixel_trend.csv`, `gws_pixel_acceleration.csv` |
| 3 | `zone_summaries.py` | Zone-level monthly summaries and trends | `gws_zone_timeseries.csv`, `gws_zone_trend.csv`, `gws_zone_acceleration.csv` |
| 4 | `validation_b5_igp.py` | Strict-criterion CGWB borehole validation (IGP only) | `validation_matched_pairs_igp.csv`, `validation_summary_igp.csv` |
| 5 | `validation_diagnostics_b5.py` | Sy sensitivity sweep, pixel-wise R², RMSE/σ | `validation_sy_sensitivity.csv`, `validation_pixelwise_r2.csv`, `validation_signal_noise.csv` |
| 6 | `gap_robustness_b8.py` | GRACE/GRACE-FO transition gap robustness test | `gap_sensitivity.csv` |
| 7 | `FIGURES_B7.py` | Publication figures (Fig 1a, 1b, 2a, 2b) | PNG figures at 300 DPI |
| 8 | `zone_level_regression.py` | Zone-level Sy regression + exact permutation + bootstrap | `zone_level_regression_results.csv`, Fig 3 |

**Important:** Steps 1–3 must complete before any downstream script is run. Steps 4–8 can be run independently of each other once Step 3 is complete.

---

## Input Data — Download These First

All input datasets are freely available from their original sources. This repository does not redistribute raw data.

### 1. GRACE/GRACE-FO CSR RL06.3 Mascon Solution
- **Source:** University of Texas Center for Space Research
- **URL:** https://www2.csr.utexas.edu/grace/RL06_mascons.html
- **File:** `CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc`
- **Version:** RL06.3 (released March 2026; includes corrected GRACE-FO accelerometer transplant data post-July 2023)
- **Size:** ~105 MB
- **Place at:** `D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER\01_raw_data\grace_fo\`

### 2. GLDAS Noah 2.1 Monthly
- **Source:** NASA GES DISC
- **URL:** https://disc.gsfc.nasa.gov/datasets/GLDAS_NOAH025_M_2.1
- **Files:** 264 monthly NetCDF4 files, January 2003 – December 2024
- **Naming:** `GLDAS_NOAH025_M.A[YYYYMM].021.nc4`
- **Size:** ~6.3 GB total
- **Place at:** `D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER\01_raw_data\gldas\`
- **Note:** NASA Earthdata Login required (free registration at https://urs.earthdata.nasa.gov)

### 3. CGWB National Groundwater Level Database
- **Source:** Central Ground Water Board, Ministry of Jal Shakti, Government of India
- **URL:** https://cgwb.gov.in (Data → Ground Water Level Data)
- **File:** `1_India_GWLs_2000_2024_wells_within_India.csv`
- **Size:** ~250 MB
- **Place at:** `D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER\01_raw_data\cgwb\`
- **Note:** No account required; freely available

### 4. Aquifer Characteristics
- **Source:** CGWB Ground Water Resource Assessment 2020, Table 2.1 (Norms Recommended for Specific Yield)
- **URL:** https://cgwb.gov.in (Publications → Dynamic Ground Water Resources of India 2020)
- **File:** `aquifer_characteristics.csv` (compiled manually from Table 2.1)
- **Place at:** `D:\INDIA_AQUIFER_STUDY\ROLE_D_INTEGRATION\00_inputs_from_other_roles\`
- **Note:** This compiled file is included in the Zenodo data deposit

---

## Required Folder Structure

```
D:\INDIA_AQUIFER_STUDY\
│
├── ROLE_B_GROUNDWATER\
│   ├── 01_raw_data\
│   │   ├── grace_fo\
│   │   │   └── CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc
│   │   ├── gldas\
│   │   │   └── GLDAS_NOAH025_M.A[YYYYMM].021.nc4  (264 files)
│   │   └── cgwb\
│   │       └── 1_India_GWLs_2000_2024_wells_within_India.csv
│   ├── 02_processed\         ← created by process_gws.py
│   ├── 03_outputs\           ← created by zone_summaries.py and validation scripts
│   └── 05_figures\           ← created by FIGURES_B7.py
│
└── ROLE_D_INTEGRATION\
    ├── 00_inputs_from_other_roles\
    │   └── aquifer_characteristics.csv
    ├── 02_panel_regression\  ← created by zone_level_regression.py
    └── 04_paper_figures\     ← created by zone_level_regression.py
```

All output folders are created automatically by the scripts. Only the raw input folders and files need to be set up manually.

---

## Python Environment

Tested on Python 3.10 (Anaconda base environment), Windows 11.

Install all dependencies with:
```bash
pip install numpy pandas scipy statsmodels matplotlib xarray netCDF4 pymannkendall rasterio
```

Or with conda:
```bash
conda install numpy pandas scipy statsmodels matplotlib xarray netcdf4 rasterio
pip install pymannkendall
```

No GPU or cloud compute required. All analysis runs locally. The most memory-intensive step is `process_gws.py` (loads 264 GLDAS files; ~6 GB RAM peak). All other scripts run in under 1 GB RAM.

**Approximate run times on a mid-range laptop (Intel i7, 8 GB RAM):**

| Script | Approx. run time |
|--------|-----------------|
| `process_gws.py` | 15–30 minutes |
| `trend_analysis.py` | 20–40 minutes |
| `zone_summaries.py` | 2–5 minutes |
| `validation_b5_igp.py` | 5–10 minutes |
| `validation_diagnostics_b5.py` | 10–20 minutes |
| `gap_robustness_b8.py` | < 1 minute |
| `FIGURES_B7.py` | 5–10 minutes |
| `zone_level_regression.py` | < 1 minute |

---

## The Six Hydrogeological Zones

| Zone ID | Name | States | Aquifer Type | Sy (%) |
|---------|------|--------|-------------|--------|
| Z1_IGP_NW | IGP North-West | Punjab, Haryana, W UP, Delhi | Unconfined alluvial | 10 |
| Z2_IGP_E | IGP East | Bihar, E UP, W Bengal | Alluvial | 10 |
| Z3_DECCAN | Deccan Basalt | Maharashtra, Telangana, N Karnataka | Fractured basalt | 2 |
| Z4_CRYST_S | Crystalline South | Tamil Nadu, S Karnataka, S Andhra | Fractured granite/gneiss | 2 |
| Z5_ARID_NW | Arid Alluvial NW | Rajasthan (W), Gujarat (N) | Confined/semi-confined alluvial | 8 |
| Z6_HIMALAYA | Himalayan Foothills | Uttarakhand, HP, Assam | Mountain-front recharge | 8 |

Zone bounding boxes are defined in each script and documented in Supplementary Table S2.

---

## Key Results (for verification)

After running all scripts, the following statistics should be reproduced exactly:

| Statistic | Expected value |
|-----------|---------------|
| Z1 mean GWS trend (2003–2024) | −4.76 ± 0.12 cm/yr |
| Z1 median pixel-wise validation R² | 0.81 (n = 73 pixels) |
| Z2 median pixel-wise validation R² | 0.06 (n = 14 pixels) — FAILED |
| Validation R² invariance across Sy 0.05–0.20 | R² = 0.27–0.28 |
| Zone-level OLS R² (Sy vs GWS trend, n = 6) | 0.44 |
| OLS slope | −0.34 cm/yr per % Sy |
| Spearman ρ | −0.96 |
| Exact permutation p (all 720 orderings) | 0.022 |
| Bootstrap 95% CI for slope | [−1.51, −0.11] |
| Gap robustness: sign flips across cutoffs | 0 (zero) |

---

## Data Deposit

Processed output files (CSVs) are deposited on Zenodo:
**DOI:** [to be added upon acceptance]

Raw input data should be downloaded directly from the original sources listed above.

---

## Citation

Kumar, V. (2026). Aquifer specific yield is associated with decadal groundwater storage trajectories across India's major hydrogeological zones. *Hydrogeology Journal*. [DOI to be added upon acceptance]

---

## License

Code: MIT License
Data: CC BY 4.0 (see Zenodo deposit)
