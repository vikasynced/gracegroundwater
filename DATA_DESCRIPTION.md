# Data Description — Zenodo Deposit

**Title:** Processed groundwater storage data supporting: "Aquifer Specific Yield Is Associated with Decadal Groundwater Storage Trajectories Across India's Major Hydrogeological Zones"

**Author:** Vikash Kumar, Independent Researcher, New Delhi, India
**Contact:** vikasynced@gmail.com
**License:** CC BY 4.0
**Related publication:** Hydrogeology Journal [DOI to be added upon acceptance]
**Related code repository:** [GitHub URL to be added]

---

## Overview

This deposit contains all processed CSV files needed to reproduce the figures, tables, and statistical results in the paper without re-downloading or re-processing the raw GRACE, GLDAS, and CGWB source data. Raw input data are available from their original sources (see README.md in the GitHub repository). All files are plain-text comma-separated values (UTF-8 encoding) with a header row.

---

## Files in This Deposit

### 1. `gws_monthly.csv`
**Size:** ~250 MB | **Rows:** 2,487,764 | **Produced by:** `process_gws.py`

Per-pixel monthly groundwater storage anomaly for India, derived by subtracting GLDAS Noah 2.1 terrestrial water storage components from GRACE/GRACE-FO CSR RL06.3 total water storage anomalies.

| Column | Unit | Description |
|--------|------|-------------|
| `date` | YYYY-MM-DD | First day of the month |
| `lat` | decimal degrees N | GRACE pixel centre latitude (0.25° grid) |
| `lon` | decimal degrees E | GRACE pixel centre longitude (0.25° grid) |
| `gws_cm` | cm equivalent water height | Groundwater storage anomaly relative to 2004–2009 pixel mean. Pixels with \|gws_cm\| > 200 cm clipped (0.05% of data; Himalayan glacier artefacts) |

**Notes:**
- Baseline period: January 2004 – December 2009, per-pixel mean subtracted
- Mission gap June 2017 – June 2018 (GRACE/GRACE-FO transition): no data rows for those months
- 22 duplicate GRACE timestamps removed; later (corrected) observation retained
- Domain: India bounding box 6–38°N, 66–100°E

---

### 2. `gws_zone_timeseries.csv`
**Rows:** 1,386 | **Produced by:** `zone_summaries.py`

Zone-level monthly mean groundwater storage anomaly, aggregated from `gws_monthly.csv` across all pixels within each of the six hydrogeological zone bounding boxes.

| Column | Unit | Description |
|--------|------|-------------|
| `zone_id` | — | Zone identifier (Z1_IGP_NW, Z2_IGP_E, Z3_DECCAN, Z4_CRYST_S, Z5_ARID_NW, Z6_HIMALAYA) |
| `date` | YYYY-MM-DD | First day of the month |
| `mean_gws_cm` | cm | Mean GWS anomaly across all pixels in zone |
| `std_gws_cm` | cm | Standard deviation of GWS anomaly across pixels |
| `n_pixels` | count | Number of GRACE pixels contributing to zone mean |

---

### 3. `gws_zone_trend.csv`
**Rows:** 12 (6 zones × 2 periods) | **Produced by:** `zone_summaries.py`

Zone-level Mann-Kendall trend statistics and Sen's slope for two periods: 2003–2018 (pre-gap) and 2018–2024 (post-gap, July 2018 cutoff placing the entire mission gap in the pre-period).

| Column | Unit | Description |
|--------|------|-------------|
| `zone_id` | — | Zone identifier |
| `zone_name` | — | Full zone name |
| `period` | — | Analysis period (2003-2018 or 2018-2024) |
| `sen_slope_cm_per_year` | cm/yr | Zone-level Sen's slope (monthly data × 12 to annualise) |
| `mk_p_value` | — | Mann-Kendall two-tailed p-value |
| `mk_trend` | — | Mann-Kendall trend direction (increasing/decreasing/no trend) |
| `n_months` | count | Number of months in analysis period (excluding gap) |
| `n_pixels_mean` | count | Mean number of pixels per monthly zone average |

---

### 4. `gws_zone_acceleration.csv`
**Rows:** 6 | **Produced by:** `zone_summaries.py`

Zone-level acceleration: difference between post-2018 and pre-2018 Sen's slope. Positive = storage recovering or depleting more slowly. Negative = depletion accelerating.

| Column | Unit | Description |
|--------|------|-------------|
| `zone_id` | — | Zone identifier |
| `zone_name` | — | Full zone name |
| `slope_2003_2018` | cm/yr | Sen's slope for pre-gap period |
| `slope_2018_2024` | cm/yr | Sen's slope for post-gap period |
| `acceleration_cm_per_year` | cm/yr | Post minus pre slope (positive = recovery/slowdown) |
| `interpretation` | — | Plain-language interpretation flag |

---

### 5. `validation_matched_pairs_igp.csv`
**Produced by:** `validation_b5_igp.py`

Annual matched pairs of GRACE-derived GWS anomaly and CGWB well-derived GWS anomaly, aggregated to GRACE pixel level. Contains only years where all four quarterly observations (Jan, May, Aug, Nov) are present.

| Column | Unit | Description |
|--------|------|-------------|
| `grace_lat` | decimal degrees N | GRACE pixel centre latitude |
| `grace_lon` | decimal degrees E | GRACE pixel centre longitude |
| `year` | YYYY | Calendar year |
| `zone` | — | Zone (Z1 or Z2) |
| `gws_cm_grace` | cm | GRACE-derived GWS anomaly (annual mean of 4 quarters) |
| `gws_cm_well_mean` | cm | Well-derived GWS anomaly (mean across wells in pixel, annual mean of 4 quarters) |
| `n_wells` | count | Mean number of wells contributing to pixel mean |
| `n_quarters` | count | Always 4 (rows with fewer quarters excluded) |

**Well filter criteria applied:** bore wells and tube wells only (dug wells excluded); unconfined aquifer; ≥20 quarterly observations; within IGP bounding boxes (Z1: 27–32°N, 73–78°E; Z2: 24–28°N, 80–89°E). Specific yield = 0.10; baseline = 2004–2009 per-well mean depth.

---

### 6. `validation_summary_igp.csv`
**Rows:** 3 (Z1, Z2, Z1+Z2) | **Produced by:** `validation_b5_igp.py`

Summary validation statistics (pooled across all pixel-years within each zone).

| Column | Unit | Description |
|--------|------|-------------|
| `zone` | — | Zone (Z1, Z2, or Z1+Z2 combined) |
| `n_pairs` | count | Number of pixel-year pairs |
| `pearson_r` | — | Pearson correlation coefficient |
| `R2` | — | Coefficient of determination (= r²). NOTE: this is the pooled R², which is mathematically distinct from the pixel-wise temporal R² reported in the manuscript. See validation_pixelwise_r2.csv. |
| `p_value` | — | Two-tailed p-value for Pearson r |
| `RMSE_cm` | cm | Root mean square error |

---

### 7. `validation_sy_sensitivity.csv`
**Rows:** 18 (6 Sy values × 3 zones) | **Produced by:** `validation_diagnostics_b5.py`

Supports Supplementary Table S3. Validation metrics repeated for six specific yield values (0.05, 0.08, 0.10, 0.12, 0.15, 0.20) to demonstrate that Pearson R² is mathematically invariant to Sy (since Sy scales the well GWS linearly, and R² is invariant under linear scaling).

| Column | Unit | Description |
|--------|------|-------------|
| `Sy` | dimensionless | Specific yield value tested |
| `zone` | — | Zone (Z1, Z2, or Z1+Z2) |
| `n_pairs` | count | Number of pixel-year pairs |
| `R2` | — | Pooled coefficient of determination |
| `RMSE_cm` | cm | Root mean square error |
| `slope` | cm/cm | OLS regression slope (GRACE ~ well) |
| `intercept` | cm | OLS regression intercept |

---

### 8. `validation_pixelwise_r2.csv`
**Rows:** 87 (one per matched GRACE pixel) | **Produced by:** `validation_diagnostics_b5.py`

Supports Supplementary Table S3 and manuscript text. Pixel-wise temporal R²: for each GRACE pixel with ≥5 matched annual pairs, R² is computed from the time series of that single pixel (GRACE vs wells). This is the primary validation metric reported in the abstract and results (median = 0.81 for Z1).

| Column | Unit | Description |
|--------|------|-------------|
| `grace_lat` | decimal degrees N | GRACE pixel centre latitude |
| `grace_lon` | decimal degrees E | GRACE pixel centre longitude |
| `zone` | — | Zone (Z1 or Z2) |
| `n_years` | count | Number of annual paired observations for this pixel |
| `R2` | — | Pixel-wise temporal R² (Pearson r² for this pixel's time series) |
| `pearson_r` | — | Pearson r for this pixel's time series |
| `RMSE_cm` | cm | Per-pixel RMSE |
| `sigma_gws_cm` | cm | Standard deviation of GRACE GWS for this pixel |
| `RMSE_over_sigma` | — | RMSE/σ signal-to-noise ratio (< 1.0 = validation has predictive skill) |

---

### 9. `validation_signal_noise.csv`
**Rows:** 3 | **Produced by:** `validation_diagnostics_b5.py`

Zone-level RMSE/σ signal-to-noise diagnostic. RMSE/σ < 1.0 indicates validation has predictive skill above climatological mean.

| Column | Unit | Description |
|--------|------|-------------|
| `zone` | — | Zone (Z1, Z2, or Z1+Z2) |
| `sigma_gws_cm` | cm | Standard deviation of pooled GRACE GWS |
| `RMSE_pooled_cm` | cm | Pooled RMSE |
| `RMSE_over_sigma_pooled` | — | Pooled RMSE / pooled σ |
| `median_RMSE_over_sigma_pixel` | — | Median of per-pixel RMSE/σ ratios |
| `n_pixels` | count | Number of GRACE pixels |
| `n_annual_pairs` | count | Total annual pixel-year pairs |

---

### 10. `gap_sensitivity.csv`
**Produced by:** `gap_robustness_b8.py`

Supports Supplementary Table S1. Sensitivity of pre/post-2018 trend split to cutoff date choice: January 2018 vs July 2018. Zero sign flips detected across all six zones, confirming the post-2018 spatial pattern is robust to placement of the GRACE/GRACE-FO mission gap.

| Column | Unit | Description |
|--------|------|-------------|
| `zone_id` | — | Zone identifier |
| `cutoff` | — | Cutoff label (Jan 2018 or Jul 2018) |
| `period` | — | Pre or post |
| `sen_slope_cm_per_year` | cm/yr | Sen's slope for this zone-cutoff-period combination |
| `mk_p_value` | — | Mann-Kendall p-value |
| `mk_trend` | — | Mann-Kendall trend direction |
| `n_months` | count | Months in analysis period |

---

### 11. `zone_level_regression_results.csv`
**Rows:** 6 | **Produced by:** `zone_level_regression.py`

Zone-level aggregated statistics used for the Sy–GWS trend regression (Fig 3 and Table 1 extension). This is the input to all zone-level regression analyses.

| Column | Unit | Description |
|--------|------|-------------|
| `zone_id` | — | Zone identifier |
| `specific_yield_typical_pct` | % | Zone-level typical specific yield from CGWB 2020 Table 2.1 |
| `mean_trend` | cm/yr | Zone mean of per-pixel Sen's slopes (2003–2024 full period) |
| `median_trend` | cm/yr | Zone median of per-pixel Sen's slopes |
| `std_trend` | cm/yr | Standard deviation of per-pixel Sen's slopes within zone |
| `n_pixels` | count | Number of GRACE pixels in zone |
| `sem_trend` | cm/yr | Standard error of the zone mean (std / √n_pixels) |

---

### 12. `aquifer_characteristics.csv`
**Rows:** 6 | **Source:** Compiled manually from CGWB Ground Water Resource Assessment 2020, Table 2.1

Zone-level aquifer properties compiled from CGWB 2020. This is the reference table for specific yield values used in the regression.

| Column | Unit | Description |
|--------|------|-------------|
| `zone_id` | — | Zone identifier |
| `zone_name` | — | Full zone name |
| `states` | — | Indian states within zone |
| `aquifer_type` | — | Dominant aquifer lithology |
| `principal_aquifer_system` | — | CGWB Principal Aquifer System code (e.g., AL01, BS01, GR01) |
| `specific_yield_min_pct` | % | Minimum Sy from CGWB 2020 Table 2.1 range |
| `specific_yield_max_pct` | % | Maximum Sy from CGWB 2020 Table 2.1 range |
| `specific_yield_typical_pct` | % | Typical (mid-range) Sy used in analysis |
| `cgwb_2020_table_ref` | — | Row reference in CGWB 2020 Table 2.1 |

---

## How These Files Relate to the Paper

| Paper element | File(s) |
|---------------|---------|
| Table 1 (zone trends) | `gws_zone_trend.csv`, `zone_level_regression_results.csv` |
| Figure 1a (trend map) | `gws_monthly.csv` (re-processed by `FIGURES_B7.py`) |
| Figure 1b (time series) | `gws_zone_timeseries.csv`, `gws_zone_trend.csv` |
| Figure 2 (validation) | `validation_matched_pairs_igp.csv`, `validation_pixelwise_r2.csv` |
| Figure 3 (Sy regression) | `zone_level_regression_results.csv`, `aquifer_characteristics.csv` |
| Supplementary Table S1 (gap robustness) | `gap_sensitivity.csv` |
| Supplementary Table S2 (zone characteristics) | `aquifer_characteristics.csv` |
| Supplementary Table S3 (Sy sensitivity + pixel R²) | `validation_sy_sensitivity.csv`, `validation_pixelwise_r2.csv` |

---

## Validation Status by Zone

| Zone | Validation outcome | Reason |
|------|--------------------|--------|
| Z1 IGP North-West | ✅ Validated (median pixel R² = 0.81) | 457 strict-criterion bore/tube wells |
| Z2 IGP East | ❌ Unvalidated (median pixel R² = 0.06) | Low GRACE signal variance relative to well noise |
| Z3 Deccan Basalt | ❌ Not attempted | Insufficient bore/tube well density |
| Z4 Crystalline South | ❌ Not attempted | Insufficient bore/tube well density |
| Z5 Arid Alluvial NW | ❌ Not attempted | Insufficient bore/tube well density |
| Z6 Himalayan Foothills | ❌ Not attempted | Predominantly dug wells, mountain signal noise |

Trend estimates for unvalidated zones should be interpreted with appropriate caution. See manuscript Discussion section for full treatment of limitations.
