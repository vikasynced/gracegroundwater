"""
validation_diagnostics_b5.py — Post-B5b Sensitivity & Diagnostics
==================================================================
Three analyses that the validation narrative requires before submission:

  1. Sy sensitivity sweep (0.05, 0.08, 0.10, 0.12, 0.15, 0.20)
     → How does R² / RMSE change with specific yield?
     → Is there an optimal Sy, or is R² stable across a reasonable range?

  2. Pixel-wise temporal R²
     → Compute R² per GRACE pixel (GWS vs wells over time)
     → Report median, IQR, and fraction of pixels with R² > 0.5
     → Tests whether pooled R² = 0.28 hides stronger within-pixel dynamics

  3. RMSE / σ_GWS ratio (signal-to-noise diagnostic)
     → RMSE relative to the standard deviation of GRACE GWS
     → RMSE/σ < 1.0 = validation has predictive skill above climatology

Input:  01_raw_data/cgwb/1_India_GWLs_2000_2024_wells_within_India.csv
        02_processed/gws_monthly.csv
Output: 03_outputs/validation_sy_sensitivity.csv
        03_outputs/validation_pixelwise_r2.csv
        03_outputs/validation_signal_noise.csv
        05_figures/validation_sy_sensitivity.png       (multi-panel)
        05_figures/validation_pixelwise_r2_hist.png    (histogram)
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
BASE = Path(r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER")
CGWB_FILE = BASE / "01_raw_data" / "cgwb" / "1_India_GWLs_2000_2024_wells_within_India.csv"
GWS_FILE = BASE / "02_processed" / "gws_monthly.csv"
OUT_DIR = BASE / "03_outputs"
FIG_DIR = BASE / "05_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Parameters
# ---------------------------------------------------------------------------
BASELINE_START = 2004
BASELINE_END = 2009
MIN_OBS = 20
GRACE_RES = 0.25
QUARTERLY_MONTHS = [1, 5, 8, 11]

SY_VALUES = [0.05, 0.08, 0.10, 0.12, 0.15, 0.20]

ZONES_BBOX = {
    "Z1": {"lat_min": 27, "lat_max": 32, "lon_min": 73, "lon_max": 78},
    "Z2": {"lat_min": 24, "lat_max": 28, "lon_min": 80, "lon_max": 89},
}

LAT_COL = "Latitude"
LON_COL = "Longitude"
WELL_TYPE_COL = "Type of Well"
AQUIFER_COL = "Aquifer Type"

MONTH_ABBR = ["Jan", "May", "Aug", "Nov"]
MONTH_MAP = {"Jan": 1, "May": 5, "Aug": 8, "Nov": 11}

# ===========================================================================
# PART A: SHARED DATA LOADING (do once, reuse for all analyses)
# ===========================================================================
print("=" * 70)
print("VALIDATION DIAGNOSTICS — Sy Sensitivity + Pixel-wise R² + RMSE/σ")
print("=" * 70)

print("\n[LOAD] Loading and filtering CGWB data ...")

df_raw = pd.read_csv(CGWB_FILE, low_memory=False)

# Apply all 4 filters (same as B5b)
df = df_raw.copy()
in_z1 = (
    (df[LAT_COL] >= ZONES_BBOX["Z1"]["lat_min"]) & (df[LAT_COL] <= ZONES_BBOX["Z1"]["lat_max"]) &
    (df[LON_COL] >= ZONES_BBOX["Z1"]["lon_min"]) & (df[LON_COL] <= ZONES_BBOX["Z1"]["lon_max"])
)
in_z2 = (
    (df[LAT_COL] >= ZONES_BBOX["Z2"]["lat_min"]) & (df[LAT_COL] <= ZONES_BBOX["Z2"]["lat_max"]) &
    (df[LON_COL] >= ZONES_BBOX["Z2"]["lon_min"]) & (df[LON_COL] <= ZONES_BBOX["Z2"]["lon_max"])
)
df = df[in_z1 | in_z2].copy()
df["zone"] = np.where(in_z1[df.index], "Z1", "Z2")

bore_mask = df[WELL_TYPE_COL].str.lower().str.contains("bore|tube", na=False)
df = df[bore_mask]

unconf_mask = df[AQUIFER_COL].str.lower().str.contains("unconfined", na=False)
df = df[unconf_mask]

def is_quarterly_col(colname):
    colname = str(colname).strip()
    for m in MONTH_ABBR:
        if colname.startswith(m + "-"):
            suffix = colname[len(m)+1:]
            if suffix.isdigit() and len(suffix) in [2, 4]:
                return True
    return False

depth_cols = [c for c in df.columns if is_quarterly_col(c)]
df["n_obs"] = df[depth_cols].notna().sum(axis=1)
df = df[df["n_obs"] >= MIN_OBS].copy()

print(f"  Filtered wells: {len(df)} (Z1: {(df['zone']=='Z1').sum()}, Z2: {(df['zone']=='Z2').sum()})")

# Parse depth columns to long format
def parse_quarterly_col(colname):
    colname = str(colname).strip()
    for m_abbr, m_int in MONTH_MAP.items():
        if colname.startswith(m_abbr + "-"):
            yr_str = colname[len(m_abbr)+1:]
            if yr_str.isdigit():
                yr = int(yr_str)
                if yr < 100:
                    yr = 2000 + yr
                return m_int, yr
    return None, None

depth_records = []
for col in depth_cols:
    month, year = parse_quarterly_col(col)
    if month is None:
        continue
    sub = df[[LAT_COL, LON_COL, "zone", col]].copy()
    sub.columns = ["latitude", "longitude", "zone", "depth_m"]
    sub["year"] = year
    sub["month"] = month
    sub = sub.dropna(subset=["depth_m"])
    if len(sub) > 0:
        depth_records.append(sub)

df_long = pd.concat(depth_records, ignore_index=True)

# Compute baseline per well
well_ids = df_long.groupby(["latitude", "longitude"])
baseline_depth = well_ids.apply(
    lambda g: g[(g["year"] >= BASELINE_START) & (g["year"] <= BASELINE_END)]["depth_m"].mean()
).reset_index()
baseline_depth.columns = ["latitude", "longitude", "baseline_depth_m"]
df_long = df_long.merge(baseline_depth, on=["latitude", "longitude"], how="left")
df_long = df_long.dropna(subset=["baseline_depth_m"])  # keep only wells with baseline

print(f"  Long-format records with baseline: {len(df_long):,}")

# Load GRACE
print("  Loading GRACE GWS ...")
gws = pd.read_csv(GWS_FILE, parse_dates=["date"])
gws["year"] = gws["date"].dt.year
gws["month"] = gws["date"].dt.month
gws = gws[gws["month"].isin(QUARTERLY_MONTHS)].copy()

# GRACE grid centres
grace_lats = np.sort(gws["lat"].unique())
grace_lons = np.sort(gws["lon"].unique())

def snap_to_grace_grid(lat, lon):
    lat_idx = np.argmin(np.abs(grace_lats - lat))
    lon_idx = np.argmin(np.abs(grace_lons - lon))
    return grace_lats[lat_idx], grace_lons[lon_idx]

df_long["grace_lat"], df_long["grace_lon"] = zip(*df_long.apply(
    lambda row: snap_to_grace_grid(row["latitude"], row["longitude"]), axis=1
))

# ===========================================================================
# PART 1: Sy SENSITIVITY SWEEP
# ===========================================================================
print("\n" + "=" * 70)
print("PART 1: SPECIFIC YIELD SENSITIVITY SWEEP")
print("=" * 70)

sy_results = []

for sy in SY_VALUES:
    print(f"\n  Testing Sy = {sy:.2f} ...")

    # Compute GWS anomaly with this Sy
    df_temp = df_long.copy()
    df_temp["gws_cm_well"] = -(df_temp["depth_m"] - df_temp["baseline_depth_m"]) * sy * 100

    # Aggregate to pixel-month
    well_pixel = (
        df_temp.groupby(["grace_lat", "grace_lon", "year", "month", "zone"])
        .agg(gws_cm_well_mean=("gws_cm_well", "mean"), n_wells=("gws_cm_well", "count"))
        .reset_index()
    )

    # Merge with GRACE
    matched = well_pixel.merge(
        gws, left_on=["grace_lat", "grace_lon", "year", "month"],
        right_on=["lat", "lon", "year", "month"], how="inner"
    )
    matched["gws_cm_grace"] = matched["gws_cm"]
    matched = matched.drop(columns=["gws_cm", "date"], errors="ignore")

    # Annual means
    annual = (
        matched.groupby(["grace_lat", "grace_lon", "year", "zone"])
        .agg(gws_cm_grace=("gws_cm_grace", "mean"),
             gws_cm_well_mean=("gws_cm_well_mean", "mean"),
             n_quarters=("month", "nunique"))
        .reset_index()
    )
    annual = annual[annual["n_quarters"] == 4].copy()

    # Compute metrics per zone
    for zid in ["Z1", "Z2", "Z1+Z2"]:
        if zid == "Z1+Z2":
            zdata = annual
        else:
            zdata = annual[annual["zone"] == zid]

        if len(zdata) < 10:
            sy_results.append({"Sy": sy, "zone": zid, "n_pairs": len(zdata),
                               "R2": np.nan, "RMSE_cm": np.nan, "slope": np.nan, "intercept": np.nan})
            continue

        mask = ~np.isnan(zdata["gws_cm_grace"]) & ~np.isnan(zdata["gws_cm_well_mean"])
        x = zdata.loc[mask, "gws_cm_grace"]
        y = zdata.loc[mask, "gws_cm_well_mean"]

        if len(x) < 10:
            sy_results.append({"Sy": sy, "zone": zid, "n_pairs": len(x),
                               "R2": np.nan, "RMSE_cm": np.nan, "slope": np.nan, "intercept": np.nan})
            continue

        r, p = stats.pearsonr(x, y)
        rmse = np.sqrt(np.mean((x - y) ** 2))
        slope, intercept, _, _, _ = stats.linregress(x, y)

        sy_results.append({
            "Sy": sy, "zone": zid, "n_pairs": len(x),
            "R2": r**2, "RMSE_cm": rmse, "slope": slope, "intercept": intercept,
        })

df_sy = pd.DataFrame(sy_results)
sy_path = OUT_DIR / "validation_sy_sensitivity.csv"
df_sy.to_csv(sy_path, index=False)
print(f"\n  Saved: {sy_path}")

# Print summary table
print(f"\n  {'Sy':<8} {'Zone':<8} {'N pairs':>10} {'R²':>10} {'RMSE':>10} {'Slope':>10}")
print(f"  {'─' * 60}")
for _, row in df_sy.iterrows():
    print(f"  {row['Sy']:<8.2f} {row['zone']:<8} {row['n_pairs']:>10.0f} "
          f"{row['R2']:>10.4f} {row['RMSE_cm']:>10.1f} {row['slope']:>10.3f}")

# ===========================================================================
# PART 2: PIXEL-WISE TEMPORAL R²
# ===========================================================================
print("\n" + "=" * 70)
print("PART 2: PIXEL-WISE TEMPORAL R²")
print("=" * 70)

# Use default Sy = 0.10 for this analysis
sy_default = 0.10
df_long["gws_cm_well"] = -(df_long["depth_m"] - df_long["baseline_depth_m"]) * sy_default * 100

well_pixel = (
    df_long.groupby(["grace_lat", "grace_lon", "year", "month", "zone"])
    .agg(gws_cm_well_mean=("gws_cm_well", "mean"), n_wells=("gws_cm_well", "count"))
    .reset_index()
)

matched = well_pixel.merge(
    gws, left_on=["grace_lat", "grace_lon", "year", "month"],
    right_on=["lat", "lon", "year", "month"], how="inner"
)
matched["gws_cm_grace"] = matched["gws_cm"]

annual = (
    matched.groupby(["grace_lat", "grace_lon", "year", "zone"])
    .agg(gws_cm_grace=("gws_cm_grace", "mean"),
         gws_cm_well_mean=("gws_cm_well_mean", "mean"),
         n_quarters=("month", "nunique"))
    .reset_index()
)
annual = annual[annual["n_quarters"] == 4].copy()

# Compute R² per pixel
pixel_r2 = []
for (glat, glon), grp in annual.groupby(["grace_lat", "grace_lon"]):
    if len(grp) < 5:
        continue
    mask = ~np.isnan(grp["gws_cm_grace"]) & ~np.isnan(grp["gws_cm_well_mean"])
    x = grp.loc[mask, "gws_cm_grace"]
    y = grp.loc[mask, "gws_cm_well_mean"]
    if len(x) < 5:
        continue
    r, p = stats.pearsonr(x, y)
    zone = grp["zone"].mode()[0] if len(grp["zone"].mode()) > 0 else "unknown"
    # Also compute per-pixel RMSE and GRACE std dev
    rmse_pix = np.sqrt(np.mean((x - y) ** 2))
    sigma_gws = np.std(x)
    pixel_r2.append({
        "grace_lat": glat, "grace_lon": glon, "zone": zone,
        "n_years": len(x), "R2": r**2, "pearson_r": r,
        "RMSE_cm": rmse_pix, "sigma_gws_cm": sigma_gws,
        "RMSE_over_sigma": rmse_pix / sigma_gws if sigma_gws > 0 else np.nan,
    })

df_pix = pd.DataFrame(pixel_r2)
pix_path = OUT_DIR / "validation_pixelwise_r2.csv"
df_pix.to_csv(pix_path, index=False)
print(f"  Saved: {pix_path}")
print(f"  Pixels with ≥5 annual pairs: {len(df_pix)}")

# Summary statistics
print(f"\n  Pixel-wise R² summary:")
for zid in ["Z1", "Z2"]:
    zdat = df_pix[df_pix["zone"] == zid]
    if len(zdat) == 0:
        continue
    median_r2 = np.median(zdat["R2"])
    q25 = np.percentile(zdat["R2"], 25)
    q75 = np.percentile(zdat["R2"], 75)
    frac_above_05 = (zdat["R2"] > 0.5).mean()
    frac_positive = (zdat["R2"] > 0).mean()
    print(f"    {zid}: n={len(zdat)}, median R²={median_r2:.3f} (IQR {q25:.3f}–{q75:.3f}), "
          f"R²>0.5: {frac_above_05:.1%}, R²>0: {frac_positive:.1%}")

# Combined
median_r2_all = np.median(df_pix["R2"])
q25_all = np.percentile(df_pix["R2"], 25)
q75_all = np.percentile(df_pix["R2"], 75)
frac_above_05_all = (df_pix["R2"] > 0.5).mean()
frac_positive_all = (df_pix["R2"] > 0).mean()
print(f"    Combined: n={len(df_pix)}, median R²={median_r2_all:.3f} "
      f"(IQR {q25_all:.3f}–{q75_all:.3f}), R²>0.5: {frac_above_05_all:.1%}, "
      f"R²>0: {frac_positive_all:.1%}")

# ===========================================================================
# PART 3: RMSE / σ_GWS RATIO
# ===========================================================================
print("\n" + "=" * 70)
print("PART 3: RMSE / σ_GWS SIGNAL-TO-NOISE DIAGNOSTIC")
print("=" * 70)

# Use the pixel-wise results already computed
sn_rows = []
for zid in ["Z1", "Z2", "Z1+Z2"]:
    if zid == "Z1+Z2":
        zdat = df_pix
    else:
        zdat = df_pix[df_pix["zone"] == zid]

    if len(zdat) == 0:
        continue

    # Pool all pixel data for this zone
    z_annual = annual[annual["zone"] == zid] if zid != "Z1+Z2" else annual
    mask = ~np.isnan(z_annual["gws_cm_grace"]) & ~np.isnan(z_annual["gws_cm_well_mean"])
    sigma_gws = np.std(z_annual.loc[mask, "gws_cm_grace"])
    rmse_pooled = np.sqrt(np.mean(
        (z_annual.loc[mask, "gws_cm_grace"] - z_annual.loc[mask, "gws_cm_well_mean"])**2
    ))

    # Also compute median per-pixel RMSE/σ
    median_rmse_over_sigma = np.nanmedian(zdat["RMSE_over_sigma"])

    sn_rows.append({
        "zone": zid,
        "sigma_gws_cm": sigma_gws,
        "RMSE_pooled_cm": rmse_pooled,
        "RMSE_over_sigma_pooled": rmse_pooled / sigma_gws if sigma_gws > 0 else np.nan,
        "median_RMSE_over_sigma_pixel": median_rmse_over_sigma,
        "n_pixels": len(zdat),
        "n_annual_pairs": len(z_annual) if zid != "Z1+Z2" else len(annual),
    })

df_sn = pd.DataFrame(sn_rows)
sn_path = OUT_DIR / "validation_signal_noise.csv"
df_sn.to_csv(sn_path, index=False)
print(f"  Saved: {sn_path}")

print(f"\n  {'Zone':<10} {'σ_GWS':>10} {'RMSE':>10} {'RMSE/σ':>10} {'Median pix RMSE/σ':>20} {'Interpretation'}")
print(f"  {'─' * 85}")
for _, row in df_sn.iterrows():
    ratio = row["RMSE_over_sigma_pooled"]
    if ratio < 0.5:
        interp = "Excellent — signal dominates"
    elif ratio < 1.0:
        interp = "Good — signal exceeds noise"
    elif ratio < 1.5:
        interp = "Marginal — noise comparable to signal"
    else:
        interp = "Poor — noise dominates signal"
    print(f"  {row['zone']:<10} {row['sigma_gws_cm']:>10.1f} {row['RMSE_pooled_cm']:>10.1f} "
          f"{ratio:>10.2f} {row['median_RMSE_over_sigma_pixel']:>20.2f}   {interp}")

# ===========================================================================
# FIGURE 1: Sy Sensitivity Multi-Panel
# ===========================================================================
print("\n[FIGURE] Generating Sy sensitivity plot ...")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

sy_pivot_r2 = df_sy.pivot(index="Sy", columns="zone", values="R2")
sy_pivot_rmse = df_sy.pivot(index="Sy", columns="zone", values="RMSE_cm")
sy_pivot_slope = df_sy.pivot(index="Sy", columns="zone", values="slope")

colors = {"Z1": "#D62728", "Z2": "#1F77B4", "Z1+Z2": "#2CA02C"}
markers = {"Z1": "o", "Z2": "s", "Z1+Z2": "D"}

# Panel 1: R² vs Sy
ax = axes[0]
for zid in ["Z1", "Z2", "Z1+Z2"]:
    if zid in sy_pivot_r2.columns:
        ax.plot(sy_pivot_r2.index, sy_pivot_r2[zid], marker=markers[zid],
                color=colors[zid], linewidth=1.5, markersize=7, label=zid)
ax.axhline(y=0.28, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, label="Sy=0.10 baseline")
ax.set_xlabel("Specific Yield (Sy)", fontsize=10)
ax.set_ylabel("R²", fontsize=10)
ax.set_title("R² Sensitivity to Specific Yield", fontsize=11, fontweight="bold")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=0)

# Panel 2: RMSE vs Sy
ax = axes[1]
for zid in ["Z1", "Z2", "Z1+Z2"]:
    if zid in sy_pivot_rmse.columns:
        ax.plot(sy_pivot_rmse.index, sy_pivot_rmse[zid], marker=markers[zid],
                color=colors[zid], linewidth=1.5, markersize=7, label=zid)
ax.set_xlabel("Specific Yield (Sy)", fontsize=10)
ax.set_ylabel("RMSE (cm)", fontsize=10)
ax.set_title("RMSE Sensitivity to Specific Yield", fontsize=11, fontweight="bold")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 3: Regression slope vs Sy (should be near 1.0 if well-calibrated)
ax = axes[2]
for zid in ["Z1", "Z2", "Z1+Z2"]:
    if zid in sy_pivot_slope.columns:
        ax.plot(sy_pivot_slope.index, sy_pivot_slope[zid], marker=markers[zid],
                color=colors[zid], linewidth=1.5, markersize=7, label=zid)
ax.axhline(y=1.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, label="Ideal slope = 1.0")
ax.set_xlabel("Specific Yield (Sy)", fontsize=10)
ax.set_ylabel("Regression Slope", fontsize=10)
ax.set_title("Slope Sensitivity to Specific Yield", fontsize=11, fontweight="bold")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig_path1 = FIG_DIR / "validation_sy_sensitivity.png"
fig.savefig(fig_path1, dpi=300, bbox_inches="tight", facecolor="white")
print(f"  Saved: {fig_path1}")
plt.close()

# ===========================================================================
# FIGURE 2: Pixel-wise R² Histogram
# ===========================================================================
print("[FIGURE] Generating pixel-wise R² histogram ...")

fig, ax = plt.subplots(figsize=(8, 5))

for zid, color in [("Z1", "#D62728"), ("Z2", "#1F77B4")]:
    zdat = df_pix[df_pix["zone"] == zid]
    if len(zdat) > 0:
        ax.hist(zdat["R2"], bins=20, alpha=0.5, color=color, label=f"{zid} (n={len(zdat)})",
                edgecolor="white", linewidth=0.5)

# Combined
ax.hist(df_pix["R2"], bins=20, alpha=0.3, color="gray",
        label=f"Combined (n={len(df_pix)})", edgecolor="black", linewidth=0.5, histtype="step")

# Median line
ax.axvline(x=median_r2_all, color="black", linestyle="--", linewidth=1.2,
           label=f"Median R² = {median_r2_all:.3f}")

ax.set_xlabel("Pixel-wise R² (temporal correlation per GRACE pixel)", fontsize=11)
ax.set_ylabel("Number of pixels", fontsize=11)
ax.set_title("Distribution of Pixel-Wise Temporal R²\n(GRACE GWS vs CGWB bore/tube wells, 525 wells across 123 pixels)",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="upper right")
ax.grid(True, alpha=0.3, axis="y")

# Add annotation for fraction > 0.5
ax.text(0.95, 0.90, f"R² > 0.5: {frac_above_05_all:.0%} of pixels",
        transform=ax.transAxes, fontsize=10, ha="right",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

plt.tight_layout()
fig_path2 = FIG_DIR / "validation_pixelwise_r2_hist.png"
fig.savefig(fig_path2, dpi=300, bbox_inches="tight", facecolor="white")
print(f"  Saved: {fig_path2}")
plt.close()

# ===========================================================================
# FINAL SUMMARY
# ===========================================================================
print("\n" + "=" * 70)
print("VALIDATION DIAGNOSTICS COMPLETE")
print("=" * 70)
print(f"\n  Files saved:")
print(f"    {sy_path}")
print(f"    {pix_path}")
print(f"    {sn_path}")
print(f"    {fig_path1}")
print(f"    {fig_path2}")
print(f"\n  Key findings (copy these into your manuscript draft):")
print(f"    - Sy sensitivity: R² varies from __ to __ across Sy 0.05–0.20")
print(f"    - Median pixel-wise R²: {median_r2_all:.3f} (IQR {q25_all:.3f}–{q75_all:.3f})")
print(f"    - Fraction of pixels with R² > 0.5: {frac_above_05_all:.1%}")
print(f"    - RMSE/σ_GWS ratio: see signal_noise.csv for interpretation")
print(f"\n  Defensibility checklist:")
print(f"    [ ] R² stable across Sy 0.08–0.15 → robustness argument")
print(f"    [ ] Median pixel-wise R² > pooled R² → within-pixel dynamics good")
print(f"    [ ] RMSE/σ_GWS < 1.0 → validation has predictive skill")
print("=" * 70)