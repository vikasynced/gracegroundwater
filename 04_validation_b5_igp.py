"""
validation_b5_igp.py — Task B5b
=================================
IGP-only groundwater validation: GRACE-derived GWS vs CGWB bore/tube wells.
Strict-criterion filters applied to isolate regional aquifer signal.

Filters:
  - Zones: Z1 (lat 27–32N, lon 73–78E) + Z2 (lat 24–28N, lon 80–89E)
  - Well type: Bore well + Tube well ONLY (no dug wells, no piezometers)
  - Aquifer type: Unconfined only
  - Minimum observations: ≥20 quarterly readings
  - Specific yield: Sy = 0.10 (literature default for IGP alluvium)
  - Baseline period: 2004–2009 per well
  - Season: Jan, May, Aug, Nov quarterly

Output:
  03_outputs/validation_matched_pairs_igp.csv  — well-pixel matched seasonal pairs
  03_outputs/validation_summary_igp.csv        — R², Pearson r, RMSE by zone
  05_figures/validation_igp_scatter.png        — publication scatter plot (300 DPI)
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import matplotlib.pyplot as plt
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
SY = 0.10                          # Specific yield (dimensionless)
BASELINE_START = 2004
BASELINE_END = 2009
MIN_OBS = 20                       # Minimum quarterly observations per well
GRACE_RES = 0.25                   # GRACE pixel resolution in degrees

# Zone bounding boxes (section 10.1)
ZONES = {
    "Z1": {"lat_min": 27, "lat_max": 32, "lon_min": 73, "lon_max": 78},
    "Z2": {"lat_min": 24, "lat_max": 28, "lon_min": 80, "lon_max": 89},
}

# Quarterly months
QUARTERLY_MONTHS = [1, 5, 8, 11]  # Jan, May, Aug, Nov

# ---------------------------------------------------------------------------
# 2. Load CGWB well data
# ---------------------------------------------------------------------------
print("=" * 70)
print("B5b: IGP-ONLY GROUNDWATER VALIDATION")
print("=" * 70)

print("\n[1/7] Loading CGWB well data ...")
df_raw = pd.read_csv(CGWB_FILE, low_memory=False)
print(f"  Total wells in CGWB database: {len(df_raw):,}")

# CGWB column names (verified)
LAT_COL = "Latitude"
LON_COL = "Longitude"
WELL_TYPE_COL = "Type of Well"
AQUIFER_COL = "Aquifer Type"

# ---------------------------------------------------------------------------
# 3. Apply strict-criterion filters
# ---------------------------------------------------------------------------
print("\n[2/7] Applying strict-criterion filters ...")

df = df_raw.copy()

# --- Filter 1: IGP bounding boxes ---
print("\n  Filter 1: IGP bounding boxes (Z1 + Z2) ...")
in_z1 = (
    (df[LAT_COL] >= ZONES["Z1"]["lat_min"]) & (df[LAT_COL] <= ZONES["Z1"]["lat_max"]) &
    (df[LON_COL] >= ZONES["Z1"]["lon_min"]) & (df[LON_COL] <= ZONES["Z1"]["lon_max"])
)
in_z2 = (
    (df[LAT_COL] >= ZONES["Z2"]["lat_min"]) & (df[LAT_COL] <= ZONES["Z2"]["lat_max"]) &
    (df[LON_COL] >= ZONES["Z2"]["lon_min"]) & (df[LON_COL] <= ZONES["Z2"]["lon_max"])
)
df = df[in_z1 | in_z2].copy()
df["zone"] = np.where(in_z1[df.index], "Z1", "Z2")
print(f"  Remaining: {len(df):,} wells")

# --- Filter 2: Bore well + Tube well only ---
print("\n  Filter 2: Bore well + Tube well only ...")
bore_mask = df[WELL_TYPE_COL].str.lower().str.contains("bore|tube", na=False)
df = df[bore_mask].copy()
print(f"  Remaining: {len(df):,} wells")

# --- Filter 3: Unconfined aquifer only ---
print("\n  Filter 3: Unconfined aquifer only ...")
unconf_mask = df[AQUIFER_COL].str.lower().str.contains("unconfined", na=False)
df = df[unconf_mask].copy()
print(f"  Remaining: {len(df):,} wells")

# --- Filter 4: >= MIN_OBS quarterly observations ---
print(f"\n  Filter 4: ≥{MIN_OBS} quarterly observations ...")

# Identify quarterly depth columns — format is 'Jan-00', 'May-00', etc.
MONTH_ABBR = ["Jan", "May", "Aug", "Nov"]

def is_quarterly_col(colname):
    """Check if a column name matches 'Jan-00' pattern (month-yy)."""
    colname = str(colname).strip()
    for m in MONTH_ABBR:
        if colname.startswith(m + "-"):
            suffix = colname[len(m)+1:]
            if suffix.isdigit() and len(suffix) in [2, 4]:
                return True
    return False

depth_cols = [c for c in df.columns if is_quarterly_col(c)]
print(f"  Found {len(depth_cols)} quarterly depth columns")
print(f"  Sample columns: {depth_cols[:5]}...")

# Count non-null observations per well
df["n_obs"] = df[depth_cols].notna().sum(axis=1)
df = df[df["n_obs"] >= MIN_OBS].copy()
print(f"  Remaining: {len(df):,} wells")

# Zone-wise breakdown
print(f"\n  Zone breakdown after all filters:")
for zid in ["Z1", "Z2"]:
    n = (df["zone"] == zid).sum()
    print(f"    {zid}: {n} wells")

# ---------------------------------------------------------------------------
# 4. Convert depth-to-water to GWS anomaly per well
# ---------------------------------------------------------------------------
print("\n[3/7] Converting depth-to-water to GWS anomaly (Sy = 0.10) ...")

# Parse column name like 'Jan-00' → month=1, year=2000
MONTH_MAP = {"Jan": 1, "May": 5, "Aug": 8, "Nov": 11}

def parse_quarterly_col(colname):
    """Return (month_int, year_int) or (None, None)."""
    colname = str(colname).strip()
    for m_abbr, m_int in MONTH_MAP.items():
        if colname.startswith(m_abbr + "-"):
            yr_str = colname[len(m_abbr)+1:]
            if yr_str.isdigit():
                yr = int(yr_str)
                if yr < 100:
                    yr = 2000 + yr if yr >= 0 else 1900 + yr
                return m_int, yr
    return None, None

# Build long-format records
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

if not depth_records:
    print("  ERROR: No depth records parsed.")
    raise SystemExit("Cannot proceed — check quarterly column format.")

df_long = pd.concat(depth_records, ignore_index=True)
print(f"  Long-format records: {len(df_long):,}")

# Compute baseline (2004–2009 mean depth) per well
well_ids = df_long.groupby(["latitude", "longitude"])

baseline_depth = well_ids.apply(
    lambda g: g[(g["year"] >= BASELINE_START) & (g["year"] <= BASELINE_END)]["depth_m"].mean()
).reset_index()
baseline_depth.columns = ["latitude", "longitude", "baseline_depth_m"]

# Merge baseline back
df_long = df_long.merge(baseline_depth, on=["latitude", "longitude"], how="left")

# Compute GWS anomaly: dGWS_cm = -(depth - baseline) × Sy × 100
df_long["gws_cm_well"] = -(df_long["depth_m"] - df_long["baseline_depth_m"]) * SY * 100

# Drop records without baseline
n_before = len(df_long)
df_long = df_long.dropna(subset=["gws_cm_well"])
print(f"  Records with valid baseline: {len(df_long):,} (dropped {n_before - len(df_long):,})")

# ---------------------------------------------------------------------------
# 5. Load GRACE GWS and snap coordinates
# ---------------------------------------------------------------------------
print("\n[4/7] Loading GRACE GWS data ...")
gws = pd.read_csv(GWS_FILE, parse_dates=["date"])
gws["year"] = gws["date"].dt.year
gws["month"] = gws["date"].dt.month

# Filter to quarterly months
gws = gws[gws["month"].isin(QUARTERLY_MONTHS)].copy()
print(f"  GRACE quarterly records: {len(gws):,}")

# Get unique GRACE grid centres
grace_lats = np.sort(gws["lat"].unique())
grace_lons = np.sort(gws["lon"].unique())
print(f"  GRACE grid: {len(grace_lats)} lats × {len(grace_lons)} lons")
print(f"  Lat grid centres sample: {grace_lats[:5]}...")
print(f"  Lon grid centres sample: {grace_lons[:5]}...")

# Snap well coordinates to NEAREST GRACE grid centre (not round-to-nearest-0.25)
def snap_to_grace_grid(lat, lon):
    """Find the nearest GRACE grid centre for a given coordinate."""
    lat_idx = np.argmin(np.abs(grace_lats - lat))
    lon_idx = np.argmin(np.abs(grace_lons - lon))
    return grace_lats[lat_idx], grace_lons[lon_idx]

print("  Snapping well coordinates to GRACE grid...")
df_long["grace_lat"], df_long["grace_lon"] = zip(*df_long.apply(
    lambda row: snap_to_grace_grid(row["latitude"], row["longitude"]), axis=1
))

# Verify a few snaps
for i in range(3):
    print(f"    Well ({df_long['latitude'].iloc[i]:.4f}, {df_long['longitude'].iloc[i]:.4f}) → "
          f"GRACE ({df_long['grace_lat'].iloc[i]:.2f}, {df_long['grace_lon'].iloc[i]:.2f})")

# ---------------------------------------------------------------------------
# 6. Match wells to GRACE pixels and compute seasonal annual means
# ---------------------------------------------------------------------------
print("\n[5/7] Matching wells to GRACE pixels ...")

# Aggregate wells per GRACE pixel per (year, month) → mean well GWS
well_pixel = (
    df_long.groupby(["grace_lat", "grace_lon", "year", "month", "zone"])
    .agg(
        gws_cm_well_mean=("gws_cm_well", "mean"),
        n_wells=("gws_cm_well", "count"),
    )
    .reset_index()
)
print(f"  Well-pixel-month combinations: {len(well_pixel):,}")

# Merge with GRACE GWS
matched = well_pixel.merge(
    gws,
    left_on=["grace_lat", "grace_lon", "year", "month"],
    right_on=["lat", "lon", "year", "month"],
    how="inner",
)

# Rename GRACE GWS column
matched["gws_cm_grace"] = matched["gws_cm"]
matched = matched.drop(columns=["gws_cm", "date"], errors="ignore")

print(f"  Matched well-GRACE pairs: {len(matched):,}")
unique_pixels = matched.groupby(["grace_lat", "grace_lon"]).ngroups
print(f"  Unique pixels matched: {unique_pixels}")

if len(matched) == 0:
    print("\n  ⚠  STILL NO MATCHES. Checking coordinate overlap:")
    print(f"  Well-pixel lats: {well_pixel['grace_lat'].min():.2f} to {well_pixel['grace_lat'].max():.2f}")
    print(f"  Well-pixel lons: {well_pixel['grace_lon'].min():.2f} to {well_pixel['grace_lon'].max():.2f}")
    print(f"  GRACE lats: {gws['lat'].min():.2f} to {gws['lat'].max():.2f}")
    print(f"  GRACE lons: {gws['lon'].min():.2f} to {gws['lon'].max():.2f}")
    # Sample a specific well-pixel entry and look for matching GRACE
    sample = well_pixel.iloc[0]
    print(f"\n  Sample well-pixel: lat={sample['grace_lat']:.4f}, lon={sample['grace_lon']:.4f}, "
          f"year={sample['year']}, month={sample['month']}")
    grace_match = gws[(gws['lat'] == sample['grace_lat']) & (gws['lon'] == sample['grace_lon'])]
    print(f"  Exact GRACE match at that lat/lon: {len(grace_match)} rows")
    # Check near match
    grace_near = gws[(np.abs(gws['lat'] - sample['grace_lat']) < 0.01) & 
                     (np.abs(gws['lon'] - sample['grace_lon']) < 0.01)]
    print(f"  Near GRACE match (±0.01°): {len(grace_near)} rows")
    raise SystemExit("Cannot proceed — coordinate mismatch unresolved.")

# Compute annual mean (across 4 quarters) per pixel per year
annual_pairs = (
    matched.groupby(["grace_lat", "grace_lon", "year", "zone"])
    .agg(
        gws_cm_grace=("gws_cm_grace", "mean"),
        gws_cm_well_mean=("gws_cm_well_mean", "mean"),
        n_wells=("n_wells", "mean"),
        n_quarters=("month", "nunique"),
    )
    .reset_index()
)

# Keep only years with all 4 quarters
annual_pairs = annual_pairs[annual_pairs["n_quarters"] == 4].copy()
print(f"  Annual pairs (4 quarters): {len(annual_pairs):,}")

# ---------------------------------------------------------------------------
# 7. Compute validation statistics
# ---------------------------------------------------------------------------
print("\n[6/7] Computing validation statistics ...")

def compute_metrics(x, y):
    """Compute Pearson r, R², RMSE for two aligned series."""
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_clean = x[mask]
    y_clean = y[mask]
    if len(x_clean) < 10:
        return {"n": len(x_clean), "r": np.nan, "R2": np.nan, "RMSE": np.nan, "p_value": np.nan}
    r, p = stats.pearsonr(x_clean, y_clean)
    rmse = np.sqrt(np.mean((x_clean - y_clean) ** 2))
    return {"n": len(x_clean), "r": r, "R2": r ** 2, "p_value": p, "RMSE": rmse}

results = {}
for zid in ["Z1", "Z2"]:
    zdata = annual_pairs[annual_pairs["zone"] == zid]
    results[zid] = compute_metrics(zdata["gws_cm_grace"], zdata["gws_cm_well_mean"])

# Z1 + Z2 combined
results["Z1+Z2"] = compute_metrics(
    annual_pairs["gws_cm_grace"], annual_pairs["gws_cm_well_mean"]
)

# Print summary
print("\n" + "=" * 70)
print("VALIDATION RESULTS")
print("=" * 70)
print(f"  {'Zone':<10} {'N pairs':>10} {'Pearson r':>12} {'R²':>10} {'p-value':>12} {'RMSE (cm)':>12}")
print(f"  {'─' * 70}")
for zid, m in results.items():
    print(f"  {zid:<10} {m['n']:>10} {m['r']:>+12.4f} {m['R2']:>10.4f} {m['p_value']:>12.2e} {m['RMSE']:>12.2f}")

# ---------------------------------------------------------------------------
# 8. Save outputs
# ---------------------------------------------------------------------------
print("\n[7/7] Saving outputs ...")

# Matched pairs
pairs_path = OUT_DIR / "validation_matched_pairs_igp.csv"
annual_pairs.to_csv(pairs_path, index=False)
print(f"  Saved: {pairs_path}")

# Summary statistics
summary_rows = []
for zid, m in results.items():
    summary_rows.append({
        "zone": zid,
        "n_pairs": m["n"],
        "pearson_r": m["r"],
        "R2": m["R2"],
        "p_value": m["p_value"],
        "RMSE_cm": m["RMSE"],
    })
df_summary = pd.DataFrame(summary_rows)
summary_path = OUT_DIR / "validation_summary_igp.csv"
df_summary.to_csv(summary_path, index=False)
print(f"  Saved: {summary_path}")

# ---------------------------------------------------------------------------
# 9. Publication scatter plot
# ---------------------------------------------------------------------------
print("\n  Generating scatter plot ...")

# Guard against empty data
if len(annual_pairs) == 0:
    print("  ⚠  No annual pairs to plot. Skipping figure.")
else:
    fig, ax = plt.subplots(figsize=(7, 6))

    colors = {"Z1": "#D62728", "Z2": "#1F77B4"}

    for zid in ["Z1", "Z2"]:
        zdata = annual_pairs[annual_pairs["zone"] == zid]
        if len(zdata) > 0:
            ax.scatter(
                zdata["gws_cm_grace"],
                zdata["gws_cm_well_mean"],
                c=colors[zid],
                alpha=0.6,
                edgecolors="none",
                s=40,
                label=f"{zid} (n={results[zid]['n']}, R²={results[zid]['R2']:.3f})",
            )

    # 1:1 line
    lims = [
        min(annual_pairs["gws_cm_grace"].min(), annual_pairs["gws_cm_well_mean"].min()),
        max(annual_pairs["gws_cm_grace"].max(), annual_pairs["gws_cm_well_mean"].max()),
    ]
    pad = (lims[1] - lims[0]) * 0.1
    lims = [lims[0] - pad, lims[1] + pad]

    ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5, label="1:1 line")

    # Combined trend line
    mask = ~np.isnan(annual_pairs["gws_cm_grace"]) & ~np.isnan(annual_pairs["gws_cm_well_mean"])
    if mask.sum() > 1:
        slope, intercept, r_val, p_val, std_err = stats.linregress(
            annual_pairs.loc[mask, "gws_cm_grace"],
            annual_pairs.loc[mask, "gws_cm_well_mean"],
        )
        x_fit = np.linspace(lims[0], lims[1], 100)
        ax.plot(x_fit, slope * x_fit + intercept, "k-", linewidth=1.0,
                label=f"Combined (R²={results['Z1+Z2']['R2']:.3f})")

    ax.set_xlabel("GRACE GWS anomaly (cm)", fontsize=11)
    ax.set_ylabel("CGWB well GWS anomaly (cm)", fontsize=11)
    ax.set_title("IGP Groundwater Validation: GRACE vs CGWB\n(Bore/Tube wells, Unconfined, Sy=0.10)",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color="gray", linewidth=0.5, alpha=0.5)
    ax.axvline(x=0, color="gray", linewidth=0.5, alpha=0.5)

    # RMSE annotation
    if not np.isnan(results['Z1+Z2']['RMSE']):
        rmse_text = f"RMSE = {results['Z1+Z2']['RMSE']:.1f} cm"
        ax.text(0.05, 0.95, rmse_text, transform=ax.transAxes, fontsize=9,
                verticalalignment="top", bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    plt.tight_layout()
    fig_path = FIG_DIR / "validation_igp_scatter.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {fig_path}")

# ---------------------------------------------------------------------------
# 10. Final summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("B5b COMPLETE")
print("=" * 70)
print(f"  Strict-criterion wells: {len(df):,} (Z1: {(df['zone']=='Z1').sum()}, Z2: {(df['zone']=='Z2').sum()})")
print(f"  Matched pixel-year pairs: {len(annual_pairs)}")
if len(annual_pairs) > 0:
    print(f"  Primary validation R² (Z1+Z2): {results['Z1+Z2']['R2']:.4f}")
print(f"\n  Target: R² = 0.50–0.70 (Asoka 2017 IGP range)")
print("=" * 70)