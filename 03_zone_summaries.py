"""
zone_summaries.py — Task B6
============================
Computes per-zone GWS time-series, trends, and acceleration.
Input:  02_processed/gws_monthly.csv
Output: 03_outputs/gws_zone_timeseries.csv
        03_outputs/gws_zone_trend.csv
        03_outputs/gws_zone_acceleration.csv

BLOCKS: Role D integration & regression
"""

import pandas as pd
import numpy as np
import pymannkendall as mk
from pathlib import Path

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
BASE = Path(r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER")
INPUT_FILE = BASE / "02_processed" / "gws_monthly.csv"
OUT_DIR = BASE / "03_outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Zone definitions (from handover section 10.1)
# ---------------------------------------------------------------------------
ZONES = {
    "Z1_IGP_NW": {
        "name": "IGP North-West",
        "lat_min": 27, "lat_max": 32,
        "lon_min": 73, "lon_max": 78,
    },
    "Z2_IGP_E": {
        "name": "IGP East",
        "lat_min": 24, "lat_max": 28,
        "lon_min": 80, "lon_max": 89,
    },
    "Z3_DECCAN": {
        "name": "Deccan Basalt",
        "lat_min": 16, "lat_max": 21,
        "lon_min": 73, "lon_max": 80,
    },
    "Z4_CRYST_S": {
        "name": "Crystalline South",
        "lat_min": 10, "lat_max": 15,
        "lon_min": 77, "lon_max": 81,
    },
    "Z5_ARID_NW": {
        "name": "Arid Alluvial NW",
        "lat_min": 21, "lat_max": 27,
        "lon_min": 68, "lon_max": 75,
    },
    "Z6_HIMALAYA": {
        "name": "Himalayan Foothills",
        "lat_min": 28, "lat_max": 33,
        "lon_min": 76, "lon_max": 93,
    },
}

PERIODS = {
    "2003-2018": ("2003-01-01", "2018-12-31"),
    "2018-2024": ("2018-01-01", "2024-12-31"),
    # The overlap at 2018 is intentional — it assigns 2018 to both periods
    # so the split reflects the policy launch year. Document this in methods.
}

# ---------------------------------------------------------------------------
# 2. Load data
# ---------------------------------------------------------------------------
print("Loading gws_monthly.csv ...")
df = pd.read_csv(INPUT_FILE, parse_dates=["date"])
print(f"  Rows: {len(df):,}")
print(f"  Date range: {df['date'].min()} to {df['date'].max()}")
print(f"  Columns: {list(df.columns)}")

# ---------------------------------------------------------------------------
# 3. Assign zone to each pixel row
# ---------------------------------------------------------------------------
print("\nAssigning zones ...")

def assign_zone(lat, lon):
    """Return zone_id if pixel falls in exactly one zone, else None."""
    hits = []
    for zid, z in ZONES.items():
        if (z["lat_min"] <= lat <= z["lat_max"]) and (z["lon_min"] <= lon <= z["lon_max"]):
            hits.append(zid)
    if len(hits) == 1:
        return hits[0]
    elif len(hits) > 1:
        # Overlap resolution (section 10.3): Z1 beats Z6 in overlap region
        if "Z1_IGP_NW" in hits:
            return "Z1_IGP_NW"
        # Any other overlap — warn and drop
        print(f"  WARNING: pixel ({lat:.4f}, {lon:.4f}) in multiple zones: {hits}")
        return None
    return None

df["zone_id"] = df.apply(lambda row: assign_zone(row["lat"], row["lon"]), axis=1)

# Report zone counts
zone_pixel_counts = df.groupby("zone_id")[["lat", "lon"]].apply(
    lambda g: g.drop_duplicates().shape[0]
)
print(f"  Pixels assigned per zone:")
for zid, count in zone_pixel_counts.items():
    print(f"    {zid}: {count:,} unique pixels")
print(f"  Unassigned pixels (outside all 6 zones): {df['zone_id'].isna().sum():,}")

# Drop unassigned
df = df.dropna(subset=["zone_id"]).copy()

# ---------------------------------------------------------------------------
# 4. Zone-level time-series (monthly mean & std)
# ---------------------------------------------------------------------------
print("\nComputing zone monthly summaries ...")

zone_ts = (
    df.groupby(["zone_id", "date"])["gws_cm"]
    .agg(["mean", "std", "count"])
    .reset_index()
)
zone_ts.columns = ["zone_id", "date", "mean_gws_cm", "std_gws_cm", "n_pixels"]

# Save
ts_path = OUT_DIR / "gws_zone_timeseries.csv"
zone_ts.to_csv(ts_path, index=False)
print(f"  Saved: {ts_path}")
print(f"  Rows: {len(zone_ts):,}")

# ---------------------------------------------------------------------------
# 5. Zone-level trends (Mann-Kendall + Sen's slope per period)
# ---------------------------------------------------------------------------
print("\nComputing per-zone, per-period trends ...")

trend_rows = []

for zid, zinfo in ZONES.items():
    zdata = zone_ts[zone_ts["zone_id"] == zid].copy()
    if zdata.empty:
        print(f"  {zid}: no data — skipping")
        continue

    for period_name, (start, end) in PERIODS.items():
        mask = (zdata["date"] >= start) & (zdata["date"] <= end)
        period_data = zdata.loc[mask, "mean_gws_cm"].dropna()

        if len(period_data) < 10:
            print(f"  {zid} {period_name}: insufficient data ({len(period_data)} months) — skipping")
            continue

        # Mann-Kendall + Sen's slope
        result = mk.original_test(period_data)
        sen_slope = mk.sens_slope(period_data)

        n_pixels_mean = zdata.loc[mask, "n_pixels"].mean()

        trend_rows.append({
            "zone_id": zid,
            "zone_name": zinfo["name"],
            "period": period_name,
            "sen_slope_cm_per_year": sen_slope.slope * 12,  # monthly → annual
            "mk_p_value": result.p,
            "mk_trend": result.trend,
            "n_months": len(period_data),
            "n_pixels_mean": int(n_pixels_mean),
        })

df_trends = pd.DataFrame(trend_rows)
trend_path = OUT_DIR / "gws_zone_trend.csv"
df_trends.to_csv(trend_path, index=False)
print(f"  Saved: {trend_path}\n")
print(df_trends.to_string(index=False))

# ---------------------------------------------------------------------------
# 6. Acceleration (post minus pre slope)
# ---------------------------------------------------------------------------
print("\n\nComputing acceleration ...")

accel_rows = []
for zid in df_trends["zone_id"].unique():
    pre_row = df_trends[(df_trends["zone_id"] == zid) & (df_trends["period"] == "2003-2018")]
    post_row = df_trends[(df_trends["zone_id"] == zid) & (df_trends["period"] == "2018-2024")]

    if pre_row.empty or post_row.empty:
        print(f"  {zid}: missing one period — skipping acceleration")
        continue

    slope_pre = pre_row["sen_slope_cm_per_year"].values[0]
    slope_post = post_row["sen_slope_cm_per_year"].values[0]
    accel = slope_post - slope_pre  # positive = less depletion / recovery

    accel_rows.append({
        "zone_id": zid,
        "zone_name": pre_row["zone_name"].values[0],
        "slope_2003_2018": slope_pre,
        "slope_2018_2024": slope_post,
        "acceleration_cm_per_year": accel,
        "interpretation": (
            "depletion accelerated" if (slope_pre < 0 and accel < 0) else
            "depletion slowed / reversed" if (slope_pre < 0 and accel > 0) else
            "check manually — unexpected pattern"
        ),
    })

df_accel = pd.DataFrame(accel_rows)
accel_path = OUT_DIR / "gws_zone_acceleration.csv"
df_accel.to_csv(accel_path, index=False)
print(f"  Saved: {accel_path}\n")
print(df_accel.to_string(index=False))

# ---------------------------------------------------------------------------
# 7. Summary print
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("B6 COMPLETE — 3 files saved to 03_outputs/")
print(f"  1. gws_zone_timeseries.csv   — {len(zone_ts):,} rows")
print(f"  2. gws_zone_trend.csv         — {len(df_trends)} rows")
print(f"  3. gws_zone_acceleration.csv  — {len(df_accel)} rows")
print("=" * 70)
print("\nNext: hand acceleration.csv to Role D, or write B5b validation script.")