import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

# ── Configuration ──────────────────────────────────────────────
GWS_FILE = r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER\02_processed\gws_monthly.csv"
OUT_DIR  = r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER\03_outputs"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
# ───────────────────────────────────────────────────────────────

print("Loading gws_monthly.csv...")
df = pd.read_csv(GWS_FILE, parse_dates=["date"])
print(f"  Loaded {len(df):,} rows")

# ── Step 1: Baseline anomaly correction (2004–2009 mean) ────────
print("\nStep 1: Baseline correction (2004–2009 mean per pixel)...")
baseline = df[df.date.dt.year.between(2004, 2009)]
pixel_mean = baseline.groupby(["lat", "lon"])["gws_cm"].mean().reset_index()
pixel_mean.rename(columns={"gws_cm": "baseline_mean"}, inplace=True)

df = df.merge(pixel_mean, on=["lat", "lon"], how="left")
df["gws_anom"] = df["gws_cm"] - df["baseline_mean"]
df.drop(columns=["baseline_mean"], inplace=True)

print(f"  GWS anomaly mean : {df.gws_anom.mean():.2f} cm  (should be near 0)")
print(f"  GWS anomaly range: {df.gws_anom.min():.2f} to {df.gws_anom.max():.2f} cm")

# ── Step 2: Mann-Kendall + Sen's Slope function ─────────────────
def mk_sens(values):
    """Returns (tau, p_value, sens_slope_cm_per_month)"""
    n = len(values)
    if n < 4:
        return np.nan, np.nan, np.nan
    # Mann-Kendall S statistic
    s = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            diff = values[j] - values[i]
            if diff > 0: s += 1
            elif diff < 0: s -= 1
    # Variance
    var_s = n * (n - 1) * (2 * n + 5) / 18
    # Z score
    if s > 0: z = (s - 1) / np.sqrt(var_s)
    elif s < 0: z = (s + 1) / np.sqrt(var_s)
    else: z = 0
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    tau = s / (n * (n - 1) / 2)
    # Sen's slope
    slopes = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            slopes.append((values[j] - values[i]) / (j - i))
    sens = np.median(slopes)
    return tau, p, sens

# ── Step 3: Per-pixel trend for two periods ─────────────────────
periods = {
    "2003_2018": (2003, 2018),
    "2018_2024": (2018, 2024),
}

all_trends = []

for period_name, (yr_start, yr_end) in periods.items():
    print(f"\nStep 3: Computing trends for {period_name}...")
    subset = df[df.date.dt.year.between(yr_start, yr_end)]
    pixels = subset.groupby(["lat", "lon"])

    results = []
    total = len(pixels)
    for i, ((lat, lon), grp) in enumerate(pixels):
        if i % 500 == 0:
            print(f"  {i}/{total} pixels...", end="\r")
        grp_sorted = grp.sort_values("date")
        vals = grp_sorted["gws_anom"].values
        tau, p, sens = mk_sens(vals)
        results.append({
            "lat": lat, "lon": lon,
            "period": period_name,
            "tau": tau,
            "p_value": p,
            "sens_slope_cm_month": sens,
            "sens_slope_cm_year": sens * 12,
            "significant": int(p < 0.05) if not np.isnan(p) else 0
        })

    period_df = pd.DataFrame(results)
    all_trends.append(period_df)
    sig = period_df["significant"].sum()
    print(f"  Done. {sig}/{len(period_df)} pixels significant (p<0.05)")
    print(f"  Median Sen's slope: {period_df.sens_slope_cm_year.median():.3f} cm/year")

# ── Step 4: Save gws_pixel_trend.csv ───────────────────────────
# NOTE: Output is named gws_PIXEL_trend.csv (per-pixel, not zone-level).
# Zone-level trends are produced by zone_summaries.py (gws_zone_trend.csv).
# These two files are distinct; do NOT rename this file to gws_zone_trend.csv.
print("\nStep 4: Saving outputs...")
trend_df = pd.concat(all_trends, ignore_index=True)
trend_path = Path(OUT_DIR) / "gws_pixel_trend.csv"
trend_df.to_csv(trend_path, index=False)
print(f"  Saved: {trend_path}")
print(f"  Rows : {len(trend_df):,}")

# ── Step 5: Acceleration (post-2018 vs pre-2018 slope) ─────────
print("\nStep 5: Computing acceleration...")
pre  = trend_df[trend_df.period == "2003_2018"][["lat","lon","sens_slope_cm_year"]].rename(
       columns={"sens_slope_cm_year": "slope_pre"})
post = trend_df[trend_df.period == "2018_2024"][["lat","lon","sens_slope_cm_year"]].rename(
       columns={"sens_slope_cm_year": "slope_post"})

accel = pre.merge(post, on=["lat","lon"])
accel["acceleration_cm_year"] = accel["slope_post"] - accel["slope_pre"]
accel["accelerating"] = (accel["acceleration_cm_year"] < 0).astype(int)

accel_path = Path(OUT_DIR) / "gws_pixel_acceleration.csv"
accel.to_csv(accel_path, index=False)
print(f"  Saved: {accel_path}")
print(f"  Pixels with accelerating depletion: {accel.accelerating.sum():,} ({100*accel.accelerating.mean():.1f}%)")
print(f"  Median acceleration: {accel.acceleration_cm_year.median():.3f} cm/year")

print("\nTask B4 complete.")
print("  Output: gws_pixel_trend.csv        — per-pixel Mann-Kendall trends (two periods)")
print("  Output: gws_pixel_acceleration.csv — per-pixel acceleration (post minus pre slope)")
print("  NOTE: Zone-level trends are in gws_zone_trend.csv produced by zone_summaries.py")