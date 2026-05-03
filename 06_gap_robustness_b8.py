"""
gap_robustness_b8.py — Task B8
================================
Tests sensitivity of pre/post-2018 trend split to the choice of cutoff date.
Compares January 2018 vs July 2018 as the boundary between periods.
If zone-level acceleration signs flip across cutoffs, the GRACE/GRACE-FO
transition gap (Jun 2017 – Jun 2018) is contaminating results.

Input:  03_outputs/gws_zone_timeseries.csv (from B6)
Output: console comparison table
        (optionally) 03_outputs/gap_sensitivity.csv
"""

import pandas as pd
import numpy as np
import pymannkendall as mk
from pathlib import Path

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
BASE = Path(r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER")
INPUT_FILE = BASE / "03_outputs" / "gws_zone_timeseries.csv"
OUT_DIR = BASE / "03_outputs"

# ---------------------------------------------------------------------------
# 1. Load zone time-series
# ---------------------------------------------------------------------------
print("Loading zone time-series ...")
zone_ts = pd.read_csv(INPUT_FILE, parse_dates=["date"])
print(f"  Rows: {len(zone_ts):,}")
print(f"  Zones: {zone_ts['zone_id'].nunique()}")
print(f"  Date range: {zone_ts['date'].min()} to {zone_ts['date'].max()}")

# ---------------------------------------------------------------------------
# 2. Define cutoff dates to test
# ---------------------------------------------------------------------------
# Jan 2018: pre = 2003-01 to 2017-12, post = 2018-01 to 2024-12
# Jul 2018: pre = 2003-01 to 2018-06, post = 2018-07 to 2024-12
# The mission gap was Jun 2017 – Jun 2018. Jul 2018 cutoff places the
# entire gap in the PRE period, giving the cleanest post-period signal.

CUTOFFS = {
    "Jan 2018": {
        "pre_end": "2017-12-31",
        "post_start": "2018-01-01",
    },
    "Jul 2018": {
        "pre_end": "2018-06-30",
        "post_start": "2018-07-01",
    },
}

PRE_START = "2003-01-01"
POST_END = "2024-12-31"

# ---------------------------------------------------------------------------
# 3. Function to compute trend for a given date range
# ---------------------------------------------------------------------------
def compute_trend(series, period_name):
    """Run Mann-Kendall + Sen's slope on a pandas Series."""
    clean = series.dropna()
    if len(clean) < 10:
        return {"sen_slope_cm_per_year": np.nan, "mk_p_value": np.nan,
                "mk_trend": "insufficient_data", "n_months": len(clean)}

    result = mk.original_test(clean)
    sen_slope = mk.sens_slope(clean)

    return {
        "sen_slope_cm_per_year": sen_slope.slope * 12,  # monthly → annual
        "mk_p_value": result.p,
        "mk_trend": result.trend,
        "n_months": len(clean),
    }

# ---------------------------------------------------------------------------
# 4. Run comparison
# ---------------------------------------------------------------------------
print("\n" + "=" * 90)
print("GRACE/GRACE-FO TRANSITION GAP ROBUSTNESS TEST (B8)")
print("=" * 90)

all_rows = []

for zone_id in sorted(zone_ts["zone_id"].unique()):
    zdata = zone_ts[zone_ts["zone_id"] == zone_id].set_index("date")["mean_gws_cm"]

    print(f"\n{'─' * 60}")
    print(f"  Zone: {zone_id}")
    print(f"{'─' * 60}")
    print(f"  {'Cutoff':<12} {'Period':<20} {'Sen slope':>12} {'p-value':>14} {'Trend':>16} {'N months':>10}")
    print(f"  {'─' * 60}")

    for cutoff_name, bounds in CUTOFFS.items():
        pre_data = zdata.loc[PRE_START:bounds["pre_end"]]
        post_data = zdata.loc[bounds["post_start"]:POST_END]

        pre_trend = compute_trend(pre_data, "pre")
        post_trend = compute_trend(post_data, "post")

        accel = (post_trend["sen_slope_cm_per_year"] - pre_trend["sen_slope_cm_per_year"]
                 if not np.isnan(pre_trend["sen_slope_cm_per_year"])
                 and not np.isnan(post_trend["sen_slope_cm_per_year"])
                 else np.nan)

        print(f"  {cutoff_name:<12} {'pre (2003–{})'.format(bounds['pre_end'][:4]):<20} "
              f"{pre_trend['sen_slope_cm_per_year']:>+12.4f} "
              f"{pre_trend['mk_p_value']:>14.2e} "
              f"{pre_trend['mk_trend']:>16} "
              f"{pre_trend['n_months']:>10}")

        print(f"  {cutoff_name:<12} {'post ({}–2024)'.format(bounds['post_start'][:4]):<20} "
              f"{post_trend['sen_slope_cm_per_year']:>+12.4f} "
              f"{post_trend['mk_p_value']:>14.2e} "
              f"{post_trend['mk_trend']:>16} "
              f"{post_trend['n_months']:>10}")

        if not np.isnan(accel):
            direction = "→ FASTER DEPLETION" if accel < 0 else "→ SLOWER DEPLETION / GAIN"
            print(f"  {'':12} {'ACCELERATION':<20} {accel:>+12.4f} {'':>14} {'':>16} {'':>10}  {direction}")

        all_rows.append({
            "zone_id": zone_id,
            "cutoff": cutoff_name,
            "period": "pre",
            "sen_slope_cm_per_year": pre_trend["sen_slope_cm_per_year"],
            "mk_p_value": pre_trend["mk_p_value"],
            "mk_trend": pre_trend["mk_trend"],
            "n_months": pre_trend["n_months"],
        })
        all_rows.append({
            "zone_id": zone_id,
            "cutoff": cutoff_name,
            "period": "post",
            "sen_slope_cm_per_year": post_trend["sen_slope_cm_per_year"],
            "mk_p_value": post_trend["mk_p_value"],
            "mk_trend": post_trend["mk_trend"],
            "n_months": post_trend["n_months"],
        })

# ---------------------------------------------------------------------------
# 5. Stability check — compare accelerations across cutoffs
# ---------------------------------------------------------------------------
print("\n")
print("=" * 90)
print("STABILITY CHECK — Does acceleration sign flip across cutoffs?")
print("=" * 90)

df_all = pd.DataFrame(all_rows)
df_all["acceleration"] = np.nan

# Compute acceleration per zone-cutoff
for (zone_id, cutoff), grp in df_all.groupby(["zone_id", "cutoff"]):
    pre_val = grp[grp["period"] == "pre"]["sen_slope_cm_per_year"].values
    post_val = grp[grp["period"] == "post"]["sen_slope_cm_per_year"].values
    if len(pre_val) == 1 and len(post_val) == 1:
        accel = post_val[0] - pre_val[0]
        df_all.loc[(df_all["zone_id"] == zone_id) & (df_all["cutoff"] == cutoff), "acceleration"] = accel

# Pivot to compare
pivot = df_all.pivot_table(
    index="zone_id",
    columns="cutoff",
    values="acceleration",
    aggfunc="first"
)

print(f"\n{'Zone':<20} {'Jan 2018 accel':>16} {'Jul 2018 accel':>16} {'Difference':>14} {'FLIP?':>8}")
print(f"{'─' * 74}")
for zone_id in pivot.index:
    jan = pivot.loc[zone_id, "Jan 2018"]
    jul = pivot.loc[zone_id, "Jul 2018"]
    diff = jul - jan
    flip = "YES ⚠" if (jan * jul < 0) else "no"
    sign = " ⚠ FLIP DETECTED" if jan * jul < 0 else ""
    print(f"{zone_id:<20} {jan:>+16.4f} {jul:>+16.4f} {diff:>+14.4f}  {flip}")

# ---------------------------------------------------------------------------
# 6. Save full sensitivity table
# ---------------------------------------------------------------------------
sensitivity_path = OUT_DIR / "gap_sensitivity.csv"
df_all.to_csv(sensitivity_path, index=False)
print(f"\nFull sensitivity table saved to: {sensitivity_path}")

# ---------------------------------------------------------------------------
# 7. Verdict
# ---------------------------------------------------------------------------
n_flips = 0
for zone_id in pivot.index:
    jan = pivot.loc[zone_id, "Jan 2018"]
    jul = pivot.loc[zone_id, "Jul 2018"]
    if jan * jul < 0:
        n_flips += 1

print("\n" + "=" * 90)
print("VERDICT")
print("=" * 90)
if n_flips == 0:
    print(f"  ✅  No acceleration sign flips detected.")
    print(f"  The post-2018 spatial pattern is robust to cutoff choice.")
    print(f"  Use Jul 2018 cutoff in the paper (places entire gap in pre-period).")
elif n_flips == 6:
    print(f"  ❌  ALL zones flip — the gap is likely dominating the signal.")
    print(f"  Flag this in methods. Consider restricting post-period to 2019–2024.")
else:
    print(f"  ⚠   {n_flips}/6 zones show sign flips.")
    print(f"  Those zones are unreliable. Report only stable zones in main text;")
    print(f"  flag unstable zones in supplementary.")

print("\nB8 COMPLETE.\n")