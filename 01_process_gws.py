import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────
GRACE_FILE = r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER\01_raw_data\grace_fo\CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc"
GLDAS_DIR  = r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER\01_raw_data\gldas"
OUT_DIR    = r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER\02_processed"

LAT_MIN, LAT_MAX = 6.0,  38.0
LON_MIN, LON_MAX = 66.0, 100.0

SM_VARS  = ["SoilMoi0_10cm_inst", "SoilMoi10_40cm_inst",
            "SoilMoi40_100cm_inst", "SoilMoi100_200cm_inst"]
CAN_VAR  = "CanopInt_inst"
SWE_VAR  = "SWE_inst"
# ───────────────────────────────────────────────────────────────

Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

# ── Step 1: Load and clip GRACE ─────────────────────────────────
print("Step 1: Loading GRACE-FO...")
grace_raw = xr.open_dataset(GRACE_FILE)

time_decoded = pd.to_datetime(
    grace_raw.time.values, unit="D", origin=pd.Timestamp("2002-01-01")
)
grace = grace_raw["lwe_thickness"].assign_coords(time=time_decoded)
grace_india = grace.sel(
    lat=slice(LAT_MIN, LAT_MAX),
    lon=slice(LON_MIN, LON_MAX)
)
print(f"  GRACE clipped: {grace_india.shape}")

# ── Step 2: Load GLDAS ──────────────────────────────────────────
print("\nStep 2: Loading GLDAS files...")
gldas_files = sorted(Path(GLDAS_DIR).glob("GLDAS_NOAH025_M.A*.nc4"))
print(f"  Found {len(gldas_files)} files")

ds_list = []
for f in gldas_files:
    ds = xr.open_dataset(f)[SM_VARS + [CAN_VAR, SWE_VAR]]
    ds_list.append(ds)

gldas = xr.concat(ds_list, dim="time")
gldas_india = gldas.sel(
    lat=slice(LAT_MIN, LAT_MAX),
    lon=slice(LON_MIN, LON_MAX)
)
print(f"  GLDAS clipped: {len(gldas_india.time)} months")

# ── Step 3: Compute GLDAS total ─────────────────────────────────
print("\nStep 3: Computing soil moisture + canopy + SWE...")

def clean(da):
    return da.where(da < 1e+14)

SM_total = sum(clean(gldas_india[v]) for v in SM_VARS) / 10.0
CAN      = clean(gldas_india[CAN_VAR]) / 10.0
SWE      = clean(gldas_india[SWE_VAR]).where(
               clean(gldas_india[SWE_VAR]) < 5000) / 10.0

GLDAS_total = SM_total + CAN + SWE

print(f"  SM   range (cm): {float(SM_total.min()):.1f} to {float(SM_total.max()):.1f}")
print(f"  CAN  range (cm): {float(CAN.min()):.1f} to {float(CAN.max()):.1f}")
print(f"  SWE  range (cm): {float(SWE.min()):.1f} to {float(SWE.max()):.1f}")
print(f"  GLDAS total (cm): {float(GLDAS_total.min()):.1f} to {float(GLDAS_total.max()):.1f}")

# ── Step 4: Align time and compute GWS ──────────────────────────
print("\nStep 4: Aligning GRACE and GLDAS...")

grace_times = pd.DatetimeIndex(grace_india.time.values).to_period("M").to_timestamp()
gldas_times = pd.DatetimeIndex(gldas_india.time.values).to_period("M").to_timestamp()

grace_india = grace_india.assign_coords(time=grace_times)
GLDAS_total = GLDAS_total.assign_coords(time=gldas_times)

# Deduplicate GRACE
_, unique_idx = np.unique(grace_india.time.values, return_index=True)
grace_india = grace_india.isel(time=unique_idx)
print(f"  GRACE after dedup: {len(grace_india.time)} months")

# Common months
common = sorted(set(pd.DatetimeIndex(grace_india.time.values)) &
                set(pd.DatetimeIndex(gldas_india.time.values)))
print(f"  Common months: {len(common)} ({common[0].strftime('%Y-%m')} to {common[-1].strftime('%Y-%m')})")

grace_common = grace_india.sel(time=common)
gldas_common = GLDAS_total.sel(time=common)

# Regrid GLDAS to GRACE grid
gldas_regrid = gldas_common.interp(
    lat=grace_common.lat,
    lon=grace_common.lon,
    method="linear"
)

GWS = grace_common - gldas_regrid
GWS.name = "gws_cm"
print("  GWS computed.")

# ── Step 5: Save CSV ─────────────────────────────────────────────
print("\nStep 5: Saving gws_monthly.csv...")

records = []
for t in GWS.time.values:
    da  = GWS.sel(time=t)
    df  = da.to_dataframe(name="gws_cm").reset_index()
    df  = df.dropna(subset=["gws_cm"])
    df["date"] = pd.Timestamp(t).strftime("%Y-%m-%d")
    records.append(df[["date", "lat", "lon", "gws_cm"]])

out_df = pd.concat(records, ignore_index=True)

# Clip physically implausible values (Himalayan glacier pixels)
before = len(out_df)
out_df = out_df[out_df.gws_cm.between(-200, 150)]
clipped = before - len(out_df)
print(f"  Clipped {clipped:,} extreme pixels ({100*clipped/before:.2f}%)")

out_path = Path(OUT_DIR) / "gws_monthly.csv"
out_df.to_csv(out_path, index=False)

print(f"  Saved: {out_path}")
print(f"  Rows : {len(out_df):,}")
print(f"  Sample:\n{out_df.head()}")
print("\nTask B3 complete.")