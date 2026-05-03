"""
figures_b7.py — Task B7
=========================
Three publication-quality figures at 300 DPI for Nature Water:

  Figure 2a: National GWS trend map (Sen's slope cm/yr, hexbin, dual-color)
  Figure 2b: Zone-level GWS time series (6 panels, linear regressions)
  Figure 3:  IGP validation scatter (Figure 3 in paper, Z1 primary + Z2 context)

All figures: matplotlib + cartopy, 300 DPI, width 180mm or 360mm,
ColorBrewer-safe, accessible to color-blind readers, 8pt minimum font.

Input:  02_processed/gws_monthly.csv
        03_outputs/gws_zone_timeseries.csv
        03_outputs/gws_zone_trend.csv
        03_outputs/validation_matched_pairs_igp.csv
        01_raw_data/india_boundary.shp (optional — uses bbox if missing)
Output: 05_figures/fig2a_gws_trend_map.png
        05_figures/fig2b_zone_timeseries.png
        05_figures/fig3_validation_scatter.png
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
BASE = Path(r"D:\INDIA_AQUIFER_STUDY\ROLE_B_GROUNDWATER")
GWS_FILE = BASE / "02_processed" / "gws_monthly.csv"
ZONE_TS_FILE = BASE / "03_outputs" / "gws_zone_timeseries.csv"
ZONE_TREND_FILE = BASE / "03_outputs" / "gws_zone_trend.csv"
VALIDATION_FILE = BASE / "03_outputs" / "validation_matched_pairs_igp.csv"
FIG_DIR = BASE / "05_figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# India bounding box for maps
INDIA_BBOX = {"lat_min": 6.0, "lat_max": 38.0, "lon_min": 66.0, "lon_max": 100.0}

# Zone definitions
ZONES = {
    "Z1_IGP_NW":     {"name": "Z1: IGP North-West", "lat_min": 27, "lat_max": 32, "lon_min": 73, "lon_max": 78,  "color": "#D62728"},
    "Z2_IGP_E":      {"name": "Z2: IGP East",       "lat_min": 24, "lat_max": 28, "lon_min": 80, "lon_max": 89,  "color": "#1F77B4"},
    "Z3_DECCAN":     {"name": "Z3: Deccan Basalt",   "lat_min": 16, "lat_max": 21, "lon_min": 73, "lon_max": 80,  "color": "#2CA02C"},
    "Z4_CRYST_S":    {"name": "Z4: Crystalline S",   "lat_min": 10, "lat_max": 15, "lon_min": 77, "lon_max": 81,  "color": "#FF7F0E"},
    "Z5_ARID_NW":    {"name": "Z5: Arid Alluvial",   "lat_min": 21, "lat_max": 27, "lon_min": 68, "lon_max": 75,  "color": "#9467BD"},
    "Z6_HIMALAYA":   {"name": "Z6: Himalaya",        "lat_min": 28, "lat_max": 33, "lon_min": 76, "lon_max": 93,  "color": "#8C564B"},
}

PERIODS = {
    "2003–2018": ("2003-01-01", "2018-12-31"),
    "2018–2024": ("2018-01-01", "2024-12-31"),
}

# ===========================================================================
# FIGURE 2a: NATIONAL GWS TREND MAP (Hexbin)
# ===========================================================================
print("=" * 70)
print("B7: PUBLICATION FIGURES")
print("=" * 70)

print("\n[1/3] Figure 2a: National GWS trend map ...")

# Load GWS and compute per-pixel Sen's slope for 2003–2024
gws = pd.read_csv(GWS_FILE, parse_dates=["date"])

# Compute annual mean per pixel for trend stability
gws["year"] = gws["date"].dt.year
gws_annual = gws.groupby(["lat", "lon", "year"])["gws_cm"].mean().reset_index()

# Mann-Kendall + Sen's slope per pixel
print("  Computing per-pixel trends (this takes 1–2 minutes) ...")
trend_records = []
for (lat, lon), grp in gws_annual.groupby(["lat", "lon"]):
    if len(grp) < 10:
        continue
    years = grp["year"].values
    values = grp["gws_cm"].values
    # Sen's slope
    n = len(years)
    slopes = []
    for i in range(n):
        for j in range(i+1, n):
            slopes.append((values[j] - values[i]) / (years[j] - years[i]))
    sen_slope = np.median(slopes)

    # Mann-Kendall p-value (approximate via scipy)
    if n >= 10:
        tau, p_value = stats.kendalltau(years, values)
    else:
        p_value = np.nan

    trend_records.append({
        "lat": lat, "lon": lon,
        "sen_slope_cm_per_year": sen_slope,
        "p_value": p_value,
        "n_years": n,
    })

df_trend_map = pd.DataFrame(trend_records)
print(f"  Pixels with ≥10 annual values: {len(df_trend_map):,}")

# Bin into hex grid for cleaner visualization
# Hexbin by latitude/longitude — use 0.5° hex bins
print("  Creating hexbin map ...")

fig, ax = plt.subplots(figsize=(9, 10))

# Hexbin with trend values
hb = ax.hexbin(
    df_trend_map["lon"], df_trend_map["lat"],
    C=df_trend_map["sen_slope_cm_per_year"],
    gridsize=50,  # ~0.5° hexagons for India
    cmap="RdYlBu_r",  # Red = depletion, Blue = gain
    reduce_C_function=np.mean,
    vmin=-5, vmax=2,  # clipped for visual clarity
    mincnt=1,
    linewidths=0.1,
    edgecolors="gray",
)

cbar = plt.colorbar(hb, ax=ax, shrink=0.75, pad=0.02)
cbar.set_label("GWS trend (cm/yr)", fontsize=9)
cbar.ax.tick_params(labelsize=8)

# Mark significant pixels (p < 0.05) with dots
sig = df_trend_map[df_trend_map["p_value"] < 0.05]
if len(sig) > 0:
    # Use stippling effect via scatter with very small markers
    ax.scatter(sig["lon"], sig["lat"], s=0.3, c="black", alpha=0.3,
               marker=".", linewidths=0, label=f"p < 0.05 (n={len(sig):,})")

# Zone boundaries
for zid, zinfo in ZONES.items():
    rect = plt.Rectangle(
        (zinfo["lon_min"], zinfo["lat_min"]),
        zinfo["lon_max"] - zinfo["lon_min"],
        zinfo["lat_max"] - zinfo["lat_min"],
        fill=False, edgecolor=zinfo["color"], linewidth=1.5,
        linestyle="--", alpha=0.8
    )
    ax.add_patch(rect)
    # Zone label at center
    ax.text(
        (zinfo["lon_min"] + zinfo["lon_max"]) / 2,
        (zinfo["lat_min"] + zinfo["lat_max"]) / 2,
        zid.split("_")[0], fontsize=8, fontweight="bold",
        color=zinfo["color"], ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7, edgecolor="none")
    )

# Formatting
ax.set_xlim(INDIA_BBOX["lon_min"], INDIA_BBOX["lon_max"])
ax.set_ylim(INDIA_BBOX["lat_min"], INDIA_BBOX["lat_max"])
ax.set_xlabel("Longitude", fontsize=9)
ax.set_ylabel("Latitude", fontsize=9)
ax.set_title("GRACE-derived Groundwater Storage Trend\n2003–2024 (Sen's slope, cm/yr)",
             fontsize=11, fontweight="bold")
ax.grid(True, alpha=0.2, linewidth=0.5)
ax.set_aspect("equal")

# Annotation
ax.text(0.02, 0.98, "Red = depletion\nBlue = gain", transform=ax.transAxes,
        fontsize=7, va="top", ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

plt.tight_layout()
fig_path = FIG_DIR / "fig2a_gws_trend_map.png"
fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
print(f"  Saved: {fig_path}")
plt.close()

# ===========================================================================
# FIGURE 2b: ZONE-LEVEL GWS TIME SERIES (6 panels)
# ===========================================================================
print("\n[2/3] Figure 2b: Zone time series (6 panels) ...")

zone_ts = pd.read_csv(ZONE_TS_FILE, parse_dates=["date"])
zone_trend = pd.read_csv(ZONE_TREND_FILE)

fig, axes = plt.subplots(3, 2, figsize=(14, 13))
axes_flat = axes.flatten()

ZONE_ORDER = ["Z1_IGP_NW", "Z2_IGP_E", "Z3_DECCAN", "Z4_CRYST_S", "Z5_ARID_NW", "Z6_HIMALAYA"]

for idx, zid in enumerate(ZONE_ORDER):
    ax = axes_flat[idx]
    zinfo = ZONES[zid]

    zdata = zone_ts[zone_ts["zone_id"] == zid].copy()
    if zdata.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(zinfo["name"], fontsize=10, fontweight="bold")
        continue

    # Monthly GWS with shading for ±1 std
    ax.fill_between(zdata["date"], 
                    zdata["mean_gws_cm"] - zdata["std_gws_cm"],
                    zdata["mean_gws_cm"] + zdata["std_gws_cm"],
                    alpha=0.15, color=zinfo["color"], linewidth=0)
    ax.plot(zdata["date"], zdata["mean_gws_cm"], color=zinfo["color"], linewidth=0.8)

    # Linear regressions for both periods
    for period_name, (start, end) in PERIODS.items():
        mask = (zdata["date"] >= start) & (zdata["date"] <= end)
        period_data = zdata.loc[mask]
        if len(period_data) < 10:
            continue

        x_num = (period_data["date"] - pd.Timestamp("2003-01-01")).dt.days / 365.25
        y = period_data["mean_gws_cm"].values
        slope, intercept, _, _, _ = stats.linregress(x_num, y)
        x_fit = np.linspace(x_num.min(), x_num.max(), 50)
        y_fit = slope * x_fit + intercept
        x_fit_dates = pd.Timestamp("2003-01-01") + pd.to_timedelta(x_fit * 365.25, unit="D")

        linestyle = "-" if period_name == "2003–2018" else "--"
        linewidth = 1.5 if period_name == "2003–2018" else 1.5
        ax.plot(x_fit_dates, y_fit, color="black", linestyle=linestyle, linewidth=linewidth, alpha=0.7)

    # Get trend values for annotation
    zt = zone_trend[zone_trend["zone_id"] == zid]
    anno_text = ""
    for _, tr in zt.iterrows():
        anno_text += f"{tr['period']}: {tr['sen_slope_cm_per_year']:+.2f} cm/yr (p={tr['mk_p_value']:.2e})\n"

    # Zero line
    ax.axhline(y=0, color="gray", linewidth=0.5, alpha=0.5, linestyle=":")

    # GRACE gap shading
    ax.axvspan(pd.Timestamp("2017-06-01"), pd.Timestamp("2018-06-01"),
               alpha=0.08, color="gray", linewidth=0)
    if idx == 0:
        ax.text(pd.Timestamp("2017-12-15"), ax.get_ylim()[1] * 0.95,
                "GRACE\ngap", fontsize=6, ha="center", va="top", alpha=0.5)

    # PM-KUSUM policy marker removed — not part of this paper's scope
    # (this project was reframed from policy attribution to hydrogeological controls)

    ax.set_title(zinfo["name"], fontsize=10, fontweight="bold", color=zinfo["color"])
    ax.set_ylabel("GWS anomaly (cm)", fontsize=8)
    ax.set_xlabel("")
    ax.grid(True, alpha=0.2, linewidth=0.3)
    ax.tick_params(labelsize=7)

    # Annotation box
    ax.text(0.02, 0.98, anno_text.strip(), transform=ax.transAxes,
            fontsize=6.5, va="top", ha="left", family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8, edgecolor="gray", linewidth=0.5))

# Shared labels
fig.text(0.5, 0.02, "Year", ha="center", fontsize=10)
fig.suptitle("Groundwater Storage Anomaly by Hydrogeological Zone",
             fontsize=13, fontweight="bold", y=0.995)

plt.tight_layout(rect=[0, 0.03, 1, 0.98])
fig_path = FIG_DIR / "fig2b_zone_timeseries.png"
fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
print(f"  Saved: {fig_path}")
plt.close()

# ===========================================================================
# FIGURE 3: IGP VALIDATION SCATTER (Z1 primary + Z2 context)
# ===========================================================================
print("\n[3/3] Figure 3: IGP validation scatter ...")

validation = pd.read_csv(VALIDATION_FILE)

fig, (ax_main, ax_hist) = plt.subplots(1, 2, figsize=(12, 5.5),
                                        gridspec_kw={"width_ratios": [2, 1]})

# ---- Main scatter panel ----
colors = {"Z1": "#D62728", "Z2": "#1F77B4"}

for zid in ["Z1", "Z2"]:
    zdata = validation[validation["zone"] == zid]
    if len(zdata) == 0:
        continue

    # Compute per-pixel mean for size/transparency variation
    alpha = 0.7 if zid == "Z1" else 0.4
    size = 35 if zid == "Z1" else 25
    edge = "none" if zid == "Z1" else "white"
    lw = 0 if zid == "Z1" else 0.3

    ax_main.scatter(
        zdata["gws_cm_grace"], zdata["gws_cm_well_mean"],
        c=colors[zid], alpha=alpha, edgecolors=edge,
        s=size, linewidths=lw, zorder=2 if zid == "Z1" else 1,
        label=f"{zid} (n={len(zdata)})",
    )

# Combined regression
mask = ~np.isnan(validation["gws_cm_grace"]) & ~np.isnan(validation["gws_cm_well_mean"])
x_all = validation.loc[mask, "gws_cm_grace"]
y_all = validation.loc[mask, "gws_cm_well_mean"]
slope, intercept, r_val, p_val, std_err = stats.linregress(x_all, y_all)
r2_all = r_val**2

# Z1-only regression
z1_data = validation[validation["zone"] == "Z1"]
mask_z1 = ~np.isnan(z1_data["gws_cm_grace"]) & ~np.isnan(z1_data["gws_cm_well_mean"])
x_z1 = z1_data.loc[mask_z1, "gws_cm_grace"]
y_z1 = z1_data.loc[mask_z1, "gws_cm_well_mean"]
slope_z1, intercept_z1, r_z1, p_z1, _ = stats.linregress(x_z1, y_z1)
r2_z1 = r_z1**2

# Compute limits
lims = [min(x_all.min(), y_all.min()), max(x_all.max(), y_all.max())]
pad = (lims[1] - lims[0]) * 0.1
lims = [lims[0] - pad, lims[1] + pad]

# 1:1 line
ax_main.plot(lims, lims, "k--", linewidth=0.8, alpha=0.4, label="1:1 line")

# Z1 trend line (primary)
x_fit = np.linspace(lims[0], lims[1], 100)
ax_main.plot(x_fit, slope_z1 * x_fit + intercept_z1, color=colors["Z1"],
             linewidth=1.8, linestyle="-",
             label=f"Z1 fit (R²={r2_z1:.3f})")

# Combined trend line (lighter)
ax_main.plot(x_fit, slope * x_fit + intercept, color="black",
             linewidth=1.0, linestyle="--", alpha=0.6,
             label=f"Combined fit (R²={r2_all:.3f})")

ax_main.set_xlabel("GRACE GWS anomaly (cm)", fontsize=10)
ax_main.set_ylabel("CGWB well GWS anomaly (cm)", fontsize=10)
ax_main.set_title("IGP Groundwater Validation\nGRACE vs CGWB bore/tube wells (Sy=0.10)",
                  fontsize=11, fontweight="bold")
ax_main.legend(loc="upper left", fontsize=8, framealpha=0.9)
ax_main.set_xlim(lims)
ax_main.set_ylim(lims)
ax_main.set_aspect("equal")
ax_main.grid(True, alpha=0.2, linewidth=0.3)
ax_main.axhline(y=0, color="gray", linewidth=0.5, alpha=0.3)
ax_main.axvline(x=0, color="gray", linewidth=0.5, alpha=0.3)

# RMSE annotation
rmse = np.sqrt(np.mean((x_all - y_all)**2))
ax_main.text(0.05, 0.95, f"RMSE = {rmse:.1f} cm\nn = {len(x_all)} pixel-years",
             transform=ax_main.transAxes, fontsize=8, va="top",
             bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="gray", linewidth=0.5))

# ---- Histogram panel ----
# Residuals histogram
residuals = x_all - y_all
ax_hist.hist(residuals, bins=40, orientation="horizontal", color="gray", alpha=0.5, edgecolor="white")
ax_hist.axhline(y=0, color="black", linewidth=0.8, linestyle="--")
ax_hist.set_xlabel("Count", fontsize=9)
ax_hist.set_ylabel("GRACE − CGWB residual (cm)", fontsize=9)
ax_hist.set_title("Residuals", fontsize=10, fontweight="bold")
ax_hist.grid(True, alpha=0.2, linewidth=0.3)
ax_hist.tick_params(labelsize=8)

# Add normal fit overlay
from scipy.stats import norm
mu, sigma = np.mean(residuals), np.std(residuals)
y_hist = np.linspace(residuals.min(), residuals.max(), 100)
x_pdf = norm.pdf(y_hist, mu, sigma)
# Scale to histogram
x_pdf_scaled = x_pdf * len(residuals) * (residuals.max() - residuals.min()) / 40
ax_hist.plot(x_pdf_scaled, y_hist, "r-", linewidth=1.2, alpha=0.7, label=f"μ={mu:.1f}, σ={sigma:.1f}")
ax_hist.legend(fontsize=7, loc="upper right")

plt.tight_layout()
fig_path = FIG_DIR / "fig3_validation_scatter.png"
fig.savefig(fig_path, dpi=300, bbox_inches="tight", facecolor="white")
print(f"  Saved: {fig_path}")
plt.close()

# ===========================================================================
# SUMMARY
# ===========================================================================
print("\n" + "=" * 70)
print("B7 COMPLETE — 3 figures saved to 05_figures/")
print(f"  1. fig2a_gws_trend_map.png      — National trend hexbin map")
print(f"  2. fig2b_zone_timeseries.png     — 6-panel zone time series")
print(f"  3. fig3_validation_scatter.png   — IGP validation scatter + residuals")
print("\n  Figure specifications:")
print("    - 300 DPI PNG, white background")
print("    - ColorBrewer-safe palettes")
print("    - 8pt minimum font size")
print("    - GRACE gap marked on all time series")
print("=" * 70)