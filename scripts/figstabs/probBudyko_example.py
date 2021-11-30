import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.25, rc={"lines.linewidth": 2})

# shps
shp_dep = gpd.read_file("data/raw/observed/EPSs/C_aporte.shp")
shp_SA = gpd.read_file("data/raw/shps/Sudamérica.shp").to_crs({"init": "epsg:4326"})
shp_lks = gpd.read_file("data/raw/shps/lago_titicaca_sideteva_puno.shp").to_crs({"init": "epsg:4326"})

omega_value = pd.read_csv("data/processed/others/omega/csv_omega_values.csv", index_col=0)
q_pre_median = xr.open_dataset("data/processed/present/PISCO/runoff/median_ensemble_runoff.nc")
q_pre_sd = xr.open_dataset("data/processed/present/PISCO/runoff/sd_ensemble_runoff.nc")
Q_shp = gpd.read_file("data/processed/present/PISCO/runoff/Q_mov15yearly_shp.shp")
cluster_areas = xr.open_dataset("data/processed/others/budyko_groups.nc")

# hist omega by region
fig, axes = plt.subplots(4, 3, figsize=(8, 40), dpi = 150)
i = 0
for triaxis in axes:
    for axis in triaxis:
        axis.hist(omega_value[omega_value.columns[i]].values.tolist(), bins = 10)
        axis.set_title(omega_value.columns[i])
        axis.xaxis.set_tick_params(labelsize = 7, pad = -3)
        axis.yaxis.set_tick_params(labelsize = 7, pad = -3)
        axis.set_xlabel("")
        axis.set_ylabel("")
        i = i+1
        if i == 11:
            break
plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.15,
                    hspace=0.4)
axes.flat[-1].set_visible(False)
fig.text(0.5, 0.05, 'Omega (adimensional)', ha='center', size = 9)
fig.text(0.05, 0.5, 'Frecuencia', va='center', rotation='vertical', size = 9)
plt.close()

# hist omega for "M" region
fig, ax = plt.subplots(figsize=(7, 4), dpi = 150)
ax.hist(omega_value["M"].values.tolist(), bins = 10)
ax.xaxis.set_tick_params(labelsize = 10, pad = -3)
ax.yaxis.set_tick_params(labelsize = 10, pad = -3)
ax.set_xlabel("Omega (adimensional)")
ax.set_ylabel("Frecuencia")

plt.savefig('output/omega_distribution_region_M.png', bbox_inches='tight',pad_inches = 0.01, dpi = 150)
plt.close()


# example of q for "M" region
for a, b in zip([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1, np.nan]):
    cluster_areas.gridded_groups.values[cluster_areas.gridded_groups.values == a] = b

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(8, 5), dpi = 150, sharey=True, gridspec_kw = {'wspace':0.1})
(cluster_areas.gridded_groups * q_pre_median.q).plot(ax=ax0)
shp_SA.geometry.boundary.plot(ax = ax0, edgecolor = "black", linewidth = .65)
ax0.set_ylim(-17.9, -12)
ax0.set_xlim(-72.5, -68.5)
ax0.set_ylabel("")
ax0.set_xlabel("")
ax0.set_title("")
ax0.xaxis.set_tick_params(labelsize = 7.5, pad = -3)
ax0.yaxis.set_tick_params(labelsize = 7.5, pad = -3)
ax0.grid(True, linestyle='--', color = "black", alpha = 0.1)

(cluster_areas.gridded_groups * q_pre_sd.q).plot(ax=ax1)
shp_SA.geometry.boundary.plot(ax = ax1, edgecolor = "black", linewidth = .65)
ax1.set_ylim(-17.9, -12)
ax1.set_xlim(-72.5, -68.5)
ax1.set_ylabel("")
ax1.set_xlabel("")
ax1.set_title("")
ax1.xaxis.set_tick_params(labelsize = 7.5, pad = -3)
ax1.yaxis.set_tick_params(labelsize = 7.5, pad = -3)
ax1.grid(True, linestyle='--', color = "black", alpha = 0.1)

plt.savefig('output/Q_median_cv_present_region_M.png', bbox_inches='tight',pad_inches = 0.01, dpi = 150)
plt.close()