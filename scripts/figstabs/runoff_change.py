import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.25, rc={"lines.linewidth": 2})

# shps
shp_dep = gpd.read_file("data/raw/observed/EPSs/C_aporte.shp")
shp_SA = gpd.read_file("data/raw/shps/Sudamérica.shp").to_crs({"init": "epsg:4326"})
shp_lks = gpd.read_file("data/raw/shps/lago_titicaca_sideteva_puno.shp").to_crs({"init": "epsg:4326"})

access_q_fut = xr.open_dataset("data/processed/future/runoff/ACCESS/median_ensemble_runoff.nc")
hadgem_q_fut = xr.open_dataset("data/processed/future/runoff/HADGEM/median_ensemble_runoff.nc")
mpiesm_q_fut = xr.open_dataset("data/processed/future/runoff/MPIESM/median_ensemble_runoff.nc")
q_pre = xr.open_dataset("data/processed/present/PISCO/runoff/median_ensemble_runoff.nc")

q_fut = xr.concat([100*(access_q_fut-(q_pre+1))/(q_pre+1),
                   100*(hadgem_q_fut-(q_pre+1))/(q_pre+1),
                   100*(mpiesm_q_fut-(q_pre+1))/(q_pre+1)], dim="models")

q_change_percent = q_fut.median(dim="models")
q_change_std = q_fut.std(dim="models")
#q_std_median = (xr.concat([access_q_fut, hadgem_q_fut, mpiesm_q_fut], dim="models").std(dim="models"))/(xr.concat([access_q_fut, hadgem_q_fut, mpiesm_q_fut], dim="models").median(dim="models") + 1)

# runoff changes
fig, (ax0, ax1) = plt.subplots(1, 2, figsize = (6, 4), dpi = 250, sharey=True, sharex=True, gridspec_kw = {'wspace':0, 'hspace':0})

plot_b = q_change_percent.q.plot(ax = ax0, cmap = "bwr_r", add_colorbar=False, levels= 6, vmax=25, vmin=-25)
axin = inset_axes(ax0, width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = ax0.transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Cambio de Q (%)', labelpad=-30, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0)

shp_SA.geometry.boundary.plot(ax = ax0, edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = ax0, edgecolor = "black", linewidth = .5)

ax0.set_ylim(-18.5, 0.5)
ax0.set_xlim(-81.75, -68)
ax0.set_ylabel("")
ax0.set_xlabel("")
ax0.set_title("")
ax0.xaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax0.yaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax0.grid(True, linestyle='--', color = "black", alpha = 0.1)

plot_b = q_change_std.q.plot(ax = ax1, cmap = "viridis", add_colorbar=False, levels =  np.arange(0, 21, 2))
axin = inset_axes(ax1, width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = ax1.transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('SD de Cambio de Q (%)', labelpad=-25, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0)

shp_SA.geometry.boundary.plot(ax = ax1, edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = ax1, edgecolor = "black", linewidth = .5)

ax1.set_ylim(-18.5, 0.5)
ax1.set_xlim(-81.75, -68)
ax1.set_ylabel("")
ax1.set_xlabel("")
ax1.set_title("")
ax1.xaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax1.yaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax1.grid(True, linestyle='--', color = "black", alpha = 0.1)

plt.savefig('output/runoff_mean_sd_change_type_02.png', bbox_inches='tight',pad_inches = 0.01, dpi = 300)
plt.close()

"""
plot_b = q_std_median.q.plot(ax = ax2, cmap = "viridis", add_colorbar=False, vmax=1, vmin=0)
axin = inset_axes(ax2, width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = ax2.transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('CV de Q', labelpad=-27, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0)

shp_SA.geometry.boundary.plot(ax = ax2, edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = ax2, edgecolor = "black", linewidth = .5)

ax2.set_ylim(-18.5, 0.5)
ax2.set_xlim(-81.75, -68)
ax2.set_ylabel("")
ax2.set_xlabel("")
ax2.set_title("")
ax2.xaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax2.yaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax2.grid(True, linestyle='--', color = "black", alpha = 0.1)
"""