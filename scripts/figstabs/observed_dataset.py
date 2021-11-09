import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

plt.rcParams["font.family"] = "Arial"
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.25, rc={"lines.linewidth": 2})

# shps
shp_peru = gpd.read_file("data/shps/Sudamérica.shp").to_crs({"init": "epsg:4326"}).iloc[11:12]
shp_peru_no_lake = gpd.read_file("data/shps/Peru_no_lake.shp").to_crs({"init": "epsg:4326"})
shp_drainages = gpd.read_file("data/shps/vertientes.shp").to_crs({"init": "epsg:4326"})
shp_dep = gpd.read_file("data/shps/DEPARTAMENTOS.shp").to_crs({"init": "epsg:4326"})
shp_SA = gpd.read_file("data/shps/Sudamérica.shp").to_crs({"init": "epsg:4326"})
shp_lks = gpd.read_file("data/shps/lago_titicaca_sideteva_puno.shp").to_crs({"init": "epsg:4326"})

# gridded
pisco_p = xr.open_dataset("data/observed/PISCOpd.nc").resample(z="1Y").sum(skipna = False)
pisco_tx = xr.open_dataset("data/observed/PISCOdtx_v1.1.nc").resample(time="1Y").mean()


fig, (ax0, ax1) = plt.subplots(1, 2, figsize = (6, 4), dpi = 250, sharey=True, sharex=True, gridspec_kw = {'wspace':0, 'hspace':0})

plot_b = pisco_p.isel(z=2).p.plot(ax = ax0, cmap='viridis', add_colorbar=False, levels= [0, 10, 50, 100, 200, 900])
axin = inset_axes(ax0, width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = ax0.transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Precipitación ($mm$)', labelpad=-30, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0)

shp_SA.geometry.boundary.plot(ax = ax0, edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = ax0, edgecolor = "black", linewidth = .25)
shp_lks.plot(ax = ax0, edgecolor = "deepskyblue", color = "deepskyblue")

ax0.set_ylim(-18.5, 0.5)
ax0.set_xlim(-81.75, -68)
ax0.set_ylabel("")
ax0.set_xlabel("")
ax0.set_title("")
ax0.xaxis.set_tick_params(labelsize = 5.5, pad = -3)
ax0.yaxis.set_tick_params(labelsize = 5.5, pad = -3)
ax0.grid(True, linestyle='--', color = "black", alpha = 0.1)

plot_b = pisco_tx.isel(time=2).tx.plot(ax = ax1, cmap = "Spectral_r", add_colorbar=False, levels= 6)
axin = inset_axes(ax1, width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = ax1.transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Temp. máxima ($^{\circ}C$)', labelpad=-27, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0)

shp_SA.geometry.boundary.plot(ax = ax1, edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = ax1, edgecolor = "black", linewidth = .25)
shp_lks.plot(ax = ax1, edgecolor = "deepskyblue", color = "deepskyblue")

ax1.set_ylim(-18.5, 0.5)
ax1.set_xlim(-81.75, -68)
ax1.set_ylabel("")
ax1.set_xlabel("")
ax1.set_title("")
ax1.xaxis.set_tick_params(labelsize = 5.5, pad = -3)
ax1.yaxis.set_tick_params(labelsize = 5.5, pad = -3)
ax1.grid(True, linestyle='--', color = "black", alpha = 0.1)

plt.savefig('output/precp_temp_study_area.png', bbox_inches='tight',pad_inches = 0.01, dpi = 300)
plt.close()

# runoff
pisco_runoff = gpd.read_file("data/observed/runoff/PEQ_1981-2020_subcuencas_GR2M_Peru_SHP/Subbasins_GR2M_Peru.shp")
pisco_runoff_values = pd.read_excel("data/observed/runoff/PEQ_1981-2020_subcuencas_GR2M_Peru.xlsx", sheet_name="Q_mm", engine='openpyxl')
pisco_runoff_values_dates = pisco_runoff_values["Fecha"]
pisco_runoff_values = pisco_runoff_values.drop(["Fecha"], axis=1)
pisco_runoff_values = pisco_runoff_values.set_index(pisco_runoff_values_dates)
pisco_runoff_values = pisco_runoff_values.resample("1Y").apply(lambda x: np.sum(x))
pisco_runoff["q_values_exp"] = pisco_runoff_values.loc["1982-12-31"].values



fig, ax2 = plt.subplots(figsize = (4, 4), dpi = 250, sharey=True, sharex=True, gridspec_kw = {'wspace':.01, 'hspace':.01})

plot_b = pisco_runoff.plot(column = "q_values_exp", ax = ax2, cmap='viridis_r', legend=False, edgecolor='none', vmin=0, vmax=1000)

cax = fig.add_axes([0.3, 0.15, 0.03, 0.25])

cmap_bar = mpl.cm.viridis_r
norm = mpl.colors.Normalize(vmin=0, vmax=1000)
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap_bar, norm=norm, orientation='vertical',extend='max')
cb1.set_label('Caudal ($mm$)',labelpad=-37, size = 6)
cb1.ax.tick_params(labelsize=5)

shp_SA.geometry.boundary.plot(ax = ax2, edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = ax2, edgecolor = "black", linewidth = .25)
shp_lks.plot(ax = ax2, edgecolor = "deepskyblue", color = "deepskyblue")

ax2.set_ylim(-18.5, 0.5)
ax2.set_xlim(-81.75, -68)
ax2.set_ylabel("")
ax2.set_xlabel("")
ax2.set_title("")
ax2.xaxis.set_tick_params(labelsize = 5.5, pad = -3)
ax2.yaxis.set_tick_params(labelsize = 5.5, pad = -3)
ax2.grid(True, linestyle='--', color = "black", alpha = 0.1)

plt.savefig('output/runoff_study_area.png', bbox_inches='tight',pad_inches = 0.01, dpi = 300)
plt.close()