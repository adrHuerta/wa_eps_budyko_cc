import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray

import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

exec(open("src/discrete_cmap.py").read())
exec(open("src/crop_mask.py").read())

plt.rcParams["font.family"] = "Arial"
plt.rcParams['hatch.linewidth'] = 0.5
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.25, rc={"lines.linewidth": 2})

# shps
shp_dep = gpd.read_file("data/raw/observed/EPSs/C_aporte.shp")
shp_SA = gpd.read_file("data/raw/shps/Sudamérica.shp").to_crs({"init": "epsg:4326"})
shp_lks = gpd.read_file("data/raw/shps/lago_titicaca_sideteva_puno.shp").to_crs({"init": "epsg:4326"})

# gridded
pisco_p = xr.open_dataset("data/raw/observed/PISCO/prec/PISCOpd.nc").resample(z="1Y").sum(skipna = False).mean(dim="z")
pisco_p = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=pisco_p)
shp_exp_grid = xr_shp_to_grid(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_array=pisco_p.p, to_drop=["spatial_ref"])
pisco_p = xr_mask(grid_mask=shp_exp_grid, netcdf_i=pisco_p)

pisco_tx = xr.open_dataset("data/raw/observed/PISCO/temp/PISCOdtx_v1.1.nc").resample(time="1Y").mean().mean(dim="time")
pisco_tx = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=pisco_tx)
shp_exp_grid = xr_shp_to_grid(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_array=pisco_tx.tx, to_drop=["spatial_ref"])
pisco_tx = xr_mask(grid_mask=shp_exp_grid, netcdf_i=pisco_tx)

# point (after pre-processing)
pisco_runoff = gpd.read_file("data/processed/present/PISCO/runoff/Q_mov30yearly_shp.shp")
pisco_runoff_values = pd.read_csv("data/processed/present/PISCO/runoff/Q_mov30yearly.csv")
pisco_runoff["q_values_exp"] = pisco_runoff_values.mean(axis=0).values
pisco_runoff = pisco_runoff.sort_values(by="Region")

# prec and temp
fig, (ax0, ax1) = plt.subplots(1, 2, figsize = (6, 4), dpi = 250, sharey=True, sharex=True, gridspec_kw = {'wspace':0, 'hspace':0})

cmap = pl.cm.viridis
my_cmap = cmap(np.arange(cmap.N))
my_cmap[:,-1] = np.repeat(1, cmap.N)
my_cmap = ListedColormap(my_cmap)

plot_b = pisco_p.p.plot(ax = ax0, cmap=my_cmap, add_colorbar=False, levels= [0, 10, 50, 100, 200, 900])
axin = inset_axes(ax0, width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = ax0.transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Precipitación ($mm$)', labelpad=-30, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0)

shp_SA.geometry.boundary.plot(ax = ax0, edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = ax0, edgecolor = "black", linewidth = .5)
# shp_lks.plot(ax = ax0, edgecolor = "deepskyblue", color = "deepskyblue")

ax0.set_ylim(-18.5, 0.5)
ax0.set_xlim(-81.75, -68)
ax0.set_ylabel("")
ax0.set_xlabel("")
ax0.set_title("")
ax0.xaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax0.yaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax0.grid(True, linestyle='--', color = "black", alpha = 0.1)

cmap = pl.cm.Spectral_r
my_cmap = cmap(np.arange(cmap.N))
my_cmap[:,-1] = np.repeat(1, cmap.N)
my_cmap = ListedColormap(my_cmap)

plot_b = pisco_tx.tx.plot(ax = ax1, cmap = my_cmap, add_colorbar=False, levels= 6)
axin = inset_axes(ax1, width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = ax1.transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Temp. máxima ($^{\circ}C$)', labelpad=-26, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0)

shp_SA.geometry.boundary.plot(ax = ax1, edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = ax1, edgecolor = "black", linewidth = .5)
# shp_lks.plot(ax = ax1, edgecolor = "deepskyblue", color = "deepskyblue")

ax1.set_ylim(-18.5, 0.5)
ax1.set_xlim(-81.75, -68)
ax1.set_ylabel("")
ax1.set_xlabel("")
ax1.set_title("")
ax1.xaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax1.yaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax1.grid(True, linestyle='--', color = "black", alpha = 0.1)

plt.savefig('output/precp_temp_study_area.png', bbox_inches='tight',pad_inches = 0.01, dpi = 300)
plt.close()


# runoff
fig, (ax2, ax1) = plt.subplots(1, 2, figsize = (6, 4), dpi = 250, sharey=True, sharex=True, gridspec_kw = {'wspace':0, 'hspace':0})

cax = fig.add_axes([0.175, 0.15, 0.015, 0.275])

cmap_bar = pl.cm.viridis_r
my_cmap = cmap_bar(np.arange(cmap_bar.N))
my_cmap[:,-1] = np.repeat(1, cmap_bar.N)
my_cmap = ListedColormap(my_cmap)

#cmap_bar = mpl.cm.viridis_r
norm = mpl.colors.Normalize(vmin=0, vmax=1000)
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=my_cmap, norm=norm, orientation='vertical',extend='max')
cb1.set_label('Caudal ($mm$)',labelpad=-33, size = 6)
cb1.ax.tick_params(labelsize=5)

plot_b = pisco_runoff.plot(column = "q_values_exp", ax = ax2, cmap=my_cmap, legend=False, edgecolor='none', vmin=0, vmax=1000)
pisco_runoff[pisco_runoff["q_values_exp"].isna()].plot(hatch="/////////", facecolor="black", edgecolor='none', ax=ax2)
shp_SA.geometry.boundary.plot(ax = ax2, edgecolor = "black", linewidth = .5)

ax2.set_ylim(-18.5, 0.5)
ax2.set_xlim(-81.75, -68)
ax2.set_ylabel("")
ax2.set_xlabel("")
ax2.set_title("")
ax2.xaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax2.yaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax2.grid(True, linestyle='--', color = "black", alpha = 0.1)


cax = fig.add_axes([0.575, 0.15, 0.015, 0.275])
cmap_bar = discrete_cmap(N=14, base_cmap="gist_rainbow")
my_cmap = cmap_bar(np.arange(cmap_bar.N))
my_cmap[:,-1] = np.repeat(1, cmap_bar.N)
my_cmap = ListedColormap(my_cmap)

norm = mpl.colors.BoundaryNorm(np.arange(0,len(pisco_runoff["Region"].unique()) + .1), my_cmap.N)

plot_c = pisco_runoff.plot(column = "Region", ax = ax1, cmap=my_cmap,norm=norm, legend=False, edgecolor='none')
shp_SA.geometry.boundary.plot(ax = ax1, edgecolor = "black", linewidth = .65)
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap_bar, norm=norm, orientation='vertical')
cb1.set_label('Regiones',labelpad=-28.5, size = 6)
cb1.ax.tick_params(labelsize=5)
labels = pisco_runoff["Region"].unique() # np.arange(1, 14 + 1, 1)
loc = (np.arange(1, len(pisco_runoff["Region"].unique()) + 1, 1)-1) + .5
cb1.set_ticks(loc)
cb1.set_ticklabels(labels)

ax1.set_ylim(-18.5, 0.5)
ax1.set_xlim(-81.75, -68)
ax1.set_ylabel("")
ax1.set_xlabel("")
ax1.set_title("")
ax1.xaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax1.yaxis.set_tick_params(labelsize = 3.5, pad = -3)
ax1.grid(True, linestyle='--', color = "black", alpha = 0.1)

plt.savefig('output/runoff_study_area.png', bbox_inches='tight',pad_inches = 0.00, dpi = 300)
plt.close()