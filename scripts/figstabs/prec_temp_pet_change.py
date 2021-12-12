import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import rioxarray

exec(open("src/crop_mask.py").read())
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.25, rc={"lines.linewidth": 2})

# shps
shp_dep = gpd.read_file("data/raw/observed/EPSs/C_aporte.shp")
shp_SA = gpd.read_file("data/raw/shps/Sudamérica.shp").to_crs({"init": "epsg:4326"})
shp_lks = gpd.read_file("data/raw/shps/lago_titicaca_sideteva_puno.shp").to_crs({"init": "epsg:4326"})

# prec
access_p_fut = xr.open_dataset("data/processed/future/prec/p_SENAMHI_pp_d12k_ACCESS1-0_rcp85_scal.nc")
hadgem_p_fut = xr.open_dataset("data/processed/future/prec/p_SENAMHI_pp_d12k_HadGEM2-ES_rcp85_scal.nc")
mpiesm_p_fut = xr.open_dataset("data/processed/future/prec/p_SENAMHI_pp_d12k_MPI-ESM-LR_rcp85_scal.nc")
p_pre = xr.open_dataset("data/processed/present/PISCO/prec/P_mean_anual.nc")

p_fut = xr.concat([100*(access_p_fut-(p_pre+1))/(p_pre+1),
                   100*(hadgem_p_fut-(p_pre+1))/(p_pre+1),
                   100*(mpiesm_p_fut-(p_pre+1))/(p_pre+1)], dim="models")

p_change_percent = p_fut.median(dim="models")
p_change_std = p_fut.std(dim="models")

p_change_percent = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=p_change_percent)
shp_exp_grid = xr_shp_to_grid(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_array=p_change_percent.p, to_drop=["spatial_ref"])
p_change_percent = xr_mask(grid_mask=shp_exp_grid, netcdf_i=p_change_percent)
p_change_std = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=p_change_std)
p_change_std = xr_mask(grid_mask=shp_exp_grid, netcdf_i=p_change_std)

# tx
access_tx_fut = xr.open_dataset("data/processed/future/temp/tx_SENAMHI_tmax_d12k_ACCESS1-0_rcp85_scal.nc")
hadgem_tx_fut = xr.open_dataset("data/processed/future/temp/tx_SENAMHI_tmax_d12k_HadGEM2-ES_rcp85_scal.nc")
mpiesm_tx_fut = xr.open_dataset("data/processed/future/temp/tx_SENAMHI_tmax_d12k_MPI-ESM-LR_rcp85_scal.nc")
tx_pre = xr.open_dataset("data/processed/present/PISCO/temp/TX_mean_anual.nc")

tx_fut = xr.concat([100*(access_tx_fut-(tx_pre+1))/(tx_pre+1),
                   100*(hadgem_tx_fut-(tx_pre+1))/(tx_pre+1),
                   100*(mpiesm_tx_fut-(tx_pre+1))/(tx_pre+1)], dim="models")

tx_change_percent = tx_fut.median(dim="models")
tx_change_std = tx_fut.std(dim="models")

tx_change_percent = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=tx_change_percent)
shp_exp_grid = xr_shp_to_grid(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_array=tx_change_percent.tx, to_drop=["spatial_ref"])
tx_change_percent = xr_mask(grid_mask=shp_exp_grid, netcdf_i=tx_change_percent)
tx_change_std = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=tx_change_std)
tx_change_std = xr_mask(grid_mask=shp_exp_grid, netcdf_i=tx_change_std)

# tn
access_tn_fut = xr.open_dataset("data/processed/future/temp/tn_SENAMHI_tmin_d12k_ACCESS1-0_rcp85_scal.nc")
hadgem_tn_fut = xr.open_dataset("data/processed/future/temp/tn_SENAMHI_tmin_d12k_HadGEM2-ES_rcp85_scal.nc")
mpiesm_tn_fut = xr.open_dataset("data/processed/future/temp/tn_SENAMHI_tmin_d12k_MPI-ESM-LR_rcp85_scal.nc")
tn_pre = xr.open_dataset("data/processed/present/PISCO/temp/TN_mean_anual.nc")

min_negative_value = np.abs(np.min([access_tn_fut.tn.min(), hadgem_tn_fut.tn.min(), mpiesm_tn_fut.tn.min(), tn_pre.tn.min()])) + 1

access_tn_fut = access_tn_fut + min_negative_value
hadgem_tn_fut = hadgem_tn_fut + min_negative_value
mpiesm_tn_fut = mpiesm_tn_fut + min_negative_value
tn_pre = tn_pre + min_negative_value

tn_fut = xr.concat([100*(access_tn_fut-(tn_pre+1))/(tn_pre+1),
                   100*(hadgem_tn_fut-(tn_pre+1))/(tn_pre+1),
                   100*(mpiesm_tn_fut-(tn_pre+1))/(tn_pre+1)], dim="models")

tn_change_percent = tn_fut.median(dim="models")
tn_change_std = tn_fut.std(dim="models")

tn_change_percent = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=tn_change_percent)
shp_exp_grid = xr_shp_to_grid(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_array=tn_change_percent.tn, to_drop=["spatial_ref"])
tn_change_percent = xr_mask(grid_mask=shp_exp_grid, netcdf_i=tn_change_percent)
tn_change_std = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=tn_change_std)
tn_change_std = xr_mask(grid_mask=shp_exp_grid, netcdf_i=tn_change_std)

# pet
access_pet_fut = xr.open_dataset("data/processed/future/pet/pet_SENAMHI_tmax_d12k_ACCESS1-0_rcp85_scal.nc")
hadgem_pet_fut = xr.open_dataset("data/processed/future/pet/pet_SENAMHI_tmax_d12k_HadGEM2-ES_rcp85_scal.nc")
mpiesm_pet_fut = xr.open_dataset("data/processed/future/pet/pet_SENAMHI_tmax_d12k_MPI-ESM-LR_rcp85_scal.nc")
pet_pre = xr.open_dataset("data/processed/present/PISCO/pet/PET_mean_anual.nc")

pet_fut = xr.concat([100*(access_pet_fut-(pet_pre+1))/(pet_pre+1),
                   100*(hadgem_pet_fut-(pet_pre+1))/(pet_pre+1),
                   100*(mpiesm_pet_fut-(pet_pre+1))/(pet_pre+1)], dim="models")

pet_change_percent = pet_fut.median(dim="models")
pet_change_std = pet_fut.std(dim="models")

pet_change_percent = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=pet_change_percent)
shp_exp_grid = xr_shp_to_grid(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_array=pet_change_percent.pet, to_drop=["spatial_ref"])
pet_change_percent = xr_mask(grid_mask=shp_exp_grid, netcdf_i=pet_change_percent)
pet_change_std = xr_crop(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=pet_change_std)
pet_change_std = xr_mask(grid_mask=shp_exp_grid, netcdf_i=pet_change_std)


# plot
fig, axs = plt.subplots(2, 4, figsize = (5, 3.5), dpi = 250, sharey=True, sharex=True, gridspec_kw = {'wspace':0, 'hspace':0})

# prec
plot_b = p_change_percent.p.plot(ax = axs[0, 0], cmap = "bwr_r", add_colorbar=False, levels= 6, vmax=25, vmin=-25)
axin = inset_axes(axs[0, 0], width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = axs[0, 0].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Cambio de P (%)', labelpad=-16.5, size = 3.5)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1.5)

shp_SA.geometry.boundary.plot(ax = axs[0, 0], edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = axs[0, 0], edgecolor = "black", linewidth = .25)

axs[0, 0].set_ylim(-18.5, 0.5)
axs[0, 0].set_xlim(-81.75, -68)
axs[0, 0].set_ylabel("")
axs[0, 0].set_xlabel("")
axs[0, 0].set_title("")
axs[0, 0].xaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[0, 0].yaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[0, 0].grid(True, linestyle='--', color = "black", alpha = 0.1)

plot_b = p_change_std.p.plot(ax = axs[1, 0], cmap = "viridis", add_colorbar=False, levels =  np.arange(0, 21, 2))
axin = inset_axes(axs[1, 0], width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = axs[1, 0].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('SD de Cambio de P (%)', labelpad=-13, size = 3.5)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1.5)

shp_SA.geometry.boundary.plot(ax = axs[1, 0], edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = axs[1, 0], edgecolor = "black", linewidth = .25)

axs[1, 0].set_ylim(-18.5, 0.5)
axs[1, 0].set_xlim(-81.75, -68)
axs[1, 0].set_ylabel("")
axs[1, 0].set_xlabel("")
axs[1, 0].set_title("")
axs[1, 0].xaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[1, 0].yaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[1, 0].grid(True, linestyle='--', color = "black", alpha = 0.1)

# tx
plot_b = tx_change_percent.tx.plot(ax = axs[0, 1], cmap = "bwr", add_colorbar=False, levels= 6, vmax=15, vmin=-15)
axin = inset_axes(axs[0, 1], width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = axs[0, 1].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Cambio de Tx (%)', labelpad=-16.5, size = 3.5)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1.5)

shp_SA.geometry.boundary.plot(ax = axs[0, 1], edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = axs[0, 1], edgecolor = "black", linewidth = .25)

axs[0, 1].set_ylim(-18.5, 0.5)
axs[0, 1].set_xlim(-81.75, -68)
axs[0, 1].set_ylabel("")
axs[0, 1].set_xlabel("")
axs[0, 1].set_title("")
axs[0, 1].xaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[0, 1].yaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[0, 1].grid(True, linestyle='--', color = "black", alpha = 0.1)

plot_b = tx_change_std.tx.plot(ax = axs[1, 1], cmap = "viridis", add_colorbar=False, levels =  np.arange(0, 21, 2))
axin = inset_axes(axs[1, 1], width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = axs[1, 1].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('SD de Cambio de Tx (%)', labelpad=-13, size = 3.5)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1.5)

shp_SA.geometry.boundary.plot(ax = axs[1, 1], edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = axs[1, 1], edgecolor = "black", linewidth = .25)

axs[1, 1].set_ylim(-18.5, 0.5)
axs[1, 1].set_xlim(-81.75, -68)
axs[1, 1].set_ylabel("")
axs[1, 1].set_xlabel("")
axs[1, 1].set_title("")
axs[1, 1].xaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[1, 1].yaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[1, 1].grid(True, linestyle='--', color = "black", alpha = 0.1)

# tn
plot_b = tn_change_percent.tn.plot(ax = axs[0, 2], cmap = "bwr", add_colorbar=False, levels= 6, vmax=15, vmin=-15)
axin = inset_axes(axs[0, 2], width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = axs[0, 2].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Cambio de Tn (%)', labelpad=-16.5, size = 3.5)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1.5)

shp_SA.geometry.boundary.plot(ax = axs[0, 2], edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = axs[0, 2], edgecolor = "black", linewidth = .25)

axs[0, 2].set_ylim(-18.5, 0.5)
axs[0, 2].set_xlim(-81.75, -68)
axs[0, 2].set_ylabel("")
axs[0, 2].set_xlabel("")
axs[0, 2].set_title("")
axs[0, 2].xaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[0, 2].yaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[0, 2].grid(True, linestyle='--', color = "black", alpha = 0.1)

plot_b = tn_change_std.tn.plot(ax = axs[1, 2], cmap = "viridis", add_colorbar=False, levels =  np.arange(0, 21, 2))
axin = inset_axes(axs[1, 2], width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = axs[1, 2].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('SD de Cambio de Tn (%)', labelpad=-13, size = 3.5)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1.5)

shp_SA.geometry.boundary.plot(ax = axs[1, 2], edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = axs[1, 2], edgecolor = "black", linewidth = .25)

axs[1, 2].set_ylim(-18.5, 0.5)
axs[1, 2].set_xlim(-81.75, -68)
axs[1, 2].set_ylabel("")
axs[1, 2].set_xlabel("")
axs[1, 2].set_title("")
axs[1, 2].xaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[1, 2].yaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[1, 2].grid(True, linestyle='--', color = "black", alpha = 0.1)

# pet
plot_b = pet_change_percent.pet.plot(ax = axs[0, 3], cmap = "bwr", add_colorbar=False, levels= 6, vmax=50, vmin=-50)
axin = inset_axes(axs[0, 3], width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = axs[0, 3].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('Cambio de PE (%)', labelpad=-17, size = 3.5)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1.5)

shp_SA.geometry.boundary.plot(ax = axs[0, 3], edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = axs[0, 3], edgecolor = "black", linewidth = .25)

axs[0, 3].set_ylim(-18.5, 0.5)
axs[0, 3].set_xlim(-81.75, -68)
axs[0, 3].set_ylabel("")
axs[0, 3].set_xlabel("")
axs[0, 3].set_title("")
axs[0, 3].xaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[0, 3].yaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[0, 3].grid(True, linestyle='--', color = "black", alpha = 0.1)

plot_b = pet_change_std.pet.plot(ax = axs[1, 3], cmap = "viridis", add_colorbar=False, levels =  np.arange(0, 21, 2))
axin = inset_axes(axs[1, 3], width='4%', height='35%', loc = 'lower left', bbox_to_anchor = (0.05, 0.025, 1 ,1), bbox_transform = axs[1, 3].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('SD de Cambio de PE (%)', labelpad=-13, size = 3.5)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1.5)

shp_SA.geometry.boundary.plot(ax = axs[1, 3], edgecolor = "black", linewidth = .5)
shp_dep.geometry.boundary.plot(ax = axs[1, 3], edgecolor = "black", linewidth = .25)

axs[1, 3].set_ylim(-18.5, 0.5)
axs[1, 3].set_xlim(-81.75, -68)
axs[1, 3].set_ylabel("")
axs[1, 3].set_xlabel("")
axs[1, 3].set_title("")
axs[1, 3].xaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[1, 3].yaxis.set_tick_params(labelsize = 2.5, pad = -3)
axs[1, 3].grid(True, linestyle='--', color = "black", alpha = 0.1)


plt.savefig('output/P_Tx_Tn_PE_changes.png', bbox_inches='tight',pad_inches = 0.01, dpi = 250)
plt.close()