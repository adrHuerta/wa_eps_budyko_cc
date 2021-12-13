import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.25, rc={"lines.linewidth": 2})

exec(open("src/calib_Budyko.py").read())

# shps
shp_dep = gpd.read_file("data/raw/observed/EPSs/C_aporte.shp")
shp_SA = gpd.read_file("data/raw/shps/Sudamérica.shp").to_crs({"init": "epsg:4326"})
shp_lks = gpd.read_file("data/raw/shps/lago_titicaca_sideteva_puno.shp").to_crs({"init": "epsg:4326"})

omega_value = pd.read_csv("data/processed/others/omega/csv_omega_values.csv", index_col=0)
q_pre_median = xr.open_dataset("data/processed/present/PISCO/runoff/median_ensemble_runoff.nc")
q_pre_sd = xr.open_dataset("data/processed/present/PISCO/runoff/sd_ensemble_runoff.nc")
q_pre_cv = (q_pre_sd+.1)/(q_pre_median+.1)
Q_shp = gpd.read_file("data/processed/present/PISCO/runoff/Q_mov30yearly_shp.shp")
cluster_areas = xr.open_dataset("data/processed/others/budyko_groups.nc")

# 10 basins by each region, each region a budyko curve based on single basins (type 01)
"""
# hist omega by region

fig, axes = plt.subplots(3, 4, figsize=(40, 6), dpi = 150)
i = 0
for triaxis in axes:
    for axis in triaxis:
        axis.hist(omega_value[omega_value.columns[i]].values.tolist(), bins = 6)
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
fig.text(0.5, 0.05, r'$\omega$', ha='center', size = 9)
fig.text(0.07, 0.5, 'Frecuencia', va='center', rotation='vertical', size = 9)
# plt.savefig('output/omega_distribution_by_region.png', bbox_inches='tight',pad_inches = 0.01, dpi = 150)
plt.close()

# hist omega for "M" region
fig, ax = plt.subplots(figsize=(7, 4), dpi = 150)
ax.hist(omega_value["M"].values.tolist(), bins = 10)
ax.xaxis.set_tick_params(labelsize = 10, pad = -3)
ax.yaxis.set_tick_params(labelsize = 10, pad = -3)
ax.set_xlabel("Omega (adimensional)")
ax.set_ylabel("Frecuencia")

# plt.savefig('output/omega_distribution_region_M.png', bbox_inches='tight',pad_inches = 0.01, dpi = 150)
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

#plt.savefig('output/Q_median_cv_present_region_M.png', bbox_inches='tight',pad_inches = 0.01, dpi = 150)
plt.close()

# building budyko data of region "M" for each time step
selected_basin = pd.read_csv("data/processed/others/omega/csv_select_basin_used_for_omega_estimation.csv", index_col=0)
prec_by_basin_df = pd.read_csv("data/processed/present/PISCO/prec/P_mov30yearly_by_basin.csv", index_col=0)
pet_by_basin_df = pd.read_csv("data/processed/present/PISCO/pet/PET_mov30yearly_by_basin.csv", index_col=0)
aet_by_basin_df = pd.read_csv("data/processed/present/PISCO/aet/AET_mov30yearly_by_basin.csv", index_col=0)

Ei_Ai_by_year = []
for sb, year in zip(np.array_split(selected_basin["L"], 6), aet_by_basin_df.index):
    Ei_by_basin = aet_by_basin_df[sb].loc[year] / prec_by_basin_df[sb].loc[year]
    Ai_by_basin = pet_by_basin_df[sb].loc[year] / prec_by_basin_df[sb].loc[year]
    res = pd.concat([Ei_by_basin, Ai_by_basin], axis=1)
    res["time"] = int(year.split("-")[0])
    res.columns = ["Ei", "Ai", "Año"]
    Ei_Ai_by_year.append(res)

Ei_Ai_by_year = pd.concat(Ei_Ai_by_year, axis=0)

# one single figure

fig, axs = plt.subplots(2, 2, figsize=(6.5, 6), dpi = 150, gridspec_kw = {'wspace':0.2, 'hspace':0.275})

# budyko plot
plot_b = sns.scatterplot(y="Ei",x="Ai", hue="Año", data=Ei_Ai_by_year, ax = axs[0, 0])
plot_b.legend(fontsize=5)
axs[0, 0].plot(np.arange(0,1.1,.1), np.arange(0,1.1,.1), "black", linewidth=1.5)
axs[0, 0].plot(np.arange(1,4.1), np.repeat(1,4), "black", linewidth=1.5)
axs[0, 0].xaxis.set_tick_params(labelsize = 7.5, pad = -3)
axs[0, 0].yaxis.set_tick_params(labelsize = 7.5, pad = -3)
axs[0, 0].set_xlabel("PE/P", size = 7.5)
axs[0, 0].set_ylabel("AE/P", size = 7.5)
axs[0, 0].set_ylim(0, 1.1)
axs[0, 0].set_xlim(0, 4)
axs[0, 0].set_title("Curva de Budyko")

# omega empirical distribution
axs[0, 1].hist(omega_value["M"].values.tolist(), bins = 6)
axs[0, 1].xaxis.set_tick_params(labelsize = 7.5, pad = -3)
axs[0, 1].yaxis.set_tick_params(labelsize = 7.5, pad = -3)
axs[0, 1].set_xlabel(r'$\omega$', size = 7.5)
axs[0, 1].set_ylabel("Frecuencia", size = 7.5)
axs[0, 1].set_title("Distribución empirica de " + r'$\omega$')

# median q
plot_b = (cluster_areas.gridded_groups * q_pre_median.q).plot(ax=axs[1, 0], cmap = "viridis", add_colorbar=False, levels =  [0, 50, 100, 500, 1000, 2000, 4000])
axin = inset_axes(axs[1, 0], width='4%', height='45%', loc = 'lower left', bbox_to_anchor = (0.105, 0.065, 1 ,1), bbox_transform = axs[1, 0].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('$mm$', labelpad=-30, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0, width=.25, length=2)

shp_SA.geometry.boundary.plot(ax = axs[1, 0], edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = axs[1, 0], edgecolor = "black", linewidth = .5)

axs[1, 0].set_ylim(-19, -12)
axs[1, 0].set_xlim(-75, -67.5)
axs[1, 0].set_ylabel("")
axs[1, 0].set_xlabel("")
axs[1, 0].set_title("")
axs[1, 0].xaxis.set_tick_params(labelsize = 5.5, pad = -3)
axs[1, 0].yaxis.set_tick_params(labelsize = 5.5, pad = -3)
axs[1, 0].grid(True, linestyle='--', color = "black", alpha = 0.1)
axs[1, 0].set_title("Mediana de Escurrimiento")

# sd median q
plot_b = (cluster_areas.gridded_groups * q_pre_cv.q).plot(ax=axs[1, 1], cmap = "viridis", add_colorbar=False, levels = [.1, .2, .3, .4, .5])
axin = inset_axes(axs[1, 1], width='4%', height='45%', loc = 'lower left', bbox_to_anchor = (0.105, 0.065, 1 ,1), bbox_transform = axs[1, 1].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 5)
cb.ax.set_ylabel('(adimensional)', labelpad=-27, size = 6)
cb.ax.tick_params(labelsize = 6, pad = 0, width=.25, length=2)

shp_SA.geometry.boundary.plot(ax = axs[1, 1], edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = axs[1, 1], edgecolor = "black", linewidth = .5)

axs[1, 1].set_ylim(-19, -12)
axs[1, 1].set_xlim(-75, -67.5)
axs[1, 1].set_ylabel("")
axs[1, 1].set_xlabel("")
axs[1, 1].set_title("")
axs[1, 1].xaxis.set_tick_params(labelsize = 5.5, pad = -3)
axs[1, 1].yaxis.set_tick_params(labelsize = 5.5, pad = -3)
axs[1, 1].grid(True, linestyle='--', color = "black", alpha = 0.1)
axs[1, 1].set_title("CV de Escurrimiento")

plt.savefig('output/omega_distribution_Q_median_cv_present_region_M_type_01.png', bbox_inches='tight',pad_inches = 0.01, dpi = 150)
plt.close()
"""

# 10 basins by each region, a single budyko curve based on all basins for each time step (type 02)

# building budyko data of region "M" for each time step
selected_basin = pd.read_csv("data/processed/others/omega/csv_select_basin_used_for_omega_estimation.csv", index_col=0)
prec_by_basin_df = pd.read_csv("data/processed/present/PISCO/prec/P_mov30yearly_by_basin.csv", index_col=0)
pet_by_basin_df = pd.read_csv("data/processed/present/PISCO/pet/PET_mov30yearly_by_basin.csv", index_col=0)
aet_by_basin_df = pd.read_csv("data/processed/present/PISCO/aet/AET_mov30yearly_by_basin.csv", index_col=0)

Ei_Ai_by_year = []
for time, sb in zip(selected_basin.columns, aet_by_basin_df.index):
    Ei_by_basin = aet_by_basin_df.loc[time].loc[selected_basin[time]] / prec_by_basin_df.loc[time].loc[selected_basin[time]]
    Ai_by_basin = pet_by_basin_df.loc[time].loc[selected_basin[time]] / prec_by_basin_df.loc[time].loc[selected_basin[time]]
    res = pd.concat([Ei_by_basin, Ai_by_basin], axis=1)
    res["time"] = int(time.split("-")[0])
    res.columns = ["Ei", "Ai", "Año"]
    Ei_Ai_by_year.append(res)

Ei_Ai_by_year = pd.concat(Ei_Ai_by_year, axis=0)


fig, axs = plt.subplots(2, 2, figsize=(6.5, 6), dpi = 150, gridspec_kw = {'wspace':0.2, 'hspace':0.275})

# budyko plot
plot_b = sns.scatterplot(y="Ei",x="Ai", hue="Año", data=Ei_Ai_by_year, ax = axs[0, 0])
plot_b.legend(fontsize=5)
axs[0, 0].plot(np.arange(0,1.1,.1), np.arange(0,1.1,.1), "black", linewidth=1.5)
axs[0, 0].plot(np.arange(1,4.1), np.repeat(1,4), "black", linewidth=1.5)
axs[0, 0].xaxis.set_tick_params(labelsize = 7.5, pad = -3)
axs[0, 0].yaxis.set_tick_params(labelsize = 7.5, pad = -3)
axs[0, 0].set_xlabel("PE/P", size = 7.5)
axs[0, 0].set_ylabel("AE/P", size = 7.5)
axs[0, 0].set_ylim(0, 1.1)
axs[0, 0].set_xlim(0, 4)
axs[0, 0].set_title("Curva de Budyko")

# omega empirical distribution
axs[0, 1].hist(omega_value, bins = 6)
axs[0, 1].xaxis.set_tick_params(labelsize = 7.5, pad = -3)
axs[0, 1].yaxis.set_tick_params(labelsize = 7.5, pad = -3)
axs[0, 1].set_xlabel(r'$\omega$', size = 7.5)
axs[0, 1].set_ylabel("Frecuencia", size = 7.5)
axs[0, 1].set_title("Distribución empirica de " + r'$\omega$')

# median q
plot_b = (q_pre_median.q).plot(ax=axs[1, 0], cmap = "viridis", add_colorbar=False, levels =  [0, 50, 100, 500, 1000, 2000, 4000])
axin = inset_axes(axs[1, 0], width='4%', height='45%', loc = 'lower left', bbox_to_anchor = (0.0105, 0.005, 1 ,1), bbox_transform = axs[1, 0].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 4)
cb.ax.tick_params(labelsize = 4, pad = 0, width=.25, length=1)
cb.ax.set_ylabel('$mm$', labelpad=-20, size = 6)

shp_SA.geometry.boundary.plot(ax = axs[1, 0], edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = axs[1, 0], edgecolor = "black", linewidth = .5)

axs[1, 1].set_ylim(-18.5, 0.5)
axs[1, 1].set_xlim(-81.75, -68)
axs[1, 0].set_ylabel("")
axs[1, 0].set_xlabel("")
axs[1, 0].set_title("")
axs[1, 0].xaxis.set_tick_params(labelsize = 5.5, pad = -3)
axs[1, 0].yaxis.set_tick_params(labelsize = 5.5, pad = -3)
axs[1, 0].grid(True, linestyle='--', color = "black", alpha = 0.1)
axs[1, 0].set_title("Mediana de Escurrimiento")

# sd median q
plot_b = (q_pre_cv.q).plot(ax=axs[1, 1], cmap = "viridis", add_colorbar=False, levels = np.arange(0, 1.1, .1))
axin = inset_axes(axs[1, 1], width='4%', height='40%', loc = 'lower left', bbox_to_anchor = (0.055, 0.025, 1 ,1), bbox_transform = axs[1, 1].transAxes)
cb = plt.colorbar(plot_b, cax=axin, orientation = "vertical", aspect = 4)
cb.ax.tick_params(labelsize = 5, pad = 0, width=.25, length=1)
cb.ax.set_ylabel('(adimensional)', labelpad=-20, size = 6)

shp_SA.geometry.boundary.plot(ax = axs[1, 1], edgecolor = "black", linewidth = .65)
shp_dep.geometry.boundary.plot(ax = axs[1, 1], edgecolor = "black", linewidth = .5)

axs[1, 1].set_ylim(-18.5, 0.5)
axs[1, 1].set_xlim(-81.75, -68)
axs[1, 1].set_ylabel("")
axs[1, 1].set_xlabel("")
axs[1, 1].set_title("")
axs[1, 1].xaxis.set_tick_params(labelsize = 5.5, pad = -3)
axs[1, 1].yaxis.set_tick_params(labelsize = 5.5, pad = -3)
axs[1, 1].grid(True, linestyle='--', color = "black", alpha = 0.1)
axs[1, 1].set_title("CV de Escurrimiento")

plt.savefig('output/omega_distribution_Q_median_cv_present_region_M_type_02.png', bbox_inches='tight',pad_inches = 0.01, dpi = 150)
plt.close()