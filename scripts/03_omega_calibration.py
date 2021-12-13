import pandas as pd
import numpy as np
import os
import glob
import xarray as xr
import geopandas as gpd
import random
from cdo import *
cdo = Cdo()

exec(open("src/calib_Budyko.py").read())

P_by_basin = pd.read_csv("data/processed/present/PISCO/prec/P_mov30yearly_by_basin.csv", index_col=0)
PET_by_basin = pd.read_csv("data/processed/present/PISCO/pet/PET_mov30yearly_by_basin.csv", index_col=0)
Q_by_basin = pd.read_csv("data/processed/present/PISCO/runoff/Q_mov30yearly.csv", index_col=0)
Q_shp = gpd.read_file("data/processed/present/PISCO/runoff/Q_mov30yearly_shp.shp")
cluster_areas = xr.open_dataset("data/processed/others/budyko_groups.nc")

shp_order = sorted(Q_shp.Region.unique())
for i in range(len(shp_order)):
    Q_shp.loc[(Q_shp.Region == shp_order[i]), "Region"] = i + 1

# 10 basins by each region, each region a budyko curve based on single basins (type_01)
"""
# getting best omega
best_omega_by_time = []
selected_basins_by_time = []
for time in Q_by_basin.index:

    P_per_year = P_by_basin.loc[time]
    PET_per_year = PET_by_basin.loc[time]
    Q_per_year = Q_by_basin.loc[time]

    Q_shp_per_year = Q_shp
    bad_basin = P_per_year[P_per_year.isna()].index
    bad_basin = [int(item.split("_")[-1]) for item in bad_basin]
    Q_shp_per_year = Q_shp_per_year[~Q_shp_per_year["GR2M_ID"].isin(bad_basin)]

    best_omega_by_area = []
    selected_basins_by_area = []
    for area in sorted(Q_shp_per_year.Region.unique()):
        random.seed(767)
        selected_basins = random.sample(["GR2M_ID_" + str(item) for item in Q_shp_per_year[Q_shp_per_year["Region"] == area]["GR2M_ID"].tolist()], 10)
        best_omega_by_area.append(pd.DataFrame([calib_Budyko(q=Q_per_year.loc[item], p=P_per_year.loc[item], pet=PET_per_year.loc[item])  for item in selected_basins]))
        selected_basins_by_area.append(pd.DataFrame(selected_basins))

    best_omega_by_time.append(pd.concat(best_omega_by_area, axis=1))
    selected_basins_by_time.append(pd.concat(selected_basins_by_area, axis=1))

# saving selected basin IDs
select_basin = pd.concat(selected_basins_by_time, axis=0)
select_basin.reset_index(inplace=True, drop=True)
select_basin.columns = shp_order
select_basin.to_csv("data/processed/others/omega/csv_select_basin_used_for_omega_estimation.csv")

# making omega to grid

best_omega = pd.concat(best_omega_by_time, axis=0)
best_omega.reset_index(inplace=True, drop=True)
best_omega.columns = shp_order
best_omega.to_csv("data/processed/others/omega/csv_omega_values.csv")

encoding = {v: {'zlib': True, 'complevel': 5} for v in ["omega"]}

for ensemble in range(best_omega.shape[0]):
    ensemble_grid = cluster_areas.gridded_groups
    ensemble_omega = best_omega.iloc[ensemble]
    for a, b in zip([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], ensemble_omega):
        ensemble_grid.values[ensemble_grid.values == a] = b
    ensemble_grid.to_dataset(name="omega").to_netcdf("data/processed/others/omega/omega_" + str(ensemble).zfill(3) + ".nc", encoding=encoding, engine='netcdf4')

cdo.cat(input=sorted(glob.glob("data/processed/others/omega/omega_*.nc")),
        output="data/processed/others/omega/ensemble_omega.nc")
[os.remove(item) for item in glob.glob("data/processed/others/omega/omega_*.nc")]
"""


# 10 basins by each region, a single budyko curve based on all basins for each time step (type_02)
best_omega_by_time = []
selected_basins_by_time = []
for time in Q_by_basin.index:

    P_per_year = P_by_basin.loc[time]
    PET_per_year = PET_by_basin.loc[time]
    Q_per_year = Q_by_basin.loc[time]

    Q_shp_per_year = Q_shp
    bad_basin = P_per_year[P_per_year.isna()].index
    bad_basin = [int(item.split("_")[-1]) for item in bad_basin]
    Q_shp_per_year = Q_shp_per_year[~Q_shp_per_year["GR2M_ID"].isin(bad_basin)]

    best_omega_by_seed = []
    selected_basins_by_seed = []
    for seed_random in [1,2,3,4,5,6,7,8,9,10]:
        random.seed(seed_random)
        selected_basins_by_area = []
        for area in sorted(Q_shp_per_year.Region.unique()):
            selected_basins = random.sample(
                ["GR2M_ID_" + str(item) for item in Q_shp_per_year[Q_shp_per_year["Region"] == area]["GR2M_ID"].tolist()],
                5)
            selected_basins_by_area.append(selected_basins)
        selected_basins = [item for sublist in selected_basins_by_area for item in sublist]
        selected_basins_by_seed.append(pd.DataFrame(selected_basins))

        nash_stat_value = []
        for omega_w in np.arange(1.1, 5, .01):
            Q_model = pd.DataFrame([get_q(p=P_per_year.loc[item], pet=PET_per_year.loc[item], w = omega_w)  for item in selected_basins])
            Q_obs = Q_per_year[selected_basins]
            Q_obs = pd.DataFrame(Q_obs)
            Q_model.index = Q_obs.index
            Q_model.columns = ["Q"]
            Q_obs.columns = ["Q"]
            nash_stat =  1 - (np.sum(np.power(Q_obs - Q_model, 2))/np.sum(np.power(Q_obs - np.mean(Q_obs), 2)))
            nash_stat_value.append(nash_stat)

        best_omega_by_seed.append(np.arange(1.1, 11, .01)[np.argmax(nash_stat_value)])

    best_omega_by_time.append(pd.DataFrame(best_omega_by_seed))
    selected_basins_by_time.append(pd.concat(selected_basins_by_seed, axis=0))


# saving selected basin IDs
select_basin = pd.concat(selected_basins_by_time, axis=1)
select_basin.reset_index(inplace=True, drop=True)
select_basin.columns = Q_by_basin.index
select_basin.to_csv("data/processed/others/omega/csv_select_basin_used_for_omega_estimation.csv")

# making omega to grid
best_omega = pd.concat(best_omega_by_time, axis=0)
best_omega.reset_index(inplace=True, drop=True)
best_omega.to_csv("data/processed/others/omega/csv_omega_values.csv")

encoding = {v: {'zlib': True, 'complevel': 5} for v in ["omega"]}

for ensemble in range(best_omega.shape[0]):
    ensemble_grid = cluster_areas.gridded_groups
    ensemble_omega = best_omega.iloc[ensemble]
    ensemble_grid.values[ensemble_grid.values >= 0] = ensemble_omega
    ensemble_grid.to_dataset(name="omega").to_netcdf("data/processed/others/omega/omega_" + str(ensemble).zfill(3) + ".nc", encoding=encoding, engine='netcdf4')

cdo.cat(input=sorted(glob.glob("data/processed/others/omega/omega_*.nc")),
        output="data/processed/others/omega/ensemble_omega.nc")
[os.remove(item) for item in glob.glob("data/processed/others/omega/omega_*.nc")]