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

P_by_basin = pd.read_csv("data/processed/present/PISCO/prec/P_mov15yearly_by_basin.csv", index_col=0)
PET_by_basin = pd.read_csv("data/processed/present/PISCO/pet/PET_mov15yearly_by_basin.csv", index_col=0)
Q_by_basin = pd.read_csv("data/processed/present/PISCO/runoff/Q_mov15yearly.csv", index_col=0)
Q_shp = gpd.read_file("data/processed/present/PISCO/runoff/Q_mov15yearly_shp.shp")
cluster_areas = xr.open_dataset("data/processed/others/budyko_groups.nc")

shp_order = sorted(Q_shp.Region.unique())
for i in range(len(shp_order)):
    Q_shp.loc[(Q_shp.Region == shp_order[i]), "Region"] = i + 1

# getting best omega
best_omega_by_time = []
for time in Q_by_basin.index:

    P_per_year = P_by_basin.loc[time]
    PET_per_year = PET_by_basin.loc[time]
    Q_per_year = Q_by_basin.loc[time]

    best_omega_by_area = []
    for area in sorted(Q_shp.Region.unique()):
        selected_basins = random.sample(["GR2M_ID_" + str(item) for item in Q_shp[Q_shp["Region"] == area]["GR2M_ID"].tolist()], 30)
        best_omega_by_area.append(pd.DataFrame([calib_Budyko(q=Q_per_year.loc[item], p=P_per_year.loc[item], pet=PET_per_year.loc[item])  for item in selected_basins]))

    best_omega_by_time.append(pd.concat(best_omega_by_area, axis=1))

# making omega to grid

best_omega = pd.concat(best_omega_by_time, axis=0)
best_omega.reset_index(inplace=True, drop=True)

encoding = {v: {'zlib': True, 'complevel': 5} for v in ["omega"]}

for ensemble in range(best_omega.shape[0]):
    ensemble_grid = cluster_areas.gridded_groups
    ensemble_omega = best_omega.iloc[ensemble]
    for a, b in zip([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], ensemble_omega):
        ensemble_grid.values[ensemble_grid.values == a] = b
    ensemble_grid.to_dataset(name="omega").to_netcdf("data/processed/others/omega/omega_" + str(ensemble).zfill(3) + ".nc", encoding=encoding, engine='netcdf4')

cdo.cat(input=sorted(glob.glob("data/processed/others/omega/omega_*.nc")),
        output="data/processed/others/omega/ensemble_omega.nc")
[os.remove(item) for item in glob.glob("data/processed/others/omega/omega_*.nc")]