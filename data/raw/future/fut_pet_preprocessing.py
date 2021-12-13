import xarray as xr
import pandas as pd
import numpy as np
import glob
import os
import urllib.request
from cdo import *

cdo = Cdo()

# HS equation
HS_equation = urllib.request.urlopen("https://git.io/J1AXN")
exec(HS_equation.read())

nc_files_tmax = glob.glob("data/raw/future/temp/SENAMHI_tmax_*.nc")
nc_files_tmin = glob.glob("data/raw/future/temp/SENAMHI_tmin_*.nc")

nc_files_tmax
nc_files_tmin

encoding = {v: {'zlib': True, 'complevel': 5} for v in ["pet"]}

PISCOsns_grid = xr.open_dataset("data/processed/present/PISCO/pet/PET_mean_anual.nc").pet

hydro_time = pd.date_range("2036-01-01", "2065-12-31", freq="D")
hydro_time = hydro_time[~((hydro_time.day == 29) & (hydro_time.month == 2))]

for nfile in range(len(nc_files_tmax)):

    piscotx = xr.open_dataset(nc_files_tmax[nfile]).sel(time=slice("2035-01-01", "2065-12-31")).rename(
        {"lat": "latitude", "lon": "longitude", "tasmax": "tx"})
    piscotn = xr.open_dataset(nc_files_tmin[nfile]).sel(time=slice("2035-01-01", "2065-12-31")).rename(
        {"lat": "latitude", "lon": "longitude", "tasmin": "tn"})

    pisco_lat = xr.DataArray(np.tile(piscotx["latitude"].values, (145, 1)).transpose(),
                             coords=[piscotx["latitude"].values, piscotx["longitude"].values],
                             dims=["latitude", "longitude"])

    for year in range(2035, 2066):
        # getting time values as Julian day
        dates = pd.date_range(str(year) + '-01-01', str(year) + "-12-31", freq='D')
        dates = np.array([int(i.strftime("%j")) for i in dates])

        xr.apply_ufunc(hargreaves_samani,  # hargreaves_samani(time_i, tmax_i, tmin_i, lat_i)
                       dates,
                       piscotx.tx.loc[str(year) + '-01-01':str(year) + "-12-31"],
                       piscotn.tn.loc[str(year) + '-01-01':str(year) + "-12-31"],
                       pisco_lat,
                       vectorize=True,
                       input_core_dims=[["time"], ["time"], ["time"], []],
                       output_core_dims=[["time"]],
                       output_dtypes=['float32']). \
            transpose("time", "latitude", "longitude"). \
            to_dataset(name="pet"). \
            to_netcdf("data/raw/future/pet/pet_" + nc_files_tmax[nfile].split("/")[
            -1] + "_" + str(year) + ".nc", encoding=encoding, engine='netcdf4')

    piscopet = xr.open_dataset(cdo.cat(input=sorted(glob.glob("data/raw/future/pet/pet_" + nc_files_tmax[nfile].split("/")[-1] + "*.nc"))))
    piscopet = piscopet.sel(time=slice("2035-09-01", "2065-08-31"))
    piscopet = piscopet.isel(time=~piscopet.time.dt.strftime('%m-%d').isin("02-29"))
    piscopet["time"] = hydro_time
    piscopet_anual = piscopet.resample(time="1Y").sum().mean(dim="time")  # 30 years
    piscopet_anual = piscopet_anual.reindex(longitude=PISCOsns_grid.longitude.values,
                                            latitude=PISCOsns_grid.latitude.values,
                                            method="nearest")
    piscopet_anual.to_netcdf("data/processed/future/pet/pet_" + nc_files_tmax[nfile].split("/")[-1], encoding=encoding, engine='netcdf4')

[os.remove(file) for file in glob.glob("data/raw/future/pet/*.nc")]