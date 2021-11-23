import xarray as xr
import numpy as np
import os
import glob
import urllib.request
from pandas import date_range as pd_date_range
from cdo import *
cdo = Cdo()

# HS equation
HS_equation = urllib.request.urlopen("https://git.io/J1AXN")
exec(HS_equation.read())

#
hydro_time = pd_date_range("1982-01-01", "2016-12-31", freq="D")
hydro_time = hydro_time[~((hydro_time.day == 29) & (hydro_time.month == 2))]

# getting data
piscotx = xr.open_dataset("data/raw/observed/PISCO/temp/PISCOdtx_v1.1.nc")
piscotn = xr.open_dataset("data/raw/observed/PISCO/temp/PISCOdtn_v1.1.nc")

# building lat grid
pisco_lat = xr.DataArray(np.tile(piscotx["latitude"].values, (145, 1)).transpose(),
                         coords=[piscotx["latitude"].values, piscotx["longitude"].values],
                         dims=["latitude", "longitude"])

# getting time values as Julian day
dates = pd_date_range('1981-01-01', "2016-12-31", freq='D')
dates = np.array([int(i.strftime("%j")) for i in dates])

# computing PET
for year in range(1981, 2017):
    dates = pd_date_range(str(year) + '-01-01', str(year) + "-12-31", freq='D')
    dates = np.array([int(i.strftime("%j")) for i in dates])

    xr.apply_ufunc(hargreaves_samani,
                   piscotx.tx.loc[str(year) + '-01-01':str(year) + "-12-31"],
                   piscotn.tn.loc[str(year) + '-01-01':str(year) + "-12-31"],
                   dates,
                   pisco_lat,
                   vectorize=True,
                   input_core_dims=[["time"], ["time"], ["time"], []],
                   output_core_dims=[["time"]],
                   output_dtypes=['float32']).\
        transpose("time", "latitude", "longitude").\
        to_dataset(name="pet").\
        to_netcdf("data/processed/present/PISCO/pet/pet_" + str(year) + ".nc")

#
piscopet = xr.open_dataset(cdo.cat(input=sorted(glob.glob("data/processed/present/PISCO/pet/pet_" "*.nc"))))
piscopet = piscopet.sel(time=slice("1981-09-01", "2016-08-31"))
piscopet = piscopet.isel(time=~piscopet.time.dt.strftime('%m-%d').isin("02-29"))
piscopet["time"] = hydro_time

piscopet_anual = piscopet.resample(time="1Y").sum()
piscopet_anual = piscopet_anual.rolling(time=15, center=True).mean().dropna("time")
piscopet_anual = piscopet_anual.reindex(longitude=np.arange(piscopet.longitude.values[1], piscopet.longitude.values[-1], 0.008),
                                        latitude=np.arange(piscopet.latitude.values[1], piscopet.latitude.values[-1], -0.008),
                                        method="nearest")

np.round(piscopet_anual, 1).to_netcdf("data/processed/present/PISCO/pet/PET_mov15yearly.nc")
[os.remove(file) for file in glob.glob("data/processed/present/PISCO/pet/pet_" "*.nc")]