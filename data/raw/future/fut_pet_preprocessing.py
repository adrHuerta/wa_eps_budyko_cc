# -*- coding: utf-8 -*-
"""fut_pet_preprocessing.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/16zK0ECHDNz1_D09QJJY7PUwACdMvNQwQ
"""

from google.colab import drive
drive.mount('/content/drive')

!apt install cdo
!pip install cdo

import xarray as xr
import pandas as pd
import numpy as np
import glob
import urllib.request
import glob
from cdo import *
cdo = Cdo()

# HS equation
HS_equation = urllib.request.urlopen("https://git.io/J1AXN")
exec(HS_equation.read())

nc_files_tmax = glob.glob("/content/drive/MyDrive/MOD_CC_SENAMHI_DAY/SENAMHI_tmax_*_rcp85_scal.nc")
nc_files_tmin = glob.glob("/content/drive/MyDrive/MOD_CC_SENAMHI_DAY/SENAMHI_tmin_*_rcp85_scal.nc")

nc_files_tmax

nc_files_tmin

encoding = {v: {'zlib': True, 'complevel': 5} for v in ["pet"]}

# getting time values as Julian day
dates = pd.date_range('2035-01-01', "2065-12-31", freq='D')
dates = np.array([int(i.strftime("%j")) for i in dates])

PISCOsns_grid = xr.open_dataset("/content/drive/MyDrive/Google_Colab_temp/others/SNS/grid_SW.nc")

hydro_time = pd.date_range("2036-01-01", "2065-12-31", freq="D")
hydro_time = hydro_time[~((hydro_time.day == 29) & (hydro_time.month == 2))]

for nfile in range(len(nc_files_tmax)):

  piscotx = xr.open_dataset(nc_files_tmax[nfile]).sel(time=slice("2035-01-01", "2065-12-31")).rename({"lat":"latitude","lon":"longitude","tasmax":"tx"})
  piscotn = xr.open_dataset(nc_files_tmin[nfile]).sel(time=slice("2035-01-01", "2065-12-31")).rename({"lat":"latitude","lon":"longitude","tasmin":"tn"})

  pisco_lat = xr.DataArray(np.tile(piscotx["latitude"].values, (145, 1)).transpose(),
                           coords=[piscotx["latitude"].values, piscotx["longitude"].values],
                           dims=["latitude", "longitude"])
  
  for year in range(2035, 2066):
    dates = pd.date_range(str(year) + '-01-01', str(year) + "-12-31", freq='D')
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
        to_netcdf("/content/drive/MyDrive/Google_Colab_temp/others/SNS/pet/pet_" + nc_files_tmax[nfile].split("/")[-1] + "_" + str(year) + ".nc", encoding=encoding, engine='netcdf4')
  
  piscopet = xr.open_dataset(cdo.cat(input=sorted(glob.glob("/content/drive/MyDrive/Google_Colab_temp/others/SNS/pet/pet_" + nc_files_tmax[nfile].split("/")[-1] + "*.nc"))))
  piscopet = piscopet.sel(time=slice("2035-09-01", "2065-08-31"))
  piscopet = piscopet.isel(time=~piscopet.time.dt.strftime('%m-%d').isin("02-29"))
  piscopet["time"] = hydro_time
  piscopet_anual = piscopet.resample(time="1Y").sum().mean(dim="time")
  piscopet_anual = piscopet_anual.reindex(longitude=np.arange(PISCOsns_grid.longitude.values[1], PISCOsns_grid.longitude.values[-1], 0.008),
                                          latitude=np.arange(PISCOsns_grid.latitude.values[1], PISCOsns_grid.latitude.values[-1], -0.008),
                                          method="nearest")
  piscopet_anual.to_netcdf("/content/drive/MyDrive/Google_Colab_temp/others/SNS/pet_" + nc_files_tmax[nfile].split("/")[-1], encoding=encoding, engine='netcdf4')