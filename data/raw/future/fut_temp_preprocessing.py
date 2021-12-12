# -*- coding: utf-8 -*-
"""fut_temp_preprocessing.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ukfUHEq5xTWSIJmqicWshK_loj-Oe96m
"""

from google.colab import drive
drive.mount('/content/drive')

!pip install netcdf4

import xarray as xr
import pandas as pd
import numpy as np
import glob

nc_files_tmax = glob.glob("/content/drive/MyDrive/MOD_CC_SENAMHI_DAY/SENAMHI_tmax_*_rcp85_scal.nc")
nc_files_tmin = glob.glob("/content/drive/MyDrive/MOD_CC_SENAMHI_DAY/SENAMHI_tmin_*_rcp85_scal.nc")

encoding_tx = {v: {'zlib': True, 'complevel': 5} for v in ["tx"]}
encoding_tn = {v: {'zlib': True, 'complevel': 5} for v in ["tn"]}

PISCOsns_grid = xr.open_dataset("/content/drive/MyDrive/Google_Colab_temp/others/SNS/grid_SW.nc")
PISCOsns_grid

hydro_time = pd.date_range("2036-01-01", "2065-12-31", freq="D")
hydro_time = hydro_time[~((hydro_time.day == 29) & (hydro_time.month == 2))]

for nc_file in nc_files_tmax:
  
  pisco_tx = xr.open_dataset(nc_file)
  pisco_tx = pisco_tx.rename({"lat":"latitude","lon":"longitude","tasmax":"tx"})
  pisco_tx = pisco_tx.sel(time=slice("2035-09-01", "2065-08-31"))
  pisco_tx = pisco_tx.isel(time=~pisco_tx.time.dt.strftime('%m-%d').isin("02-29"))
  pisco_tx["time"] = hydro_time
  pisco_tx_anual = pisco_tx.resample(time="1Y").mean().mean(dim="time")
  pisco_tx_anual = pisco_tx_anual.reindex(longitude=PISCOsns_grid.longitude.values,
                                          latitude=PISCOsns_grid.latitude.values,
                                          method="nearest")
  pisco_tx_anual.to_netcdf("/content/drive/MyDrive/Google_Colab_temp/others/SNS/tx_" + nc_file.split("/")[-1], encoding=encoding_tx, engine='netcdf4')

for nc_file in nc_files_tmin:
  
  pisco_tn = xr.open_dataset(nc_file)
  pisco_tn = pisco_tn.rename({"lat":"latitude","lon":"longitude","tasmin":"tn"})
  pisco_tn = pisco_tn.sel(time=slice("2035-09-01", "2065-08-31"))
  pisco_tn = pisco_tn.isel(time=~pisco_tn.time.dt.strftime('%m-%d').isin("02-29"))
  pisco_tn["time"] = hydro_time
  pisco_tn_anual = pisco_tn.resample(time="1Y").mean().mean(dim="time")
  pisco_tn_anual = pisco_tn_anual.reindex(longitude=PISCOsns_grid.longitude.values,
                                          latitude=PISCOsns_grid.latitude.values,
                                          method="nearest")
  pisco_tn_anual.to_netcdf("/content/drive/MyDrive/Google_Colab_temp/others/SNS/tn_" + nc_file.split("/")[-1], encoding=encoding_tn, engine='netcdf4')

