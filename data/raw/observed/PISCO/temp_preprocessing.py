import xarray as xr
import numpy as np
from pandas import date_range as pd_date_range

#
hydro_time = pd_date_range("1982-01-01", "2016-12-31", freq="D")
hydro_time = hydro_time[~((hydro_time.day == 29) & (hydro_time.month == 2))]

#
piscopet = xr.open_dataset("data/processed/present/PISCO/pet/PET_mov30yearly.nc")

#
piscotx = xr.open_dataset("data/raw/observed/PISCO/temp/PISCOdtx_v1.1.nc")
piscotx = piscotx.sel(time=slice("1981-09-01", "2016-08-31"))
piscotx = piscotx.isel(time=~piscotx.time.dt.strftime('%m-%d').isin("02-29"))
piscotx["time"] = hydro_time
#
piscotn = xr.open_dataset("data/raw/observed/PISCO/temp/PISCOdtn_v1.1.nc")
piscotn = piscotn.sel(time=slice("1981-09-01", "2016-08-31"))
piscotn = piscotn.isel(time=~piscotn.time.dt.strftime('%m-%d').isin("02-29"))
piscotn["time"] = hydro_time

piscotx_anual = piscotx.resample(time="1Y").mean()
piscotn_anual = piscotn.resample(time="1Y").mean()

# for analysis
piscotx_mean_anual = piscotx_anual.sel(time=slice("1982-01-01", "2011-12-31")).mean(dim="time") # 30 years
piscotx_mean_anual = piscotx_mean_anual.reindex(longitude=piscopet.longitude.values,
                                                latitude=piscopet.latitude.values,
                                                method="nearest")
np.round(piscotx_mean_anual, 1).to_netcdf("data/processed/present/PISCO/temp/TX_mean_anual.nc") # 30 years

piscotn_mean_anual = piscotn_anual.sel(time=slice("1982-01-01", "2011-12-31")).mean(dim="time")
piscotn_mean_anual = piscotn_mean_anual.reindex(longitude=piscopet.longitude.values,
                                                latitude=piscopet.latitude.values,
                                                method="nearest")
np.round(piscotn_mean_anual, 1).to_netcdf("data/processed/present/PISCO/temp/TN_mean_anual.nc")