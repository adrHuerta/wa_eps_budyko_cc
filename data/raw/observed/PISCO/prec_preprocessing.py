import xarray as xr
import numpy as np
from pandas import date_range as pd_date_range

#
hydro_time = pd_date_range("1982-01-01", "2016-12-31", freq="D")
hydro_time = hydro_time[~((hydro_time.day == 29) & (hydro_time.month == 2))]

#
piscopet = xr.open_dataset("data/processed/present/PISCO/pet/PET_mov30yearly.nc")

#
piscop = xr.open_dataset("data/raw/observed/PISCO/prec/PISCOpd.nc")
piscop = piscop.rename({"z":"time"})
piscop = piscop.sel(time=slice("1981-09-01", "2016-08-31"))
piscop = piscop.isel(time=~piscop.time.dt.strftime('%m-%d').isin("02-29"))
piscop["time"] = hydro_time

piscop_anual = piscop.resample(time="1Y").sum()

# for analysis
piscop_mean_anual = piscop_anual.sel(time=slice("1982-01-01", "2011-12-31")).mean(dim="time") # 30 years
piscop_mean_anual = piscop_mean_anual.reindex(longitude=piscopet.longitude.values,
                                              latitude=piscopet.latitude.values,
                                              method="nearest")
np.round(piscop_mean_anual, 1).to_netcdf("data/processed/present/PISCO/prec/P_mean_anual.nc")

# for omega calibration
piscop_anual = piscop_anual.rolling(time=30, center=True).mean().dropna("time")
piscop_anual = piscop_anual.reindex(longitude=piscopet.longitude.values,
                                    latitude=piscopet.latitude.values,
                                    method="nearest")

np.round(piscop_anual, 1).to_netcdf("data/processed/present/PISCO/prec/P_mov30yearly.nc")