import xarray as xr
import numpy as np
from pandas import date_range as pd_date_range

#
hydro_time = pd_date_range("1982-01-01", "2016-12-31", freq="D")
hydro_time = hydro_time[~((hydro_time.day == 29) & (hydro_time.month == 2))]

#
piscopet = xr.open_dataset("data/processed/present/PISCO/pet/PET_mov15yearly.nc")

#
piscop = xr.open_dataset("data/raw/observed/PISCO/prec/PISCOpd.nc")
piscop = piscop.rename({"z":"time"})
piscop = piscop.sel(time=slice("1981-09-01", "2016-08-31"))
piscop = piscop.isel(time=~piscop.time.dt.strftime('%m-%d').isin("02-29"))
piscop["time"] = hydro_time

piscop_anual = piscop.resample(time="1Y").sum()
piscop_anual = piscop_anual.rolling(time=15, center=True).mean().dropna("time")
piscop_anual = piscop_anual.reindex(longitude=np.arange(piscopet.longitude.values[1], piscopet.longitude.values[-1], 0.008),
                                    latitude=np.arange(piscopet.latitude.values[1], piscopet.latitude.values[-1], -0.008),
                                    method="nearest")

np.round(piscop_anual, 1).to_netcdf("data/processed/present/PISCO/prec/P_mov15yearly.nc")