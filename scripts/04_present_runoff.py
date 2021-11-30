import xarray as xr
import numpy as np
import glob
from joblib import Parallel, delayed
from cdo import *
cdo = Cdo()

exec(open("src/calib_Budyko.py").read())

omega = xr.open_dataset("data/processed/others/omega/ensemble_omega.nc")
p_pre = xr.open_dataset("data/processed/present/PISCO/prec/P_mean_anual.nc")
pet_pre = xr.open_dataset("data/processed/present/PISCO/pet/PET_mean_anual.nc")

def apply_q2(omega_i):
    encoding = {v: {'zlib': True, 'complevel': 5} for v in ["q"]}
    gridded_res = xr.apply_ufunc(get_q2,
                                 p_pre.p,
                                 pet_pre.pet,
                                 omega.omega.isel(time=omega_i),
                                 vectorize=True,
                                 output_dtypes=['float32'])
    gridded_res.to_dataset(name="q").to_netcdf("data/processed/present/PISCO/runoff/ensemble/runoff_gridded_" + str(omega_i).zfill(3) + ".nc", encoding=encoding, engine='netcdf4')


Parallel(n_jobs=8, verbose=50)(
  delayed(apply_q2)(i) for i in range(len(omega.time))
)

#
q_pre = cdo.cat(input=sorted(glob.glob("data/processed/present/PISCO/runoff/ensemble/runoff_gridded_*.nc")))
q_pre = xr.open_dataset(q_pre, chunks={"latitude":100,"longitude":100})
q_pre.median(dim="time").to_netcdf("data/processed/present/PISCO/runoff/median_ensemble_runoff.nc")
q_pre.std(dim="time").to_netcdf("data/processed/present/PISCO/runoff/sd_ensemble_runoff.nc")
