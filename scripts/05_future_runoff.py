import xarray as xr
import numpy as np
import glob
from joblib import Parallel, delayed
from cdo import *
cdo = Cdo()

exec(open("src/calib_Budyko.py").read())

omega = xr.open_dataset("data/processed/others/omega/ensemble_omega.nc")
p_fut = xr.open_dataset("data/processed/future/prec/p_SENAMHI_pp_d12k_ACCESS1-0_rcp85_scal.nc")
pet_fut = xr.open_dataset("data/processed/future/pet/pet_SENAMHI_tmax_d12k_ACCESS1-0_rcp85_scal.nc")

def apply_q2(omega_i):
    encoding = {v: {'zlib': True, 'complevel': 5} for v in ["q"]}
    gridded_res = xr.apply_ufunc(get_q2,
                                 p_fut.p,
                                 pet_fut.pet,
                                 omega.omega.isel(time=omega_i),
                                 vectorize=True,
                                 output_dtypes=['float32'])
    gridded_res.to_dataset(name="q").to_netcdf("data/processed/future/runoff/ACCESS/ensemble/runoff_gridded_" + str(omega_i).zfill(3) + ".nc", encoding=encoding, engine='netcdf4')


Parallel(n_jobs=8, verbose=50)(
  delayed(apply_q2)(i) for i in range(len(omega.time))
)

q_fut = cdo.cat(input=sorted(glob.glob("data/processed/future/runoff/ACCESS/ensemble/runoff_gridded_*.nc")))
q_fut = xr.open_dataset(q_fut, chunks={"latitude":100,"longitude":100})
q_fut.median(dim="time").to_netcdf("data/processed/future/runoff/ACCESS/median_ensemble_runoff.nc")
q_fut.std(dim="time").to_netcdf("data/processed/future/runoff/ACCESS/sd_ensemble_runoff.nc")

###

omega = xr.open_dataset("data/processed/others/omega/ensemble_omega.nc")
p_fut = xr.open_dataset("data/processed/future/prec/p_SENAMHI_pp_d12k_HadGEM2-ES_rcp85_scal.nc")
pet_fut = xr.open_dataset("data/processed/future/pet/pet_SENAMHI_tmax_d12k_HadGEM2-ES_rcp85_scal.nc")

def apply_q2(omega_i):
    encoding = {v: {'zlib': True, 'complevel': 5} for v in ["q"]}
    gridded_res = xr.apply_ufunc(get_q2,
                                 p_fut.p,
                                 pet_fut.pet,
                                 omega.omega.isel(time=omega_i),
                                 vectorize=True,
                                 output_dtypes=['float32'])
    gridded_res.to_dataset(name="q").to_netcdf("data/processed/future/runoff/HADGEM/ensemble/runoff_gridded_" + str(omega_i).zfill(3) + ".nc", encoding=encoding, engine='netcdf4')


Parallel(n_jobs=8, verbose=50)(
  delayed(apply_q2)(i) for i in range(len(omega.time))
)

q_fut = cdo.cat(input=sorted(glob.glob("data/processed/future/runoff/HADGEM/ensemble/runoff_gridded_*.nc")))
q_fut = xr.open_dataset(q_fut, chunks={"latitude":100,"longitude":100})
q_fut.median(dim="time").to_netcdf("data/processed/future/runoff/HADGEM/median_ensemble_runoff.nc")
q_fut.std(dim="time").to_netcdf("data/processed/future/runoff/HADGEM/sd_ensemble_runoff.nc")


###

omega = xr.open_dataset("data/processed/others/omega/ensemble_omega.nc")
p_fut = xr.open_dataset("data/processed/future/prec/p_SENAMHI_pp_d12k_MPI-ESM-LR_rcp85_scal.nc")
pet_fut = xr.open_dataset("data/processed/future/pet/pet_SENAMHI_tmax_d12k_MPI-ESM-LR_rcp85_scal.nc")

def apply_q2(omega_i):
    encoding = {v: {'zlib': True, 'complevel': 5} for v in ["q"]}
    gridded_res = xr.apply_ufunc(get_q2,
                                 p_fut.p,
                                 pet_fut.pet,
                                 omega.omega.isel(time=omega_i),
                                 vectorize=True,
                                 output_dtypes=['float32'])
    gridded_res.to_dataset(name="q").to_netcdf("data/processed/future/runoff/MPIESM/ensemble/runoff_gridded_" + str(omega_i).zfill(3) + ".nc", encoding=encoding, engine='netcdf4')


Parallel(n_jobs=8, verbose=50)(
  delayed(apply_q2)(i) for i in range(len(omega.time))
)

q_fut = cdo.cat(input=sorted(glob.glob("data/processed/future/runoff/MPIESM/ensemble/runoff_gridded_*.nc")))
q_fut = xr.open_dataset(q_fut, chunks={"latitude":100,"longitude":100})
q_fut.median(dim="time").to_netcdf("data/processed/future/runoff/MPIESM/median_ensemble_runoff.nc")
q_fut.std(dim="time").to_netcdf("data/processed/future/runoff/MPIESM/sd_ensemble_runoff.nc")