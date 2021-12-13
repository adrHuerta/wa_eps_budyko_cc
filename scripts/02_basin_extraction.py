import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import rioxarray
from joblib import Parallel, delayed

exec(open("src/crop_mask.py").read())

pisco_pet = xr.open_dataset("data/processed/present/PISCO/pet/PET_mov30yearly.nc")
pisco_p = xr.open_dataset("data/processed/present/PISCO/prec/P_mov30yearly.nc")
pisco_q_values = pd.read_csv("data/processed/present/PISCO/runoff/Q_mov30yearly.csv", index_col=0)
pisco_q_values_shp = gpd.read_file("data/processed/present/PISCO/runoff/Q_mov30yearly_shp.shp")

def extract_prec(i_basin):
    import rioxarray
    subx = pisco_q_values_shp.iloc[i_basin: i_basin + 1]
    # prec
    pisco_p_cropped = xr_crop(shp_i=subx, netcdf_i=pisco_p)
    shp_exp_grid = xr_shp_to_grid(shp_i=subx, netcdf_array=pisco_p_cropped.isel(time=0).p)
    pisco_p_masked = xr_mask(grid_mask=shp_exp_grid, netcdf_i=pisco_p_cropped)
    return pisco_p_masked.mean(dim=["latitude", "longitude"]).to_dataframe()

def extract_pet(i_basin):
    import rioxarray
    subx = pisco_q_values_shp.iloc[i_basin: i_basin + 1]
    # pet
    pisco_pet_cropped = xr_crop(shp_i=subx, netcdf_i=pisco_pet)
    shp_exp_grid = xr_shp_to_grid(shp_i=subx, netcdf_array=pisco_pet_cropped.isel(time=0).pet)
    pisco_pet_masked = xr_mask(grid_mask=shp_exp_grid, netcdf_i=pisco_pet_cropped)
    return pisco_pet_masked.mean(dim=["latitude", "longitude"]).to_dataframe()

num_cores = 6
prec_by_basin = Parallel(n_jobs=num_cores, verbose=50)(delayed(extract_prec)(i) for i in range(len(pisco_q_values_shp.GR2M_ID)))
pet_by_basin = Parallel(n_jobs=num_cores, verbose=50)(delayed(extract_pet)(i) for i in range(len(pisco_q_values_shp.GR2M_ID)))

prec_by_basin_df = pd.concat(prec_by_basin, axis=1)
prec_by_basin_df.columns = pisco_q_values.columns
pet_by_basin_df = pd.concat(pet_by_basin, axis=1)
pet_by_basin_df.columns = pisco_q_values.columns
aet_by_basin = prec_by_basin_df - pisco_q_values # AE = P - Q

Ei_by_basin = aet_by_basin / prec_by_basin_df
Ai_by_basin = pet_by_basin_df / prec_by_basin_df
enery_limit = aet_by_basin / pet_by_basin_df
water_limit = Ei_by_basin

aet_by_basin_df = aet_by_basin
# removing "bad" basins (outside limits of Budyko's space)
for i in aet_by_basin.index:
    bad_basins_00_i = enery_limit.loc[i][enery_limit.loc[i] > 1].index.tolist() # max(AE/PE) = 1
    bad_basins_01_i = water_limit.loc[i][water_limit.loc[i] > 1].index.tolist() # max(AE/P) = 1
    bad_basins_02_i = Ei_by_basin.loc[i][Ei_by_basin.loc[i] < 0].index.tolist() # min(Ei) = 0
    bad_basins_i = list(set(bad_basins_00_i + bad_basins_01_i + bad_basins_02_i))
    prec_by_basin_df.loc[i, bad_basins_i] = np.nan
    pet_by_basin_df.loc[i, bad_basins_i] = np.nan
    aet_by_basin_df.loc[i, bad_basins_i] = np.nan
    # print(len(bad_basins_i))

# deleted basin: 770, 775, 775, 777, 779, 776
prec_by_basin_df.to_csv("data/processed/present/PISCO/prec/P_mov30yearly_by_basin.csv")
pet_by_basin_df.to_csv("data/processed/present/PISCO/pet/PET_mov30yearly_by_basin.csv")
aet_by_basin_df.to_csv("data/processed/present/PISCO/aet/AET_mov30yearly_by_basin.csv")
