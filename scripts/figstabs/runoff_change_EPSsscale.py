import geopandas as gpd
import pandas as pd
import xarray as xr
import rioxarray
import numpy as np

exec(open("src/crop_mask.py").read())

shp_EPSs_01 = gpd.read_file("data/raw/observed/EPSs/cuencasEPS1.shp")
shp_EPSs_02 = gpd.read_file("data/raw/observed/EPSs/cuencasEPS2.shp")
access_q_fut = xr.open_dataset("data/processed/future/runoff/ACCESS/median_ensemble_runoff.nc")
hadgem_q_fut = xr.open_dataset("data/processed/future/runoff/HADGEM/median_ensemble_runoff.nc")
mpiesm_q_fut = xr.open_dataset("data/processed/future/runoff/MPIESM/median_ensemble_runoff.nc")
q_pre = xr.open_dataset("data/processed/present/PISCO/runoff/median_ensemble_runoff.nc")

basins_EPSs_01 = []
for basin in shp_EPSs_01.OBJECTID.unique():
    basin_i = shp_EPSs_01[shp_EPSs_01["OBJECTID"] == basin]
    #
    q_pre_i = xr_crop(shp_i=basin_i, netcdf_i=q_pre)
    q_pre_i = xr_mask(grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=q_pre_i.q, to_drop=["spatial_ref"]), netcdf_i=q_pre_i)
    q_pre_i = float(q_pre_i.q.mean())
    #
    access_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=access_q_fut)
    access_q_fut_i = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=access_q_fut_i.q, to_drop=["spatial_ref"]),
        netcdf_i=access_q_fut_i)
    access_q_fut_i = float(access_q_fut_i.q.mean())
    #
    hadgem_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=hadgem_q_fut)
    hadgem_q_fut_i = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=hadgem_q_fut_i.q, to_drop=["spatial_ref"]),
        netcdf_i=hadgem_q_fut_i)
    hadgem_q_fut_i = float(hadgem_q_fut_i.q.mean())
    #
    mpiesm_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=mpiesm_q_fut)
    mpiesm_q_fut_i = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=mpiesm_q_fut.q, to_drop=["spatial_ref"]),
        netcdf_i=mpiesm_q_fut_i)
    mpiesm_q_fut_i = float(mpiesm_q_fut_i.q.mean())

    q_change_i = [100 * (access_q_fut_i - (q_pre_i + 1)) / (q_pre_i + 1),
                  100 * (hadgem_q_fut_i - (q_pre_i + 1)) / (q_pre_i + 1),
                  100 * (mpiesm_q_fut_i - (q_pre_i + 1)) / (q_pre_i + 1)]

    q_change_median_i = np.round(np.median(q_change_i),1)
    q_change_std_i = np.round(np.std(q_change_i), 1)

    res = pd.DataFrame({"Code":[int(basin_i.CODIGO)],
                        "Region": basin_i.NOMB_UH_N1.tolist()[0].split(" ")[-1],
                        "Q_present":[np.round(q_pre_i, 1)],
                        "ACCESS1-0":[np.round(access_q_fut_i, 1)],
                        "HadGEM2-ES":[np.round(hadgem_q_fut_i, 1)],
                        "MPI-ESM-LR":[np.round(mpiesm_q_fut_i, 1)],
                        "Q_change_median":[q_change_median_i],
                        "Q_change_std":[q_change_std_i]})
    basins_EPSs_01.append(res)


basins_EPSs_02 = []
for basin in shp_EPSs_02.OBJECTID.unique():
    basin_i = shp_EPSs_02[shp_EPSs_02["OBJECTID"] == basin]
    #
    q_pre_i = xr_crop(shp_i=basin_i, netcdf_i=q_pre)
    q_pre_i = xr_mask(grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=q_pre_i.q, to_drop=["spatial_ref"]), netcdf_i=q_pre_i)
    q_pre_i = float(q_pre_i.q.mean())
    #
    access_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=access_q_fut)
    access_q_fut_i = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=access_q_fut_i.q, to_drop=["spatial_ref"]),
        netcdf_i=access_q_fut_i)
    access_q_fut_i = float(access_q_fut_i.q.mean())
    #
    hadgem_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=hadgem_q_fut)
    hadgem_q_fut_i = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=hadgem_q_fut_i.q, to_drop=["spatial_ref"]),
        netcdf_i=hadgem_q_fut_i)
    hadgem_q_fut_i = float(hadgem_q_fut_i.q.mean())
    #
    mpiesm_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=mpiesm_q_fut)
    mpiesm_q_fut_i = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=mpiesm_q_fut.q, to_drop=["spatial_ref"]),
        netcdf_i=mpiesm_q_fut_i)
    mpiesm_q_fut_i = float(mpiesm_q_fut_i.q.mean())

    q_change_i = [100 * (access_q_fut_i - (q_pre_i + 1)) / (q_pre_i + 1),
                  100 * (hadgem_q_fut_i - (q_pre_i + 1)) / (q_pre_i + 1),
                  100 * (mpiesm_q_fut_i - (q_pre_i + 1)) / (q_pre_i + 1)]

    q_change_median_i = np.round(np.median(q_change_i),1)
    q_change_std_i = np.round(np.std(q_change_i), 1)

    res = pd.DataFrame({"Code": basin_i.NOMB_UH_N5.tolist(),
                        "Region": basin_i.NOMB_UH_N1.tolist()[0].split(" ")[-1],
                        "Q_present":[np.round(q_pre_i, 1)],
                        "ACCESS1-0":[np.round(access_q_fut_i, 1)],
                        "HadGEM2-ES":[np.round(hadgem_q_fut_i, 1)],
                        "MPI-ESM-LR":[np.round(mpiesm_q_fut_i, 1)],
                        "Q_change_median":[q_change_median_i],
                        "Q_change_std":[q_change_std_i]})
    basins_EPSs_02.append(res)

###
basins_EPSs = pd.concat([pd.concat(basins_EPSs_01, axis=0),
                         pd.concat(basins_EPSs_02, axis=0)], axis=0)

basins_EPSs[basins_EPSs.Region == "Amazonas"].drop(["Region"], axis=1).to_csv("output/table_of_q_changes_EPSs_amazon.csv")
basins_EPSs[basins_EPSs.Region == "Pacifico"].drop(["Region"], axis=1).to_csv("output/table_of_q_changes_EPSs_pacific.csv")
basins_EPSs[basins_EPSs.Region == "Titicaca"].drop(["Region"], axis=1).to_csv("output/table_of_q_changes_EPSs_titicaca.csv")

basins_EPSs.sort_values(by="Region").to_csv("output/table_of_q_changes_EPSs.csv")