import geopandas as gpd
import pandas as pd
import xarray as xr
import rioxarray
import numpy as np

exec(open("src/crop_mask.py").read())

shp_EPSs = gpd.read_file("data/raw/observed/EPSs/C_aporte.shp")
access_q_fut = xr.open_dataset("data/processed/future/runoff/ACCESS/median_ensemble_runoff.nc")
access_q_fut0 = xr.open_dataset("/home/adrian/Downloads/median_ensemble_runoff_ACCESS1-0_rcp85_2036-2065 (1).nc") # from figshare [just to check]
hadgem_q_fut = xr.open_dataset("data/processed/future/runoff/HADGEM/median_ensemble_runoff.nc")
hadgem_q_fut0 = xr.open_dataset("/home/adrian/Downloads/median_ensemble_runoff_HadGEM2-ES_rcp85_2036-2065.nc") # from figshare
mpiesm_q_fut = xr.open_dataset("data/processed/future/runoff/MPIESM/median_ensemble_runoff.nc")
mpiesm_q_fut0 = xr.open_dataset("/home/adrian/Downloads/median_ensemble_runoff_MPI-ESM-LR_rcp85_2036-2065.nc") # from figshare
q_pre = xr.open_dataset("data/processed/present/PISCO/runoff/median_ensemble_runoff.nc")
q_pre0 = xr.open_dataset("/home/adrian/Downloads/median_ensemble_runoff_1982-2011.nc") # from figshare

basins_EPSs = []
for basin in shp_EPSs.nomeps.unique():
    basin_i = shp_EPSs[shp_EPSs["nomeps"] == basin]
    #
    q_pre_i = xr_crop(shp_i=basin_i, netcdf_i=q_pre)
    q_pre_i_g = xr_mask(grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=q_pre_i.q, to_drop=["spatial_ref"]), netcdf_i=q_pre_i)
    q_pre_mean_i = np.round(float(q_pre_i_g.q.mean()), 1)
    q_pre_median_i = np.round(float(q_pre_i_g.q.median()), 1)
    #
    access_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=access_q_fut)
    access_q_fut_i_g = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=access_q_fut_i.q, to_drop=["spatial_ref"]),
        netcdf_i=access_q_fut_i)
    access_q_fut_mean_i = np.round(float(access_q_fut_i_g.q.mean()), 1)
    access_q_fut_median_i = np.round(float(access_q_fut_i_g.q.median()), 1)
    #
    hadgem_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=hadgem_q_fut)
    hadgem_q_fut_i_g = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=hadgem_q_fut_i.q, to_drop=["spatial_ref"]),
        netcdf_i=hadgem_q_fut_i)
    hadgem_q_fut_mean_i = np.round(float(hadgem_q_fut_i_g.q.mean()), 1)
    hadgem_q_fut_median_i = np.round(float(hadgem_q_fut_i_g.q.median()), 1)
    #
    mpiesm_q_fut_i = xr_crop(shp_i=basin_i, netcdf_i=mpiesm_q_fut)
    mpiesm_q_fut_i_g = xr_mask(
        grid_mask=xr_shp_to_grid(shp_i=basin_i, netcdf_array=mpiesm_q_fut.q, to_drop=["spatial_ref"]),
        netcdf_i=mpiesm_q_fut_i)
    mpiesm_q_fut_mean_i = np.round(float(mpiesm_q_fut_i_g.q.mean()), 1)
    mpiesm_q_fut_median_i = np.round(float(mpiesm_q_fut_i_g.q.median()), 1)

    q_change_i = [100 * (access_q_fut_i_g - (q_pre_i_g + 1)) / (q_pre_i_g + 1),
                  100 * (hadgem_q_fut_i_g - (q_pre_i_g + 1)) / (q_pre_i_g + 1),
                  100 * (mpiesm_q_fut_i_g - (q_pre_i_g + 1)) / (q_pre_i_g + 1)]

    q_change_mean_i = [np.round(float(i.q.mean()), 1) for i in q_change_i]
    q_change_median_i = [np.round(float(i.q.median()), 1) for i in q_change_i]

    res_present0 = pd.DataFrame({"EPS":[basin], "Información":["Presente"], "Variable":["Q"], "Promedio":[q_pre_mean_i], "Mediana":[q_pre_median_i]})
    res_present1 = pd.DataFrame({"EPS":[basin], "Información":["ACCESS 1.0"], "Variable":["Q"], "Promedio":[access_q_fut_mean_i], "Mediana":[access_q_fut_median_i]})
    res_present2 = pd.DataFrame({"EPS":[basin], "Información":["HadGEM2-ES"], "Variable":["Q"], "Promedio":[hadgem_q_fut_mean_i], "Mediana":[hadgem_q_fut_median_i]})
    res_present3 = pd.DataFrame({"EPS":[basin], "Información":["MPI-ESM-LR"], "Variable":["Q"], "Promedio":[mpiesm_q_fut_mean_i], "Mediana":[mpiesm_q_fut_median_i]})
    res_present4 = pd.DataFrame({"EPS":[basin], "Información":["ACCESS 1.0"], "Variable":["ΔQ/Q"], "Promedio":[q_change_mean_i[0]], "Mediana":[q_change_median_i[0]]})
    res_present5 = pd.DataFrame({"EPS":[basin], "Información":["HadGEM2-ES"], "Variable":["ΔQ/Q"], "Promedio":[q_change_mean_i[1]], "Mediana":[q_change_median_i[1]]})
    res_present6 = pd.DataFrame({"EPS":[basin], "Información":["MPI-ESM-LR"], "Variable":["ΔQ/Q"], "Promedio":[q_change_mean_i[2]], "Mediana":[q_change_median_i[2]]})

    res = pd.concat([res_present0, res_present1, res_present2, res_present3, res_present4, res_present5, res_present6], axis=0)
    basins_EPSs.append(res)

table_of_res = pd.concat(basins_EPSs, axis=0)
table_of_res = table_of_res.sort_values("EPS")
table_of_res.to_csv("output/Tabla_cambios_en_Q_por_EPS.csv", index=False)

def check_sing_values(x):
   x_sign = [i > 0 for i in x]
   if np.sum(x_sign) == 3:
       return "Incremento"
   elif np.sum(x_sign) == 0:
       return "Decrecimiento"
   else:
       return "No coincidencia"

Promedio_val = table_of_res[table_of_res["Variable"] == "ΔQ/Q"].groupby("EPS").apply(lambda grp: check_sing_values(x=grp["Promedio"]))
Mediana_val = table_of_res[table_of_res["Variable"] == "ΔQ/Q"].groupby("EPS").apply(lambda grp: check_sing_values(x=grp["Mediana"]))

changes_summary = pd.concat([Promedio_val, Mediana_val], 1)
changes_summary.columns = ["Promedio", "Mediana"]
changes_summary["EPS"] = changes_summary.index
changes_summary.to_csv("output/Resumen_cambios_en_Q_por_EPS.csv", index=False)