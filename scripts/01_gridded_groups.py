import xarray as xr
import numpy as np
import geopandas as gpd
import rioxarray

exec(open("src/crop_mask.py").read())

shp_obs = gpd.read_file("data/processed/present/PISCO/runoff/Q_mov30yearly_shp.shp")
pisco_grid = xr.open_dataset("data/processed/present/PISCO/prec/P_mov30yearly.nc").p.isel(time=0)
pisco_grid = pisco_grid.rio.set_crs("epsg:4326")
shp_SA = gpd.read_file("data/raw/shps/Sudamérica.shp").to_crs({"init": "epsg:4326"})

shp_order = sorted(shp_obs.Region.unique())
for i in range(len(shp_order)):
    shp_obs.loc[(shp_obs.Region == shp_order[i]), "Region"] = i + 1

gridded_groups = []
for group in shp_obs.Region.unique():
    subx = shp_obs.loc[shp_obs.Region == group]
    grid_value = pisco_grid.rio.clip(subx.geometry, drop=False)
    grid_value.values[~np.isnan(grid_value.values)] = int(subx.Region.unique())
    grid_value = grid_value.fillna(0)
    gridded_groups.append(grid_value)

gridded_groups = (gridded_groups[0] + gridded_groups[1] + gridded_groups[2] + gridded_groups[3] +
                  gridded_groups[4] + gridded_groups[5] + gridded_groups[6] + gridded_groups[7] +
                  gridded_groups[8] + gridded_groups[9] + gridded_groups[10])

gridded_groups.values[gridded_groups.values == 0] = np.nan
gridded_groups = gridded_groups.rio.write_nodata(np.nan)
gridded_groups = gridded_groups.rio.interpolate_na(method='nearest')
shp_SA = xr_shp_to_grid(shp_i=shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_array=xr_crop(shp_i = shp_SA[shp_SA["PAÍS"] == "Perú"], netcdf_i=gridded_groups))
gridded_groups = xr_mask(grid_mask=shp_SA, netcdf_i=gridded_groups).reindex_like(pisco_grid)
gridded_groups = gridded_groups.to_dataset(name="gridded_groups").drop(["spatial_ref", "time"])
gridded_groups.to_netcdf("data/processed/others/budyko_groups.nc")