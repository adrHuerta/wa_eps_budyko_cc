import geopandas as gpd
import pandas as pd
import numpy as np

hydro_time = pd.date_range("1982-01-01", "2016-12-31", freq="M")

pisco_runoff = gpd.read_file("data/raw/observed/PISCO/runoff/PEQ_1981-2020_subcuencas_GR2M_Peru_SHP/Subbasins_GR2M_Peru.shp")
pisco_runoff_values = pd.read_excel("data/raw/observed/PISCO/runoff/PEQ_1981-2020_subcuencas_GR2M_Peru.xlsx",
                                    sheet_name = "Q_mm", engine = 'openpyxl')

pisco_runoff_values_dates = pisco_runoff_values["Fecha"]
pisco_runoff_values = pisco_runoff_values.drop(["Fecha"], axis = 1)
pisco_runoff_values = pisco_runoff_values.set_index(pisco_runoff_values_dates)
pisco_runoff_values = pisco_runoff_values.loc["1981-09-01":"2016-08-31"]
pisco_runoff_values = pisco_runoff_values.set_index(hydro_time)
pisco_runoff_values = pisco_runoff_values.resample("1Y").apply(lambda x: np.sum(x))
pisco_runoff["q_values_exp"] = pisco_runoff_values.loc["2000-12-31"].values
pisco_runoff = pisco_runoff.sort_values(by="Region")

# detecting bad basins
pisco_runoff[(pisco_runoff["q_values_exp"] < 50) & (pisco_runoff["Region"].isin(['A','E','G','H','I']))].plot("q_values_exp", cmap="viridis_r")
pisco_runoff[(pisco_runoff["q_values_exp"] < 500) & (pisco_runoff["Region"].isin(['D','B']))].plot("q_values_exp", cmap="viridis_r")
pisco_runoff[(pisco_runoff["q_values_exp"] < 10) & (pisco_runoff["Region"].isin(['C']))].plot("q_values_exp", cmap="viridis_r")

ID_bad_basins = [pisco_runoff[(pisco_runoff["q_values_exp"] < 50) & (pisco_runoff["Region"].isin(['A','E','G','H','I']))].GR2M_ID.values.tolist(),
                 pisco_runoff[(pisco_runoff["q_values_exp"] < 500) & (pisco_runoff["Region"].isin(['D','B']))].GR2M_ID.values.tolist(),
                 pisco_runoff[(pisco_runoff["q_values_exp"] < 10) & (pisco_runoff["Region"].isin(['C']))].GR2M_ID.values.tolist()]
ID_bad_basins = [item for sublist in ID_bad_basins for item in sublist]
ID_bad_basins = ["GR2M_ID_" + str(item) for item in ID_bad_basins]

# making runoff data from ID_bad_basin in np.nan
new_pisco_runoff_values = pisco_runoff_values
new_pisco_runoff_values[ID_bad_basins] = np.nan
new_pisco_runoff_values = new_pisco_runoff_values.rolling(window=15, center=True).mean()
new_pisco_runoff_values = np.round(new_pisco_runoff_values.dropna(axis=0, thresh=1), 2)
new_pisco_runoff_values.to_csv("data/processed/present/PISCO/runoff/Q_mov15yearly.csv")

# merging groups
new_pisco_runoff = gpd.read_file("data/raw/observed/PISCO/runoff/PEQ_1981-2020_subcuencas_GR2M_Peru_SHP/Subbasins_GR2M_Peru.shp")
new_pisco_runoff["Region"] = new_pisco_runoff["Region"].apply(lambda x: "BD" if (x == "D") or (x == "B") else x)
new_pisco_runoff["Region"] = new_pisco_runoff["Region"].apply(lambda x: "AE" if (x == "A") or (x == "E") else x)
new_pisco_runoff["Region"] = new_pisco_runoff["Region"].apply(lambda x: "GH" if (x == "G") or (x == "H") else x)
new_pisco_runoff.to_file("data/processed/present/PISCO/runoff/Q_mov15yearly_shp.shp")