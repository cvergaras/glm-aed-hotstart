#%%
from glmpy import plots
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import os
#%%

# lake = plots.LakePlotter("output/lake.csv")
# fig, ax = plt.subplots(figsize=(10, 5))
# # lake.lake_level(ax=ax)
Data_GLM = nc.Dataset("nc/output_8.nc")
#%%
ph = Data_GLM.variables["temp"][:, :, 0, 0]
print(ph[0])

#%%

# lake = plots.LakePlotter("output/lake.csv")
# fig, ax = plt.subplots(figsize=(10, 5))
# # lake.lake_level(ax=ax)
# Data_GLM = nc.Dataset("nc/output_2.nc")
# ph = Data_GLM.variables["temp"][:, :, 0, 0]
# print(ph[8])

# #%%
# fig, ax = plt.subplots(figsize=(10, 5))
# nc1 = plots.NCProfile("output/output.nc")
# out = nc1.plot_var(ax=ax, var="temp")#, min_diff=0.01)
# col_bar = fig.colorbar(out)
# col_bar.set_label("Salinity (mg/L)")
# plt.show()
# #%%
# # #%%
# # fig, ax = plt.subplots(figsize=(10, 5))
# # nc1 = plots.NCProfile("output/output.nc")
# # out = nc1.plot_var(ax=ax, var="PHQ_PHQ_pe", min_diff=0.0001)
# # col_bar = fig.colorbar(out)
# # col_bar.set_label("pe")
# # plt.show()

# # #%%
# # fig, ax = plt.subplots(figsize=(10, 5))
# # nc1 = plots.NCProfile("output/output.nc")
# # out = nc1.plot_var(ax=ax, var="PHQ_KI_Pyrite")#, min_diff=0.01)
# # col_bar = fig.colorbar(out)
# # col_bar.set_label("Pyrite (mmol/m3)")
# # plt.show()

# # #%%
# # fig, ax = plt.subplots(figsize=(10, 5))
# # nc1 = plots.NCProfile("output/output.nc")
# # out = nc1.plot_var(ax=ax, var="OXY_oxy")#,max=400)#, min_diff=0.01)
# # col_bar = fig.colorbar(out)
# # col_bar.set_label("DO (mmol/m3)")
# # plt.show()

# # #%%
# # fig, ax = plt.subplots(figsize=(10, 5))
# # nc1 = plots.NCProfile("output/output.nc")
# # out = nc1.plot_var(ax=ax, var="temp")#,max=400)#, min_diff=0.01)
# # col_bar = fig.colorbar(out)
# # col_bar.set_label("Temperature (°C)")
# # plt.show()

# # #%%
# # fig, ax = plt.subplots(figsize=(10, 5))
# # nc1 = plots.NCProfile("output/output.nc")
# # out = nc1.plot_var(ax=ax, var="PHQ_EP_Calcite", min_diff=10)
# # col_bar = fig.colorbar(out)
# # col_bar.set_label("Temperature (°C)")
# # plt.show()

# # # %%
# # nc1 = plots.NCProfile("output/output.nc")
# # vars = nc1.get_vars()
# # print(vars)

# # # %%
# # Data_GLM = nc.Dataset("output/output.nc")
# # ph = Data_GLM.variables["PHQ_preDO"][:, :, 0, 0]
# # for i in range(0,1000,100):
# #     print(i)
# #     print(ph[i].max())
# # #%%
# # Data_GLM = nc.Dataset("output/output.nc")
# # ph1 = Data_GLM.variables["OXY_oxy"][:, :, 0, 0]
# # print(ph1[200*24+1])

# # #%%
# # # Data_GLM = nc.Dataset("output/output.nc")
# # ph = Data_GLM.variables["PHQ_preDO"][:, :, 0, 0]
# # print(ph[3])
# # # print(ph.min(), ph.max())
# # #%%
# # ph = Data_GLM.variables["PHQ_postDO"][:, :, 0, 0]
# # print(ph[200])

# # # %%
