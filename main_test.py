#%%
#Public modules
import matplotlib.pyplot as plt
import math
import pandas as pd
import numpy as np

#Personal modules
import polymer_dynamics.dna_functions as dn

#Reload modules
import importlib

importlib.reload(dn)
C_list = np.array([400])
n = len(C_list)
lengths = np.array([48502]*n)
df = dn.create_DNA_df(N_bp=lengths, C_list=C_list, perform_rounding=False)

df_lambda = df[(df['staining'] == '1:200')]
df_lambda = df_lambda.reset_index()

pd.set_option('display.float_format', lambda x: '%.4g' % x)
df_lambda[['conc_DNA','L','I', 'l_p','w_eff_T23C','R_e_T23C_Flory','R_g_T23C_Flory','c_overlap_T23C_Flory', 'tau_zimm_low_c','C_ratio_Flory']]
#%%

importlib.reload(dn)
C_list = np.array([400])
n = len(C_list)
lengths = np.array([48502]*n)
I_list = np.array([30.993*1e-3]) #I_5xTE_0BME
df = dn.create_DNA_df(N_bp=lengths, C_list=C_list, I_list=I_list)

df_lambda = df[(df['staining'] == '1:200')]
df_lambda = df_lambda.reset_index()
# indices = [0,1,2,3,4,5,22,10,28]
# df_lambda = df_lambda.iloc[indices]

pd.set_option('display.float_format', lambda x: '%.4g' % x)
df_lambda[['conc_DNA','L','I', 'l_p','w_eff_T23C','R_e_T23C_Flory','R_g_T23C_Flory','c_overlap_T23C_Flory', 'tau_zimm_low_c','C_ratio_Flory']]

#%%