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

#Ionic Strength of 5x TE:
I_5xTE_0BME = 30.993*1e-3 #Molar

#Experiment temperature in Celsius
T_deg_C = 22 

#Concentrations list
C_list = np.array([400]) #DNA concentations in uM, useful if one wants the C/C_overlap ratio    

#Lengths list
N_bp = np.array([1000, 2000, 5000, 10000, 20000,  48502, 165600])

#Ionic Strengths list
I_list = np.array([I_5xTE_0BME])

#Create dataframe
df = dn.create_DNA_df(N_bp=N_bp, 
                      C_list=C_list, 
                      I_list = I_list,
                      staining_list = ['1:200'],
                      perform_rounding=False,
                      T_deg_C=T_deg_C)

df = df.reset_index()

pd.set_option('display.float_format', lambda x: '%.4g' % x)
df[['conc_DNA','L','I', 'l_p','w_eff_T23C','R_g','R_e_T23C_Flory','R_g_T23C_Flory','c_overlap_T23C_Flory', 'tau_zimm_low_c','C_ratio_Flory']]
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

N_bp = np.array([1000, 1200, 1500, 2000, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 10000, 15000, 17000, 20000, 40000, 48502, 165600])