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


#%%
import polymer_dynamics.dna_functions as dn

N_bp = 48500 #number of base pairs
staining_ratio = 200 #1:200 staining ratio
L = dn.get_contour_length(N_bp, staining_ratio)
print(f'{L*10**6:.1f} um')
#%%
import polymer_dynamics.dna_functions as dn
import numpy as np

N_bp = 48500 #number of base pairs
staining_ratio = 200 #1:200 staining ratio

#Ionic Strength of 5x TE:
I_5xTE_0BME = 30.993*1e-3 #Molar

#The persistence length
l_p = dn.get_l_p_OSF(I_5xTE_0BME)

#The Kuhn length
b = 2*l_p

#The number of Kuhn segments
N = dn.get_Number_of_Kuhn_Segments(staining_ratio, b=b, N_bp = N_bp)

#The ideal end-to-end distance
R_e_ideal = b * np.sqrt(N)

#The ideal radius of gyration
R_g_ideal = R_e_ideal/np.sqrt(6)

print('The persistence length is {:.1f} nm'.format(l_p*10**9))
print('The number of Kuhn segments is {:.1f}'.format(N))
print('The ideal end-to-end distance is {:.1f} um'.format(R_e_ideal*10**6))
print('The ideal radius of gyration is {:.2f} um'.format(R_g_ideal*10**6))

# %%
import polymer_dynamics.dna_functions as dn

N_bp = 48500 #number of base pairs

#Ionic Strength of 5x TE:
I_5xTE_0BME = 30.993*1e-3 #Molar

#Experiment temperature in Celsius
T_deg_C = 22 

#The persistence length
l_p = dn.get_l_p_OSF(I_5xTE_0BME)

#The Kuhn length
b = 2*l_p

#The effective width
w_eff = dn.calc_w_eff(I_5xTE_0BME, T_deg_C=T_deg_C)

#The number of Kuhn segments
N = dn.get_Number_of_Kuhn_Segments(staining_ratio, b=b, N_bp = N_bp)

#The Flory radius
R_flory = dn.calc_Flory_radius(w_eff, l_p,b,N)

print('The Florey radius is {:.2f} um'.format(R_flory*10**6))
# %%
import polymer_dynamics.dna_functions as dn
import numpy as np

N_bp = 48500 #number of base pairs
staining_ratio = 200 #1:200 staining ratio

#Ionic Strength of 5x TE:
I_5xTE_0BME = 30.993*1e-3 #Molar

#The persistence length
l_p = dn.get_l_p_OSF(I_5xTE_0BME)

#The Kuhn length
b = 2*l_p

#The number of Kuhn segments
N = dn.get_Number_of_Kuhn_Segments(staining_ratio, b=b, N_bp = N_bp)

#The ideal end-to-end distance
R_e_ideal = b * np.sqrt(N)

#The ideal radius of gyration
R_g_ideal = R_e_ideal/np.sqrt(6)

c_overlap = dn.RadiusToExpOverlapConc(N_bp, R_g_ideal)

print('The overlap concentration is {:.2f} ug/mL'.format(c_overlap))
# %%
