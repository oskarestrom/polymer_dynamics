# -*- coding: utf-8 -*-
"""
Created on Wed May 12 14:11:26 2021

@author: Oskar Erik Ström
"""
#Import libraries:
import matplotlib.pyplot as plt #for plotting
plt.rcParams['figure.dpi'] = 200 #updating the dpi
import numpy as np
import math  #For importing the value of Pi
import os
import pandas as pd
import pickle
import seaborn as sns
import polymer_dynamics.fun_plotting as fpl
# %%  Constants

# DNA properties
L_singleBP = 3.40e-10  # Length of one base pair [m]
l_p = 50e-9 # Persistence Length at excess salt [m]
b = 2*l_p # Kuhn Length at excess salt [m]

# Weight and Molar Mass calculation
m_bp_avg_Da = 650 # Avg. mass 1 base pair [Da] or [g/mol]
M_bp_avg = m_bp_avg_Da/1000 # Molar mass 1 base pair, [kg/mol]

# Misc. constants
T_25 = 25+273.15 #Temperature = 25C in [K]
FloryExp = 0.5877 #The Flory Exponent

# Properties of water
eta_20C = 1.002e-3 #Dynamic viscosity of water at 20C [Pa*s]
rho = 997 #kg/m^3

# YOYO-1 properties
L_YOYO1 = 5.1e-10 # DNA strand length addition per intercalating YOYO-1 molecule

#Constants for calculating the effective width
k_B = 1.38064852*1e-23 # Boltzmann's constant [m^2 kg s^-2 K^-1]
T = 23+273.15 #Temperature [K]
N_A = 6.02214076*1e23 # Avogadro's number [mol^-1]
epsilon = 78.3 #Dielectric constant of water, assuming T=25C, https://nvlpubs.nist.gov/nistpubs/jres/56/jresv56n1p1_a1b.pdf
epsilon_0 = 8.85418782*1e-12 #Vacuum permittivity [m^-3 kg^-1 s^4 A^2]
e = 1.60217662 * 1e-19 #Elementary charge [C]
nu_eff = -0.593*1e-10 #[e/m] -0.593 e/Å, Appendix A, Wall effects on DNA stretch and relaxation, Stigter 2001

#Overlap Concentration of lambda-DNA at T=23C
c_overlap_lambda_T23 = 119.7 # For staining 1:200, buffer 5x TE, 3% BME, at 23C. Based on the Flory radius.

#Ionic Strengths for various buffers [M]
#with 3% BME
I_002xTE_003BME = 0.239*1e-3
I_005xTE_003BME = 0.609*1e-3
I_01xTE_003BME = 1.22*1e-3
I_05xTE_003BME = 5.668382*1e-3
I_1xTE_003BME = 0.010612
I_5xTE_003BME = 43.62*1e-3

#without BME
I_01xTE_0BME = 0.608481*1e-3
I_0125xTE_0BME = 0.000760
I_1xTE_0BME = 6.103*1e-3
I_5xTE_0BME = 30.993*1e-3

I_013xTE_0BME = 0.000791 	
I_015xTE_0BME = 0.000912
I_02xTE_0BME = 0.001216
I_5xTE_0BME = 31.0*1e-3


# %%  Functions Ionic Strength Calculation 

def calc_dyeExtensionRatio(n):
    if isinstance(n, str):
        raise ValueError(f'n should be a number and not a str: {n}')
    if n == 0:
        return 1
    else:
        return (L_singleBP + L_YOYO1/n)/L_singleBP  
     
def get_Number_of_Kuhn_Segments(staining_ratio=0, b=b, N_bp = 48502):
    L_raw_DNA = N_bp * L_singleBP # Raw Contour Length [m]
    L = L_raw_DNA*calc_dyeExtensionRatio(staining_ratio)
    # N = np.round(L/b).astype(int)
    N = L/b
    return N

# %%  Functions Ionic Strength Calculation 
def get_l_p_Dobrynin(I):
    """returns the persistence length in [m] based on the ionic strength in [M] based on the formula of Dobrynin"""
    return (46.1 + 1.9195/np.sqrt(I))*1e-9
def get_l_p_OSF(I):
    """returns the persistence length in [m] based on the ionic strength in [M] based on the formula of Odijk−Skolnick−Fixman (OSF) theory"""
    return (50 + 0.0324/I)*1e-9

# %%  Functions to calculate the Debye Length, Effective Width and Flory Radius
def calc_Debye_length(T,I):
    """[Calculate the Debye length]

    Args:
        T ([float]): [Temperature in Kelvin [K]]
        I ([float]): [Ionic Strength in [M]]

    Returns:
        Debye_length (float): [Debye length in [m]]
    """
    kappa_squared = (2000 * N_A * (e**2) * I)/(epsilon_0 * epsilon * k_B * T) #[m^-2]
    Debye_length = 1/np.sqrt(kappa_squared) # [m]
    return Debye_length

def calc_w_eff(I, T_deg_C=23):
    """[Calculation of the effective width. For the source, see Reisner, 2007, https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.058302
    ]

    Args:
        I ([float]): [Ionic Strength in [M]]
        T_deg_C (float, optional): [Temperature in degrees Celcius]. Defaults to 23.

    Returns:
        w_eff (float): [effective width in [m]]
    """
    T = T_deg_C+273.15
    debye_length = calc_Debye_length(T,I) # [m]
    kappa = 1/debye_length # [m^-1]    
    w_eff = (1/kappa)*(0.7704 + np.log(2*np.pi*nu_eff**2/(epsilon*epsilon_0*k_B*T*kappa))) # effective width [m]
    return w_eff

def calc_Flory_radius(w_eff, l_p,b,N):
    """[Calculates the End-end Flory Radius]

    Args:
        w_eff ([float]): [effective width in [m]]
        l_p ([float]): [persistence length in [m]]
        b ([float]): [Kuhn length in [m]]
        N ([float]): [Number of Kuhn segments]

    Returns:
       R_e ([float]): [End-end Flory Radius in [m]]
    """
    R_e = (w_eff * l_p)**(1/5) * (b*N)**FloryExp
    return R_e


# %%  Functions related to overlap conc. calculations
def RadiusToExpOverlapConc(L_BP_DNA, R_g):
    """Converts a DNA sample with a length [bp] and Radius of gyration [m] into the overlap concentration"""
    if R_g > 0:
        # Weight and Molar Mass calculation
        m_bp_avg_Da = 650 # Avg. mass 1 base pair [Da] or [g/mol]
        N_A = 6.02E+23 # Avogadro constant [mol^-1]	
        M_bp_avg = m_bp_avg_Da/1000 # Molar mass 1 base pair, [kg/mol]
        M = L_BP_DNA * M_bp_avg # Molar mass DNA, [kg/mol]
        c_overlap_SI = 3*M / (4*math.pi*N_A*R_g**3) # Overlap Concentration (SI-units) [kg/m^3]
        c_overlap = (1e9/1e6) * c_overlap_SI # Overlap Concentration [ug/mL]
        return c_overlap
    else:
        return 0
#Functions
def ContourLengthToOverlapConc(L_kBP_DNA):
    # Converts a contour length [kbp] to the literature (Pan 2018) overlap concentration
    
    if L_kBP_DNA == 48.5: #Lambda-DNA
        c_overlap = 43 #[Pan 2018, 25C] T4-DNA
        ref = '[Pan 2018, 25$^\circ$C]'
    elif L_kBP_DNA == 165.6:
        c_overlap = 18 #[Pan 2018, 25C] 
        ref = '[Pan 2018, 25$^\circ$C]'
    elif L_kBP_DNA == 25: 
        c_overlap = 68 #[Pan 2018, 25C]   
        ref = '25 kbp, [Pan 2018, 25$^\circ$C]'
    elif L_kBP_DNA == 2.9: 
        c_overlap = 278 #[Pan 2014, 25C]  
        ref = '2.9 kbp, [Pan 2014, 25$^\circ$C]'
    else:
        c_overlap = 0
    return c_overlap

def ContourLengthToIdealOverlapConc(L_kBP_DNA):
    # Converts a contour length [in kbp] into the ideal chain overlap concentration
    
    L_BP_DNA = L_kBP_DNA * 1e3
    c = 422 #DNA Concentration [ng/uL]

    # DNA properties
    L_singleBP = 3.40e-10  # Length of one base pair [m]
    L_p = 50e-9 # Persistence Length [m]
    b = 2*L_p # Kuhn Length [m]
    
    # DNA Length calculation
    L_raw_DNA = L_BP_DNA*L_singleBP # Contour Length [m]
    # L_DNA = L_raw_DNA*dyeExtensionRatio # DNA Length because of added dyes (extending the chain length)
    N = np.round(L_raw_DNA/b) #Number of Kuhn Segments
    
    # Weight and Molar Mass calculation
    m_bp_avg_Da = 650 # Avg. mass 1 base pair [Da] or [g/mol]
    N_A = 6.02E+23 # Avogadro constant [mol^-1]	
    M_bp_avg = m_bp_avg_Da/1000 # Molar mass 1 base pair, [kg/mol]
    M = L_BP_DNA * M_bp_avg # Molar mass DNA, [kg/mol]
    
    # End-End distance and Radius of gyration calculations (Ideal chain)
    R_e = b*np.sqrt(N) # End-end distance for long molecules with the WLC model [m]
    R_g_ideal = R_e/np.sqrt(6) # Radius of gyration, ideal chain [m]
    c_overlap_SI_ideal = 3*M / (4*math.pi*N_A*R_g_ideal**3) # Overlap Concentration (SI-units) for an ideal chain [kg/m^3]
    c_overlap_ideal = (1e9/1e6) * c_overlap_SI_ideal # Overlap Concentration for an ideal chain [ug/mL]
    return c_overlap_ideal

def calc_length_with_overlap_conc_unity(conc_DNA=422, staining_ratio = 50, salt='5x TE, 3% BME', x_overlap=1, theory='Flory'):
    """[Calculates the DNA length where the input concentration equals the overlap concentration ]

    Args:
        staining_ratio (int, optional): [description]. Defaults to 50.
        salt (str, optional): [description]. Defaults to '5x TE, 3% BME'.

    Raises:
        ValueError: [description]

    Returns:mmm
        [type]: [description]
    """
    
    if salt != '5x TE, 3% BME':
        raise ValueError('Incorrect salt')

    #Generate an array of lengths
    N_bp = np.arange(1000,100000,10)

    # DNA properties
    b = 2*get_l_p_OSF(I_5xTE_003BME) # Kuhn Length [m] 
    N = get_Number_of_Kuhn_Segments(staining_ratio, b=b, N_bp = N_bp)

    if theory == 'Prakash':
        R_e_ideal = b* np.sqrt(N)
        R_g_ideal = R_e_ideal/np.sqrt(6)
        R_g_T23C = R_g_ideal*calc_alpha(T_deg_C=23, N_bp = N_bp)
        R_g = R_g_T23C
    elif theory == 'Flory':
        I = I_5xTE_003BME
        w_eff_T23C = calc_w_eff(I, T_deg_C=23)
        R_e_T23C_Flory = calc_Flory_radius(w_eff_T23C, b/2,b,N)
        R_g_T23C_Flory = R_e_T23C_Flory/np.sqrt(6)
        R_g = R_g_T23C_Flory
    else:
        raise ValueError(f'Invalid theory input ({theory})')

    c_overlaps = np.zeros(len(N_bp))
    for i in range(len(N_bp)):
        c_overlaps[i]= RadiusToExpOverlapConc(N_bp[i], R_g[i])

    #Calculate the number of times it should multiply the overlap conc. Should it be 1C* or perhaps 4C*?
    c_overlaps = c_overlaps * x_overlap

    # Find the overlap conc that is closest to the given concentration
    inx = find_closest_value_inx(c_overlaps, conc_DNA)
    # c_overlaps[inx]

    L_overlap = N_bp[inx]/1000
    return L_overlap

def calc_I_with_overlap_conc_unity(conc_DNA=50, staining_ratio = 200, N_bp = 48502, x_overlap=1, theory='Flory'):
    """[Calculates the DNA length where the input concentration equals the overlap concentration ]

    Args:
        staining_ratio (int, optional): [description]. Defaults to 50.
        salt (str, optional): [description]. Defaults to '5x TE, 3% BME'.

    Raises:
        ValueError: [description]

    Returns:mmm
        [type]: [description]
    """
    
    if N_bp != 48502:
        raise ValueError('Incorrect N_bp')

    #Generate an array of lengths
    I_list = np.logspace(-5,-1,100)

    # DNA properties
    b = 2*get_l_p_OSF(I_list) # Kuhn Length [m] 

    #Number of Kuhn Segments
    N = get_Number_of_Kuhn_Segments(staining_ratio, b=b, N_bp = N_bp) 

    if theory == 'Prakash':
        R_e_ideal = b* np.sqrt(N)
        R_g_ideal = R_e_ideal/np.sqrt(6)
        R_g_T23C = R_g_ideal*calc_alpha(T_deg_C=23, N_bp = N_bp)
        R_g = R_g_T23C
    elif theory == 'Flory':        
        w_eff_T23C = calc_w_eff(I_list, T_deg_C=23)
        R_e_T23C_Flory = calc_Flory_radius(w_eff_T23C, b/2, b, N)
        R_g_T23C_Flory = R_e_T23C_Flory/np.sqrt(6)
        R_g = R_g_T23C_Flory
    else:
        raise ValueError(f'Invalid theory input ({theory})')

    c_overlaps = np.zeros(len(I_list))
    for i in range(len(I_list)):
        c_overlaps[i]= RadiusToExpOverlapConc(N_bp, R_g[i])

    #Calculate the number of times it should multiply the overlap conc. Should it be 1C* or perhaps 4C*?
    c_overlaps = c_overlaps * x_overlap

    # Find the overlap conc that is closest to the given concentration
    inx = find_closest_value_inx(c_overlaps, conc_DNA)

    I_overlap = I_list[inx]
    print(f'inx={inx}, I_overlap = {I_overlap:.2e} M')
    return I_overlap

def create_DNA_df(N_bp=np.array([]), 
                  C_list=np.array([]), 
                  I_list=np.array([]), 
                  staining_list = [],
                  eta_s_list=np.array([]), 
                  perform_rounding=True, 
                  eta_s = 1*1e-3, 
                  T_deg_C = 22):
    """[Creates a date frame that can be used to double-check the DNA physics theory.]

    Args:
        N_bp ([type], optional): [Decide what DNA length should be used. If set to an empty array, a pre-defined array will be used.
        ]. Defaults to np.array([0]).

        possible stainings:
        - 'raw'
        - '1:N' where N is a real number from 4 to infinity

        Possible current values of 'salt' in the 'I_list'
        '0.2x TE'
        '0.1x TE'
        '0.13x TE'
        '1x TE'
        '5x TE'
        '5x TE, 3% BME'
        'excess salt' - 1000*1e-3 #[M]



    Returns:
        [type]: [description]
    """

    #Solvent viscosity
    if len(eta_s_list) == 0:
        # If no solvent viscosities are given, use the default value for water
        eta_s_list = [eta_s] * len(N_bp)
        
    # Set Concentrations to NaN if not given
    if len(C_list) != len(N_bp):
        if len(C_list) == 0:
            C_list = [np.nan] * len(N_bp)
        elif len(C_list) == 1:
            C_list = [C_list[0]] * len(N_bp)
        else:
            raise ValueError(f'Invalid length of C_list ({len(C_list)})')

    # Create the dataframe
    df = create_the_dataframe(N_bp, C_list, I_list, staining_list, eta_s_list, perform_rounding, T_deg_C)

    #Add variables
    df = add_variables(df, T_deg_C, C_list)

    #Rounding
    if perform_rounding:
        df = round_variables(df)

    #Reset indexing and remove old indices.
    df = df.reset_index()
    df = df.drop(columns=['index'])

    print(f'Created a Data frame with length {len(df)}')  
    return df

def create_the_dataframe(N_bp, C_list, I_list, staining_list, eta_s_list, perform_rounding, T_deg_C):
   # Create the dataframes:
    df_base = pd.DataFrame(data={'N_bp':N_bp, 'conc_DNA':C_list, 'eta_s':eta_s_list}, columns=['N_bp', 'conc_DNA', 'eta_s'])
    df_list = []

    # Create multiple dataframes for different staining ratios
    for s in staining_list:
        staining_ratio = int(str.split(s, ':')[1])
        df = df_base.copy()
        df['L'] = df.apply(lambda row:row.N_bp * L_singleBP * calc_dyeExtensionRatio(staining_ratio), axis=1)
        df['N'] = df.apply(lambda row: row.L/b, axis=1)
        df['staining'] = s
        df['staining_ratio'] = staining_ratio
        df_list.append(df)

    #Combine the data frames based on the different stainings
    df1 = pd.concat(df_list)

    #Create multiple dataframes for different ionic strengths
    df_list = []
    for I in I_list:
        df_temp = df1.copy()
        df_temp['I'] = I
        df_list.append(df_temp)
    
    #Combine the data frames based on the different ionic strengths
    df = pd.concat(df_list)
    return df

def round_variables(df):
    df['R_e_ideal'] = np.round(df['R_e_ideal']*1e6,2)
    df['R_g_ideal'] = np.round(df['R_g_ideal']*1e6,2)
    df['R_g_T23C'] = np.round(df['R_g_T23C']*1e6,2)
    df['R_g_T25C'] = np.round(df['R_g_T25C']*1e6,2)
    df['R_e_T23C_Flory'] = np.round(df['R_e_T23C_Flory']*1e6,2)
    df['R_g_T23C_Flory'] = np.round(df['R_g_T23C_Flory']*1e6,2)
    # df['R_g_Flory_solvent_factor_T23C'] = np.round(df['R_g_Flory_solvent_factor_T23C']*1e6,2)
    df['w_eff_T23C'] = np.round(df['w_eff_T23C']*1e9,2)
    df['c_overlap_ideal'] = np.round(df['c_overlap_ideal'],1)
    df['c_overlap_T23C'] = np.round(df['c_overlap_T23C'],1)
    df['c_overlap_T25C'] = np.round(df['c_overlap_T25C'],1)
    df['c_overlap_T23C_Flory'] = np.round(df['c_overlap_T23C_Flory'],1)
    df['L'] = np.round(df['L']*1e6,2) 
    df['b'] = np.round(df['b']*1e6,3) 
    df['l_p'] = np.round(df['l_p']*1e9,1) 
    df['alpha_g_T23C'] = np.round(df['alpha_g_T23C'],2) 
    return df

def add_variables(df, T_deg_C, C_list):

    #Length
    df['l_p'] = df.apply(lambda row: get_l_p_OSF(row.I), axis=1) #Persistence length
    df['b'] = df.apply(lambda row: row.l_p*2, axis=1)  #Kuhn length
    df["N"] = df.apply(lambda row: get_Number_of_Kuhn_Segments(row.staining_ratio, b=row.b, N_bp = row.N_bp), axis=1)

    #Temperature
    df['T_deg_C'] = T_deg_C

    #Radii of gyration
    df['R_e_ideal'] = df.apply(lambda row: row.b* np.sqrt(row.N), axis=1)
    df['R_g_ideal'] = df.apply(lambda row: row.R_e_ideal/np.sqrt(6), axis=1)   

    #Effective width
    df['w_eff_T23C'] = df.apply(lambda row: calc_w_eff(row.I, T_deg_C=23), axis=1)

    #Debye Length
    df['debye_length'] = df.apply(lambda row: calc_Debye_length(T=23+273.15,I=row.I), axis=1)
    
    #Radius of Gyration
    df['R_e_T23C_Flory'] = df.apply(lambda row: calc_Flory_radius(row.w_eff_T23C, row.l_p,row.b,row.N), axis=1)
    df['R_g_T23C_Flory'] = df.apply(lambda row: row.R_e_T23C_Flory/np.sqrt(6), axis=1)

    #Overlap concentrations
    df['c_overlap_ideal'] = df.apply(lambda row: RadiusToExpOverlapConc(row.N_bp, row.R_g_ideal), axis=1)
    df['c_overlap_T23C_Flory'] = df.apply(lambda row: RadiusToExpOverlapConc(row.N_bp, row.R_g_T23C_Flory), axis=1)

    if len(C_list) > 0:
        df['C_ratio_Flory'] = df.apply(lambda row: row.conc_DNA/row.c_overlap_T23C_Flory, axis=1)  

    # Relaxation times
    df['tau_zimm_low_c'] = df.apply(lambda row: calc_zimm_relaxation_time(row.b, row.N, eta_s=row.eta_s,  T_deg_C=row.T_deg_C), axis=1)

    #Zimm Diffusion Coefficient
    df['D_Z'] = df.apply(lambda row: calc_D_Z_from_R_g(row.R_g_T23C_Flory, eta_s = row.eta_s, T_degC = row.T_deg_C), axis=1)

    return df

def rename_columns_df(df):
    df_renamed = df.rename(columns={'l_p': r'$l_p$ [nm]',
                            'R_g_ideal': r'$R_g^{\theta}$ [$\mu m$]', 
                            'R_g_T23C': r'$R_g$[$\mu m$]$(23^{\circ}C)$', 
                            'c_overlap_ideal': r'$c^*_{\theta}$', 
                            'c_overlap_T23C': r'$c^*$ ($23^{\circ}C$)'})
    return df_renamed

# %%  misc. functions
def find_closest_value_inx(arr, val):
    """[Finds the index of the closest value of val in an array]

    Args:
        arr ([np.array]): [Array]
        val ([float]): [value]

    Returns:
        index ([int]): [index of the closest value of val in an array]
    """
    difference_array = np.absolute(arr-val)
    index = difference_array.argmin()
    return index

def plot_DNA_df(df, x_col, y_col, kwargs,save_fig = False, add_legend = False, y_logscale=False, x_logscale=False, split_legend=False):
    
    #Main plot
    f, ax = plt.subplots()
    
    extra_string = '' #Add note to the file name if saved
    g = sns.scatterplot(x=x_col, y=y_col, data=df, **kwargs)
    
    if y_logscale:
        ax.set(yscale="log") #Make the x-scale logarithmic 
        fpl.add_minor_ticklabels_yaxis(ax) 
        extra_string = extra_string+'_ylog'

    if x_logscale:
        ax.set(xscale="log") #Make the x-scale logarithmic
        fpl.add_minor_ticklabels_xaxis(ax)
        extra_string = extra_string+'_xlog'
        
    ylabel = fpl.axis_labels[y_col]
    xlabel = fpl.axis_labels[x_col]
    xlabel, ylabel = fpl.replace_brackets(xlabel, ylabel) 
    plt.ylabel(ylabel) 
    plt.xlabel(xlabel)  

    if y_col == 'c_overlap_ideal' or y_col == 'c_overlap_T23C':  
        c=422
        fpl.add_straight_plot_line(ax, c, axis_dir='h')

        L_overlap = calc_length_with_overlap_conc_unity(staining_ratio = 50, salt='5x TE, 3% BME')*1000 # Lambda-DNA Overlap Concentration for an ideal chain [ug/mL]
        print(f'DNA length where the input concentration equals the overlap conc: {L_overlap} bp')
        #Calculate the DNA lengths where the input concentration equals the entanglement conc. limits
        L_lower_lim_entanglement = calc_length_with_overlap_conc_unity(staining_ratio = 50, salt='5x TE, 3% BME',x_overlap=3)*1000
        L_lower_high_entanglement =  calc_length_with_overlap_conc_unity(staining_ratio = 50, salt='5x TE, 3% BME',x_overlap=10)*1000
        print(f'DNA length where the input concentration equals the lower entanglement conc: {L_lower_lim_entanglement} bp')
        fpl.add_straight_plot_line(ax, val=L_overlap, axis_dir='v')  
        fpl.add_straight_plot_line(ax, val=L_lower_lim_entanglement, axis_dir='v')  
    if add_legend:
        if split_legend:
            fpl.split_legend_into_2_columns(ax, g)
        else:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left') #Set the legend out of the graph area
    else:
        plt.legend([],[], frameon=False)
        
    if save_fig: #Save figure
        fpl.save_figure(x_col, y_col, extra_string)
    plt.show()

def calc_zimm_relaxation_time(b, N, eta_s=1e-3,  T_deg_C=23):
    T_Kelvin = T_deg_C + 273.15
    tau_z = (eta_s * b**3 * N**(3*FloryExp)) / (k_B*T_Kelvin)
    return tau_z

def get_rho(solvent):
    if solvent == 'water':
        rho = 1000 #kg/m3
    elif 'sucrose' in solvent:
        if solvent.split(' ')[1] == '40%':
            rho = 1191 #kg/m3 #Based on Sucrose
        else:
            raise ValueError(f'{solvent} is not a valid input for the solvent')
    else:
        raise ValueError(f'{solvent} is not a valid input for the solvent')
    return rho

def calc_low_C_De_based_on_relaxation_time_U_and_w(u_m_per_s, tau_zimm, device_type = 'Q'):
    """[Returns the estimated Weissenberg number in the device. 
    Input parameters is:
        flow velocity [m/s]
        Concentration [ng/uL]
    Note that the array dimensions have to be set in the function
    ]

    Args:
        u ([type]): [description]

    Returns:
        [type]: [description]
    """
    if (device_type == 'Q') or (device_type == 'R'):
        L = 11*1e-6 #Pitch
    # elif device_type == 'H':
    #     L = 7*1e-6
    # elif device_type == 'C':
    #     L = 7*1e-6 #Radius of the post
    elif device_type == 'DLD':
        L = 17.5*1e-6 #[um] Pitch = D + G = 14.5 um + 2.8 um
    else:
        raise ValueError('Device type not calculated for yet...')

    shear_rate = u_m_per_s / L

    De =  shear_rate * tau_zimm
    if De == 0:
        De = np.nan 
    return De

def calc_D_Z_from_R_g(R_g, eta_s = 1e-3, T_degC = 23, prefactor = 0.196/np.sqrt(6)):
    """Calculate the Zimm Diffusion coefficient
    R_g - radius of gyration [m]    
    eta_s - solvent viscosity [Pa s]
    T_degC - Temperature [Degrees C]
    prefactor - prefactor, e.g. for good solvent conditions
    """    
    T_K = T_degC + 273.15 #Temperature [K]    
    D_Z = prefactor * k_B * T_K / (eta_s * R_g)
    return D_Z

def calc_diffusion_coefficient(R, eta_s = 1e-3, T_degC = 23):
    T_K = T_degC + 273.15 #Temperature [K]    
    D =  k_B * T_K / (eta_s * R)
    return D

