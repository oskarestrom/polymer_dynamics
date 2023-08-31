import matplotlib.pyplot as plt
import matplotlib

#Labels for the different parameters
axis_labels = {

    #Flow rate-related
    'pressure_mbar': r'Pressure drop, $\Delta p$ [mbar]',
    'u_DLD_mm_per_s': r'Flow speed, <u> [mm/s]',
    'u_DLD_m_per_s': r'$u$ [$\rm m/s$]',
    'u_DLD_um_per_s': r'$u$ [$\rm \mu m/s$]',
    'shear_rate': r'$\dot{\gamma}$ $\rm [s^{-1}]$',
    'Q_uL_per_min': r'Flow rate, $Q$ [$\rm \mu L/min$]',
    'Q_nL_per_s': r'$Q$ [$nL/s$]',
    'R': r'Hydraulic Resistance$R_h$ [$\rm Pa sm^{-3}$]',
    'R_inv': r'$1/R_h$ [$\rm Pa^{-1} s^{-1} m^{3}$]',

    #Dimensionless numbers
    'Re_DLD_25C_high_salt': '$Re$',
    'Re_low_c': '$Re$',  
    'Wi_DLD_25C_high_salt': r'$Wi$',
    'Wi_low_c': r'$Wi$',
    'El_DLD_25C_high_salt': r'$El =Wi/Re$',
    'El_low_c': r'$El =Wi/Re$',
    'De': r'$De$ ($\lambda_{Zimm}/\tau_{res}$)',
    'De_low_c': r'$De$',
    'Ma_DLD_25C_high_salt': r'$Ma = \sqrt{Re\cdot Wi}$',
    'Ma_low_c': r'$Ma = \sqrt{Re\cdot Wi}$',
    'Omega_25C_high_salt': r'$\Omega = \dot{\gamma}\cdot \eta_R$',

    #Viscosity
    'eta_25C_high_salt': r'$\eta_0$ [mPas]',
    'eta_25C_high_salt_pol': r'$\eta_{p0}$ [mPas]',
    'eta_s': r'$\eta_s$ [mPas]',
    'eta_from_R': r'$\eta_{app}$ [mPas]',
    # 'solvent_viscosity': r'$\eta_s$ [mPas]',

    #Time scale-related
    'large_scale_rel_time_25C_high_salt': r'Relaxation time, $\lambda_{\eta}$ [s]',
    'tau_zimm':r'$\tau_{Zimm} (Conc. independent)$',
    't_t': 'Transition time [s]',
    't_ratio': r'Relaxation Time / Transition time',
    'res_time_pore': 'Residence time [s]',

    'L': 'Wave Appearance',
    'I': r'Ionic Strength, $I$ [M]',    

    #Conc.
    'conc_DNA': r'$\lambda$-DNA Conc. [$\rm \mu g/mL$]',

    #Overlap conc.
    'c_overlap_ratio': r'$C/C^{*}$',
    'c_overlap_T23C': r'$C^{*}$ ($23^{\circ}C)$',
    'c_overlap_T23C_Flory': r'$C^{*}$ ($\rm \mu g/mL$)',
    'c_overlap_ideal': r'$c^*_{\theta}$',

    #Length related
    'length': 'Polymer Length [kbp]',
    'N_kbp': 'Polymer Length [kbp]',
    'N_bp': r'Number of base pairs [bp]',
    'b':'Kuhn length [m]',
    'l_p':'Persistence length [m]',
    'R_g_ideal': r'$R_g^{\theta}$ [$\rm \mu m$]', 
    'R_g_T23C': r'$R_g$[$\mu m$]$(23^{\circ}C)$', 
    'R_g_T23C_Flory': r'$R_g$ (Flory) [$\rm m$]',
    'w_eff_T23C': r'$w_{eff}$ (nm)',

    'T_ng_per_h': r'Throughput [ng/h]',
    
    #Diffusion
    'D_Z_T23C': r'$D_Z$ ($\mu m^2$)',

    #High viscosity solvent related
    'pressure_mbar_eta_s_adjusted': r'$\eta_s$-adjusted pressure drop, $\Delta p^{*}$ [mbar]',
}

def save_figure(x_col, y_col, extra_string='', dir_save='', save_as_transparent=True):
    if len(dir_save) == 0:
        dir_save = os.path.join(os.getcwd(), 'plots')
        if not os.path.exists(dir_save):
            os.mkdir(dir_save)
    file_name = f'plot_{x_col}-vs-{y_col}{extra_string}_.png'
    path_save = os.path.join(dir_save, file_name)
    plt.savefig(path_save, bbox_inches="tight", transparent=save_as_transparent) #Save figure 
    print(f'Saving figure at {path_save}')
    os.startfile(dir_save) #open file directory in windows



def split_legend_into_2_columns(ax, plot):
    """[Copied from https://stackoverflow.com/questions/56575756/how-to-split-seaborn-legend-into-multiple-columns]

    Args:
        plot ([type]): [description]
    """
    handles, labels =ax.get_legend_handles_labels()
    ax.get_legend().remove()
    plt.legend(handles, labels, ncol=2, loc='upper left', 
    bbox_to_anchor=(0.4, 1.35), frameon=False)

def replace_brackets(xlabel, ylabel):
    ylabel = ylabel.replace('[','(')
    ylabel = ylabel.replace(']',')')
    xlabel = xlabel.replace('[','(')
    xlabel = xlabel.replace(']',')')
    return xlabel, ylabel

def add_minor_ticklabels_yaxis(ax1):
    #Set minor ticklabels (from https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large)
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12) 
    ax1.yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9),numticks=12)
    ax1.yaxis.set_minor_locator(locmin)
    ax1.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

def add_minor_ticklabels_xaxis(ax1):
    #Set minor ticklabels (from https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large)
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12) 
    ax1.xaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9),numticks=12)
    ax1.xaxis.set_minor_locator(locmin)
    ax1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

def add_straight_plot_line(ax, val, axis_dir='h', linestyle='--'):
    
    # Horizontal line 2
    if axis_dir == 'h':
        ax.axhline(y=val, color="gray", linestyle=linestyle, zorder=0)
    elif axis_dir == 'v':
        ax.axvline(x=val, color="gray", linestyle=linestyle, zorder=0)
        #ax.axvline(x=9, color="gray", linestyle="-.", zorder=0)
    print(f'Adding conc {axis_dir} line at {val:.2f}')
