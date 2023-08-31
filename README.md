# polymer_dynamics
This package can perform simple polymer dynamics calculations.

Based on the number of base pairs for a given DNA strand, we can extract the following:
- Raw Contour Length, $L$
- Contour length dependent on salt and staining levels, $L$
- End-end Distance, $R_e$
- Ideal radius of Gyration, $R_g^{\theta}$ 
- Debye length, $1/\kappa$ or $\lambda_D$
- Effective width of a DNA molecule, $w_{eff}$
- Flory radius, $R_{\nu}$
- Estimated radius of gyration, $R_g$

In the theory section below, the equations and assumptions for these variables are given.

## Using the script
The function create_DNA_df() will create a pandas dataframe containing all the variables listed above.
```python
import polymer_dynamics.dna_functions as dn
import numpy as np
import pandas as pd

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


```
|conc_DNA|L    |I      |l_p|w_eff_T23C|R_e_T23C_Flory|R_g_T23C_Flory|c_overlap_T23C_Flory|tau_zimm_low_c|C_ratio_Flory|
|--------|-----|-------|---|----------|--------------|--------------|--------------------|--------------|-------------|
|400     |3.426e-07|0.03099|5.105e-08|5.79e-09  |1.243e-07     |5.075e-08     |1973                |0.002207      |0.2028       |
|400     |6.851e-07|0.03099|5.105e-08|5.79e-09  |1.868e-07     |7.626e-08     |1162                |0.00749       |0.3442       |
|400     |1.713e-06|0.03099|5.105e-08|5.79e-09  |3.201e-07     |1.307e-07     |577.6               |0.03768       |0.6925       |
|400     |3.426e-06|0.03099|5.105e-08|5.79e-09  |4.81e-07      |1.964e-07     |340.3               |0.1279        |1.175        |
|400     |6.851e-06|0.03099|5.105e-08|5.79e-09  |7.229e-07     |2.951e-07     |200.5               |0.4341        |1.995        |
|400     |1.661e-05|0.03099|5.105e-08|5.79e-09  |1.217e-06     |4.967e-07     |102                 |2.07          |3.921        |
|400     |5.673e-05|0.03099|5.105e-08|5.79e-09  |2.504e-06     |1.022e-06     |39.96               |18.04         |10.01        |


## Theory

### Raw contour length
The length of a single base pair, $L_{single\space bp} = 34.0$ nm. The raw contour length is given by multiplying the length of a single base pair with the number of base pairs:

$$L_{raw} = N_{bp}\cdot L_{single\space bp}$$


## Freely-jointed chain model
In the freely-jointed chain (FJC) model, a polymer is seen as a chain of N stiff, rod-like segments, with no limitations in the bond angles. With this model, DNA is seen as a random-walk where the segments have equal probability to fluctuate in all directions. The segment length is called the Kuhn length, b, and is two times the persistence length:
$$b = 2l_p$$

The contour length, L, is the length at full extension and is given by the product of the Kuhn length and the number of Kuhn segments: 

$$L = b N$$

The number of Kuhn segments, $N$ for a raw DNA strand is then calculated as $N = L_{raw}/b$.

## Persistence length

We assume the persistence length of raw DNA to be 50 nm [Dorfman, 2013].

#### Ionic strength
The ionic strength, $I$, is calculated by summing the product of the concentration and squared charge of all ions in the solution:

$$I = \dfrac{1}{2} \sum_{i=1}^{n} c_i z_i^2$$

where the one half is added to include both anions and cations, $c_i$ is the molar concentration of the ion i (in M) and $z_i$ is the charge of the ion i. 

### Persistence length increase with decreasing salt level
The persistence length, $l_p$ changes with the ionic strength as the following (according to Odijk−Skolnick−Fixman (OSF) theory (see Beyondseq, Dorfman)): 

$$l_p = l_p^{'} + \dfrac{0.0324 M}{ I} nm$$

where $ l_p^{'} \approx 50$  nm is the bare persistence length.

An alternative, later model is given by Dobrynin. The persistence length, $l_p$ changes with the ionic strength as the following (according to Dobrynin (see Beyondseq, Dorfman)): 

$$l_p = 46.1 + \dfrac{1.9195 M}{ \sqrt{I}} nm$$

In this work, however, OSF-theory is used to model the length-increase by ionic strength.

Reisner et al., (2007) showed how the length of DNA-strands doubled when going a medium with high salt to low salt. They used the OSF-theory to explain their results [[Reisner, 2007]](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.058302). 

#### Staining effect on the Persistence length 
Kundukad et al. (Doyle group)showed that the persistence length remains unaffected with the addition of YOYO-1 intercalating dyes [Kundakad 2014]. 

### Effect of staining ratio on the contour length
The contour length is increased with 0.51±0.14 nm for every incorporated dye molecule [Günther 2010]. For a staining ratio of 1 dye molecule per 10 base pairs, we get a 15% increase in the contour length. For a ratio of 1:4 we get a 37.5% increase.

Kundukad et al. (Doyle group) showed in 2014 using atomic force microscopy measurements that the persistence length of DNA does not change when YOYO-1 intercalates. Their measurements indicates that at the maximum staining ratio of 1 YOYO per 4 base pairs, the length of DNA increases by 38%. This is equivalent to 0.5 nm increase per bound YOYO molecule.

$$L_{1:4} = 1.38\cdot L_{raw}$$

The results of Kundukad et al. [Kundukad, 2014] contend the claim by Reisner et al. [[Reisner, 2007]](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.058302) where they argued that the persistence length should increase with the same factor as $L$. However, Reisner et al. provided no evidence for their claim.

The number of Kuhn segments, $N$, which is used to define the radius, is then calculated as $N = L/b = L_{raw}*s_f/b$, where $s_f$ is the staining factor.

#### Table: Number of kuhn segments, N, for various staining ratios at excess salt
df_N = df[df['salt']=='excess'].copy()
df_N=pd.pivot_table(df,index=['N_bp'],columns='staining',values='N')
df_N = df_N[['raw', '1:200', '1:50', '1:10', '1:4']]
df_N.head(20)

### Flory Radius
### Effect of ionic environment on the radius
It is important to note that most experiments conducted on the polymer physics of DNA use excess salt to screen the electrostatic interactions. 

When the other effects are also included, the following expression for the end-end distance is given:

$$R_e \approx (w_{eff} l_p)^{1/5}(bN)^{\nu} = (w_{eff} l_p)^{1/5}L^{\nu}$$

where $w_{eff}$ is the effective width of the polymer. Reisner et al. described the $w_{eff}$ for strongly charged chains as [[Reisner, 2007]](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.058302):

$$w_{eff} = \dfrac{1}{\kappa}[0.7704+log(\dfrac{\nu^2_{eff}}{2\epsilon\epsilon_0 k_B T \kappa})]$$

where $1/\kappa$ is the Debye length (also sometimes denoted $\lambda_D$), $\epsilon$ the dielectric constant of water, $\epsilon_0$ the permittivity of free space, $\nu_{eff}$ is an effective DNA line charge density, $k_B$ is Boltzman's constant and $T$ the temperature. Note that [[Hsieh, 2008 (Doyle group)]](http://web.mit.edu/doylegroup/pubs/Nanolett_hsieh_08.pdf) used a constant of $2\pi$ instead of $0.5$ inside the log-expression. If I do this, I can replicate their values of $w_{eff}$.

$\kappa$  in turn is defined as:

$$\kappa^2  = \dfrac{2000N_Ae^2I}{\epsilon_0\epsilon k_BT}$$

where $N_A$ is Avogadro's number, $e$ the electronic charge and $I$ the ionic strength of the buffer. The Debye length should be on the order of tens of nanometer for low ionic strengths (~1 mM) and Å for higher ionic strenth buffers (see calculated values at [[Montes, 2019]](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6932854/)).

Iarko et al. calculated the effective width of DNA at 5 x TBE to be 4.6 nm. [[Iarko, 2015]](https://publications.lib.chalmers.se/records/fulltext/225267/local_225267.pdf). This value is the same as I end up with in the calculation below.