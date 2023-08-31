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

## Installation

### Dependencies
You need to install the following libraries to run this package:
```shell
- pip install numpy > 1.21.5
- pip install matplotlib > 3.5.2
- pip install pandas >= 1.3.3, < 2
- pip install pickle
- pip install seaborn
```
### Install this script
Find this package on [test.pypi.org](https://test.pypi.org/project/polymer-dynamics)
```shell
pip install -i https://test.pypi.org/simple/ polymer-dynamics
```

## Using the script
The function create_DNA_df() will create a pandas dataframe containing all the variables listed above. The DNA lengths 1,2,5,10,20,48.5 and 165.6 kbp are slected as long as an ionic strength corresponding to 5x TE with a dye to base pair staining ratio at 1:200. The temperature is set to 22 C (which the effective width and Zimm relaxation time depends on).
```python
import polymer_dynamics.dna_functions as dn
import numpy as np
import pandas as pd

#Ionic Strength of 5x TE:
I_5xTE_0BME = 30.993*1e-3 #Molar

#Experiment temperature in Celsius
T_deg_C = 22 

#Concentrations list
C_list = np.array([400]) #DNA concentations in ug/mL, useful if one wants the C/C_overlap ratio    

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
df[['N_bp','L','I', 'l_p','w_eff_T23C','R_e_ideal','R_g_ideal','R_e_T23C_Flory','R_g_T23C_Flory','c_overlap_ideal','c_overlap_T23C_Flory', 'tau_zimm_low_c','C_ratio_Flory']]
```
The following table is generated:

|N_bp|L    |I      |l_p|w_eff_T23C|R_e_ideal|R_g_ideal|R_e_T23C_Flory|R_g_T23C_Flory|c_overlap_ideal|c_overlap_T23C_Flory|tau_zimm_low_c|C_ratio_Flory|
|----|-----|-------|---|----------|---------|---------|--------------|--------------|---------------|--------------------|--------------|-------------|
|1000|3.426e-07|0.03099|5.105e-08|5.79e-09  |1.87e-07 |7.634e-08|1.243e-07     |5.075e-08     |579.3          |1973                |0.002207      |0.2028       |
|2000|6.851e-07|0.03099|5.105e-08|5.79e-09  |2.645e-07|1.08e-07 |1.868e-07     |7.626e-08     |409.6          |1162                |0.00749       |0.3442       |
|5000|1.713e-06|0.03099|5.105e-08|5.79e-09  |4.182e-07|1.707e-07|3.201e-07     |1.307e-07     |259.1          |577.6               |0.03768       |0.6925       |
|10000|3.426e-06|0.03099|5.105e-08|5.79e-09  |5.914e-07|2.414e-07|4.81e-07      |1.964e-07     |183.2          |340.3               |0.1279        |1.175        |
|20000|6.851e-06|0.03099|5.105e-08|5.79e-09  |8.363e-07|3.414e-07|7.229e-07     |2.951e-07     |129.5          |200.5               |0.4341        |1.995        |
|48502|1.661e-05|0.03099|5.105e-08|5.79e-09  |1.302e-06|5.317e-07|1.217e-06     |4.967e-07     |83.18          |102                 |2.07          |3.921        |
|165600|5.673e-05|0.03099|5.105e-08|5.79e-09  |2.406e-06|9.824e-07|2.504e-06     |1.022e-06     |45.02          |39.96               |18.04         |10.01        |



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

### Calculating the Contour Length with this script
We can then use the function get_contour_length() to extract the contour length for any DNA length and staining ratio
```python
import polymer_dynamics.dna_functions as dn

N_bp = 48500 #number of base pairs
staining_ratio = 200 #1:200 staining ratio
L = dn.get_contour_length(N_bp, staining_ratio)
print(f'{L*10**6:.1f} um')
```
#### Output:
```
16.6 um
```

### Ideal Chain Radius Calculation
The worm-like chain (WLC) model (also called the Kratky-Porod model) describes the DNA dynamics more more accurately compared to the FJC-model by defining the polymer as semi-flexible. It means that the chain is stiff at the length scale of a monomer but flexible at the length scale of the entire polymer. The angle of a single segment affects the angle of its neighbour compared to FJC where the angles of segments are completely unrelated. The WLC model can be said to describe a polymer as a continuously flexible isotropic rod. The end-end distance of the polymer, $R_e$, is a measure of the size of the coiled-up polymer. The WLC model describes $R_e$ for long chains ($L≫L_p$) with [Colby & Rubinstein, 2003]:

$$R_e=b \cdot \sqrt{N} = b \cdot \sqrt{L/b}$$

Another way to characterize the size of a polymer is with the radius of gyration, $R_G$. The radius of gyration of a polymer is the root mean square distance from its center of mass. For an ideal chain, it is given by [Colby & Rubinstein, 2003]:

$$R_G^{\theta}=\dfrac{R_e}{\sqrt{6}}=b\sqrt{\dfrac{N}{6}}$$

$$ => $$

$$R_G^{\theta} \propto \sqrt{L}$$

$$R_G^{\theta} \propto l_p\sqrt{L}$$

### Calculating the ideal DNA Radii with this script

```python
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
```
#### Output:
```
The persistence length is 51.0 nm
The number of Kuhn segments is 162.7
The ideal end-to-end distance is 1.3 um
The ideal radius of gyration is 0.53 um
```

### Flory Radius and the effect of ionic environment on the radius
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

### Calculating the Flory radius with this script
```python
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

print('The Flory radius is {:.2f} um'.format(R_flory*10**6))
```
#### Output:
```
The Flory radius is 1.22 um
```
### The overlap concentration
The relation between the volume fraction of a polymer solution, i.e. the ratio of the volume of the polymer in the solution to the volume of the solution is related to the mass concentration of the polymer, the ratio of the polymer mass to the volyme of the solution is given by:

$$\phi= \dfrac{c}{\rho}= c \dfrac{v_{mon}N_A}{M_{mon}}$$

where $\rho = \dfrac{M_{mon}}{v_{mon}N_A}$

The volume fraction of **a single molecule inside its pervaded volume** (pervaded volume = the volume spanned by the polymer chain, $V=R^3$ [Colby & Rubinstein, 2003]]) is called the overlap volume fraction $\phi^*$:

$$\phi^* = \dfrac{\text{Volume taken up by the monomers}}{\text{Volume taken up by the polymer blob}}=\dfrac{Nv_{mon}}{V}$$

where $N$ is the number of monomers, $v_{mon}$ is the monomer volume and $V$ is the pervaded volume of a single molecule chain. The corresponding concentration at this volume fraction is valled the overlap concentration, $c^*$:

$$c^*=\dfrac{\text{Base pair concentration of the monomers}}{\text{Base pair concentration of the polymer blob}}\dfrac{\rho Nv_{mon}}{V}=\dfrac{M}{V N_{A}}$$

where $M$ is the polymer molecular weight, $N_A$ is Avogadro’s constant. For a polymer with a volume, $V=4πR_g^3/3 $, the overlap concentration is written as [Doi 1996]:

$$c^{*}= \dfrac{M}{\dfrac{4π}{3} R_g^3 N_A}$$

where $M$ is the polymer molecular weight, $N_A$ is Avogadro’s constant, and $R_g$ is the radius of gyration. Using this equation, $c^*$ equals 64 $\mu g/mL$ for lambda DNA (48.5 kbp) and using $R_g$=0.58 $\mu m$.

$M$ is calculated by multiplying the number of base pairs to the average weight of a base pair (650 Da or 0.65 kg/mol):

$$M = N_{bp} \cdot M_{bp \ avg} = N_{bp} \cdot 0.65 \ kg/mol $$

In other words, the concentration when the solution concentration is equal to that of the pervaded volume of a single polymer.
$$c^*=\dfrac{\text{C solution}}{\text{C inside a blob}}$$

### Calculating the Overlap Concentration with this script
```python
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
```
#### Output:
```
The overlap concentration is 83.18 uM
```
### References

Colby, R. H., & Rubinstein, M. (2003). Polymer physics. New-York: Oxford University, 100, 274.

Doi, M. (1996). Introduction to polymer physics. Oxford university press.

Hsieh, C. C., Balducci, A., & Doyle, P. S. (2008). Ionic effects on the equilibrium dynamics of DNA confined in nanoslits. Nano letters, 8(6), 1683-1688.

Iarko, V., Werner, E., Nyberg, L. K., Müller, V., Fritzsche, J., Ambjörnsson, T., ... & Mehlig, B. (2015). Extension of nanoconfined DNA: Quantitative comparison between experiment and theory. Physical review E, 92(6), 062701.

Kundukad, B., Yan, J., & Doyle, P. S. (2014). Effect of YOYO-1 on the mechanical properties of DNA. Soft matter, 10(48), 9721-9728.

Montes, R. J., Ladd, A. J., & Butler, J. E. (2019). Transverse migration and microfluidic concentration of DNA using Newtonian buffers. Biomicrofluidics, 13(4).

Reisner, W., Beech, J. P., Larsen, N. B., Flyvbjerg, H., Kristensen, A., & Tegenfeldt, J. O. (2007). Nanoconfinement-enhanced conformational response of single DNA molecules to changes in ionic environment. Physical review letters, 99(5), 058302.

