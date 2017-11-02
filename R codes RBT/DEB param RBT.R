##### DEB parameters for RBT

z = 4.5221        # -, zoom factor
v = 0.032453      # cm/d,  energy conductance 

E_G = 5267.5666    # J/cm^3, spec cost for structure
p_M = 343.88     # J/d.cm^3, vol-spec somatic maint

kap = 0.619      # -, allocation fraction to soma

E_Hb = 4.329e+01   # J, maturity at birth
E_Hj = 8.541e+02   # J, maturity at metam
E_Hp = 3.881e+06   # J, maturity at puberty 


k_J = 0.002        # 1/d, maturity maint rate coefficient
kap_R = 0.95       # -, reproduction efficiency

T_ref = 293.15     # K, Reference temperature
T_A = 8000   # K, Arrhenius temp

del_M = 0.10482      # 


##### Chemical potentials and densities
mu_V = 500000      # J/ mol, chemical potential of structure
mu_E = 550000      # J/ mol, chemical potential of reserve
d_V = 0.2          # g/cm^3, specific density of structure for fishes
d_E = 0.2          # g/cm^3, specific density of reserve for fishes
w_V = 23.9         # g/mol, mol-weights for (unhydrated) structure
w_E = 23.9         # g/mol, mol-weights for (unhydrated) reserve

M_V = d_V/ w_V     # mol/cm^3, volume-specific mass of structure

kap_G = mu_V * M_V / E_G     # -, growth efficiency



#### Return parameters in list
param=list(z, v, E_G, p_M, kap, p_Am, E_Hb, E_Hj, E_Hp, k_J, kap_R, kap_G, T_ref, T_A, d_E, w_E, del_M)




