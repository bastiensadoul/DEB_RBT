%% core DEB model parameters from "the blank"
% sets (initial values for) parameters

function [par, txtPar] = pars_blank(metaData) 

%% reference temp: don't change!
par.T_ref = 293.15; free.T_ref = 0; units.T_ref = 'K';      label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.z = 7.256;      free.z = 0;     units.z = '-';          label.z = 'zoom factor'; 
par.F_m = 6.5;      free.F_m = 0;   units.F_m = 'l/d.cm^2'; label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;    free.kap_X = 0; units.kap_X = '-';      label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;    free.kap_P = 0; units.kap_P = '-';      label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.0379;     free.v = 0;     units.v = 'cm/d';       label.v = 'energy conductance'; 
par.kap = 0.6935;   free.kap = 0;   units.kap = '-';        label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;   free.kap_R = 0; units.kap_R = '-';      label.kap_R = 'reproduction efficiency'; 
par.p_M = 133.8;    free.p_M = 0;   units.p_M = 'J/d.cm^3'; label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;        free.p_T = 0;   units.p_T = 'J/d.cm^2'; label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.001991; free.k_J = 0;   units.k_J = '1/d';      label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5025;     free.E_G = 0;   units.E_G = 'J/cm^3';   label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 5.295e1; free.E_Hb = 0;  units.E_Hb = 'J';       label.E_Hb = 'maturity at birth'; 
par.E_Hj = 4.210e2; free.E_Hj = 0;  units.E_Hj = 'J';       label.E_Hj = 'maturity at metam'; 
par.E_Hp = 1.708e6; free.E_Hp = 0;  units.E_Hp = 'J';       label.E_Hp = 'maturity at puberty'; 
par.h_a = 3.195e-47;free.h_a = 0;   units.h_a = '1/d^2';    label.h_a = 'Weibull aging acceleration'; 
par.s_G = 10;       free.s_G = 0;   units.s_G = '-';        label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.E_Hh = 26.55;   free.E_Hh = 0;  units.E_Hh = 'J';       label.E_Hh = 'maturity at hatch'; 
par.T_A = 6930;     free.T_A = 0;   units.T_A = 'K';        label.T_A = 'Arrhenius temp'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
