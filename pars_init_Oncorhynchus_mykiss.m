%% pars_init_Oncorhynchus_mykiss
% sets (initial values for) parameters

function [par, metaPar, txtPar] = pars_init_Oncorhynchus_mykiss(metaData)

metaPar.model = 'abj'; 

%% core primary parameters 
par.z = 5.0134;  free.z = 1;  units.z = '-';  label.z = 'zoom factor'; 
par.F_m = 6.5;  free.F_m = 0;  units.F_m = 'l/d.cm^2';  label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;  free.kap_X = 0;  units.kap_X = '-';  label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;  free.kap_P = 0;  units.kap_P = '-';  label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.066299;  free.v = 1;  units.v = 'cm/d';  label.v = 'energy conductance'; 
par.kap = 0.47073;  free.kap = 1;  units.kap = '-';  label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;  free.kap_R = 0;  units.kap_R = '-';  label.kap_R = 'reproduction efficiency'; 
par.p_M = 50.1853;  free.p_M = 1;  units.p_M = 'J/d.cm^3';  label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;  free.p_T = 0;  units.p_T = 'J/d.cm^2';  label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.001;  free.k_J = 0;  units.k_J = '1/d';  label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5080.8366;  free.E_G = 1;  units.E_G = 'J/cm^3';  label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 332.5223 * 0.8;  free.E_Hb = 1;  units.E_Hb = 'J';  label.E_Hb = 'maturity at birth'; 
par.E_Hj = 30297.7007 * 0.3;  free.E_Hj = 1;  units.E_Hj = 'J';  label.E_Hj = 'maturity at metam'; 
par.E_Hp = 41135.0663 * 3;  free.E_Hp = 1;  units.E_Hp = 'J';  label.E_Hp = 'maturity at puberty'; 
par.h_a = 5.7305e-72;  free.h_a = 1;  units.h_a = '1/d^2';  label.h_a = 'Weibull aging acceleration'; 
par.s_G = 10;  free.s_G = 0;  units.s_G = '-';  label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.E_Hh = 245.1778*0.8;  free.E_Hh = 1;  units.E_Hh = 'J';  label.E_Hh = 'maturity at hatch'; 
par.T_A = 8000;  free.T_A = 1;  units.T_A = 'K';  label.T_A = 'Arrhenius temp'; 
par.T_ref = 293.15;  free.T_ref = 0;  units.T_ref = 'K';  label.T_ref = 'Reference temperature'; 
par.del_M = 0.17912;  free.del_M = 1;  units.del_M = '-';  label.del_M = 'shape coefficient'; 
par.f = 1;  free.f = 0;  units.f = '-';  label.f = 'scaled functional response for 0-var data'; 
par.f_LW = 1;  free.f_LW = 0;  units.f_LW = '-';  label.f_LW = 'scaled functional response for LW data'; 
par.f_tW = 0.56211;  free.f_tW = 1;  units.f_tW = '-';  label.f_tW = 'scaled functional response for YaniHisa2002'; 
par.f_tWL = 1;  free.f_tWL = 1;  units.f_tWL = '-';  label.f_tWL = 'scaled functional response for DaviKenn2014'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
