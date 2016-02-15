%% pars_init_Oncorhynchus_mykiss_BPA0
% sets (initial values for) parameters

function [par, metaPar, txtPar] = pars_init_Oncorhynchus_mykiss_BPA0(metaData)

metaPar.model = 'abj'; 

%% core primary parameters 
par.z = 4.3938;  free.z = 0;  units.z = '-';  label.z = 'zoom factor'; 
par.F_m = 6.5;  free.F_m = 0;  units.F_m = 'l/d.cm^2';  label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;  free.kap_X = 0;  units.kap_X = '-';  label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;  free.kap_P = 0;  units.kap_P = '-';  label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.030864;  free.v = 0;  units.v = 'cm/d';  label.v = 'energy conductance'; 
par.kap = 0.4416;  free.kap = 0;  units.kap = '-';  label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;  free.kap_R = 0;  units.kap_R = '-';  label.kap_R = 'reproduction efficiency'; 
par.p_M = 173.2605;  free.p_M = 0;  units.p_M = 'J/d.cm^3';  label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;  free.p_T = 0;  units.p_T = 'J/d.cm^2';  label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.00013194;  free.k_J = 0;  units.k_J = '1/d';  label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5025.245;  free.E_G = 0;  units.E_G = 'J/cm^3';  label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 106.1359;  free.E_Hb = 0;  units.E_Hb = 'J';  label.E_Hb = 'maturity at birth'; 
par.E_Hj = 2771.7214;  free.E_Hj = 0;  units.E_Hj = 'J';  label.E_Hj = 'maturity at metam'; 
par.E_Hp = 5796157.9753;  free.E_Hp = 0;  units.E_Hp = 'J';  label.E_Hp = 'maturity at puberty'; 

%% other parameters 
par.T_A = 6237.8303;  free.T_A = 0;  units.T_A = 'K';  label.T_A = 'Arrhenius temp'; 
par.T_ref = 293.15;  free.T_ref = 0;  units.T_ref = 'K';  label.T_ref = 'Reference temperature'; 
par.del_M = 0.11347;  free.del_M = 0;  units.del_M = '-';  label.del_M = 'shape coefficient DaviKlemm2014'; 
par.f_tW_gw124fin = 0.73906;  free.f_tW_gw124fin = 0;  units.f_tW_gw124fin = '-';  label.f_tW_gw124fin = 'scaled functional response for gw124fin'; 
par.f_tW_gw124iniA = 0.65928;  free.f_tW_gw124iniA = 0;  units.f_tW_gw124iniA = '-';  label.f_tW_gw124iniA = 'scaled functional response for gw124ini tank A'; 
par.f_tW_gw124iniB = 0.61494;  free.f_tW_gw124iniB = 0;  units.f_tW_gw124iniB = '-';  label.f_tW_gw124iniB = 'scaled functional response for gw124ini tank B'; 
par.f_tW_gw124iniC = 0.65547;  free.f_tW_gw124iniC = 0;  units.f_tW_gw124iniC = '-';  label.f_tW_gw124iniC = 'scaled functional response for gw124ini tank C'; 
par.f_tW_gw150A = 0.54634;  free.f_tW_gw150A = 0;  units.f_tW_gw150A = '-';  label.f_tW_gw150A = 'scaled functional response for gw150 tank A'; 
par.f_tW_gw150B = 0.55039;  free.f_tW_gw150B = 0;  units.f_tW_gw150B = '-';  label.f_tW_gw150B = 'scaled functional response for gw150 tank B'; 
par.f_tW_gw150C = 0.55059;  free.f_tW_gw150C = 1;  units.f_tW_gw150C = '-';  label.f_tW_gw150C = 'scaled functional response for gw150 tank C'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
