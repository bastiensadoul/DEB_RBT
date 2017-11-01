function [par, metaPar, txtPar] = pars_init_Oncorhynchus_mykiss(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.z = 4.5221;       free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.032453;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.61916;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 343.8835;   free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5267.556;   free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 4.329e+01; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 8.541e+02; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 3.881e+06; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 3.004e-24;  free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 10;         free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.E_Hh = 1.650e+01; free.E_Hh  = 0;   units.E_Hh = 'J';         label.E_Hh = 'maturity at hatch'; 
par.T_A = 8000;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temp'; 
par.del_M = 0.10482;  free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coef (forked length)'; 
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.f_tW1 = 1;        free.f_tW1 = 1;   units.f_tW1 = '-';        label.f_tW1 = 'scld. fctl. res.  DaviKenn2014 (first part of experiment)'; 
par.f_tW2 = 1;        free.f_tW2 = 1;   units.f_tW2 = '-';        label.f_tW2 = 'scld. fctl. res.  DaviKenn2014 (second part of experiment)'; 
par.f_tW3 = 0.36165;  free.f_tW3 = 1;   units.f_tW3 = '-';        label.f_tW3 = 'scld. fctl. res. YaniHisa2002'; 
par.f_u12 = 0.57052;  free.f_u12 = 1;   units.f_u12 = '-';        label.f_u12 = 'scld. fctl. res. univar-10 F - 3-3wk 7'; 
par.f_u13 = 0.5368;   free.f_u13 = 1;   units.f_u13 = '-';        label.f_u13 = 'scld. fctl. res. u13 - ad libitum'; 
par.f_u13_3 = 0.13792;  free.f_u13_3 = 1;   units.f_u13_3 = '-';      label.f_u13_3 = 'scld. fctl. res. u13 - 3'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 


%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
