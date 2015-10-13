function [par, metaPar, txtPar] = pars_init_Oncorhynchus_mykiss(metaData)

metaPar.model = 'abj'; % see online manual for explanation and alternatives

% reference parameter (not to be changed)
par.T_ref = C2K(20); free.T_ref = 0; units.T_ref = 'K';     label.T_ref = 'Reference temperature';

%% Core parameters
par.z = 17.07;      free.z     = 1; units.z = '-';          label.z = 'zoom factor';
par.F_m = 6.5;      free.F_m   = 0; units.F_m = 'l/d.cm^2'; label.F_m = '{F_m}, max spec searching rate';
par.kap_X = 0.8;    free.kap_X = 0; units.kap_X = '-';      label.kap_X = 'digestion efficiency of food to reserve';
par.kap_P = 0.1;    free.kap_P = 0; units.kap_P = '-';      label.kap_P = 'faecation efficiency of food to faeces';
par.v = 0.04431;    free.v     = 1; units.v = 'cm/d';       label.v = 'energy conductance';
par.kap = 0.6875;   free.kap   = 1; units.kap = '-';        label.kap = 'allocation fraction to soma';
par.kap_R = 0.95;   free.kap_R = 0; units.kap_R = '-';      label.kap_R = 'reproduction efficiency';
par.p_M = 18.9;     free.p_M   = 1; units.p_M = 'J/d.cm^3'; label.p_M = '[p_M], vol-spec somatic maint';
par.p_T = 0;        free.p_T   = 0; units.p_T = 'J/d.cm^2'; label.p_T = '{p_T}, surf-spec somatic maint';
par.k_J = 0.001;    free.k_J   = 0; units.k_J = '1/d';      label.k_J = 'maturity maint rate coefficient';
par.E_G = 5228;     free.E_G   = 1; units.E_G = 'J/cm^3';   label.E_G = '[E_G], spec cost for structure';
par.E_Hh = 7.496e0; free.E_Hh  = 1; units.E_Hh = 'J';       label.E_Hh = 'maturity at hatch';
par.E_Hb = 1.091e1; free.E_Hb  = 1; units.E_Hb = 'J';       label.E_Hb = 'maturity at birth';
par.E_Hj = 1.630e2; free.E_Hj  = 1; units.E_Hj = 'J';       label.E_Hj = 'maturity at metam';
par.E_Hp = 2.200e4; free.E_Hp  = 1; units.E_Hp = 'J';       label.E_Hp = 'maturity at puberty';
par.h_a = 1.211e-13;free.h_a   = 1; units.h_a = '1/d^2';    label.h_a = 'Weibull aging acceleration';
par.s_G = 10;       free.s_G   = 0; units.s_G = '-';        label.s_G = 'Gompertz stress coefficient';

%% Auxiliary parameters
par.T_A = 8000;     free.T_A   = 0; units.T_A = 'K';        label.T_A = 'Arrhenius temp';
par.del_M = 0.1353; free.del_M = 1; units.del_M = '-';      label.del_M = 'shape coefficient';

%% Environmental parameters (temperatures are in data)
par.f    = 1;     free.f     = 1; units.f    = '-';       label.f = 'scaled functional response for 0-var data';
par.f_tW = 1;   free.f_tW  = 1; units.f_tW = '-';       label.f_tW = 'scaled functional response for tW data';
par.f_tWL_Davidson2014 = 1;   free.f_tWL_Davidson2014  = 0; units.f_tWL_Davidson2014 = '-';       label.f_tWL_Davidson2014 = 'scaled functional response for tW data';
par.f_tW_gw150meancontrol = 1; free.f_tW_gw150meancontrol  = 0; units.f_tW_gw150meancontrol = '-';       label.f_tW_gw150meancontrol = 'scaled functional response for tW_gw150meancontrol data';
par.f_tW_gw124bvarmeancontrol = 1; free.f_tW_gw124bvarmeancontrol  = 0; units.f_tW_gw124bvarmeancontrol = '-';       label.f_tW_gw124bvarmeancontrol = 'scaled functional response for tW_gw150meancontrol data';
par.W_0 = 1.471;    free.W_0   = 0; units.W_0 = 'g';        label.W_0 = 'wet weight at t = 0 for YaniHisa2002';
par.WDavidson2014_0 = 0.1;    free.WDavidson2014_0   = 0; units.WDavidson2014_0 = 'g';        label.WDavidson2014_0 = 'wet weight at t = 0 for tWDavidson2014';
par.W150meancontrol_0 = 0.1282200;    free.W150meancontrol_0   = 0; units.W150meancontrol_0 = 'g';        label.W150meancontrol_0 = 'wet weight at t = 0 for tW150meancontrol';
par.W124bvarmeancontrol_0 = 0.1098937;    free.W124bvarmeancontrol_0   = 0; units.W124bvarmeancontrol_0 = 'g';        label.W124bvarmeancontrol_0 = 'wet weight at t = 0 for tW124bvarmeancontrol';

% par.f_Tah = 1.0;         free.f_Tah  = 0;        units.f_Tah = '-';        label.f_Tah = 'scaled functional response for Tah data';
% par.f_Tab = 1.0;         free.f_Tab  = 0;        units.f_Tab = '-';        label.f_Tab = 'scaled functional response for Tab data';

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output:
txtPar.units = units; txtPar.label = label; par.free = free; 
