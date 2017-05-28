function [par, metaPar, txtPar] = pars_init_reconstruction(metaData)

global tpM ;

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.z = 4.456;        free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.085867;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.52204;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 362.8475;   free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5201.7976;  free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 9.467e+01; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 6.988e+02; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 5.186e+06; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 2.902e-24;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 10;         free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.E_Hh = 1.288e+01; free.E_Hh  = 1;   units.E_Hh = 'J';         label.E_Hh = 'maturity at hatch'; 
par.T_A = 9037.3025;  free.T_A   = 1;   units.T_A = 'K';          label.T_A = 'Arrhenius temp'; 
par.del_M = 0.1002;   free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coef (forked length)'; 
par.del_M2 = 0.077925;  free.del_M2 = 1;   units.del_M2 = '-';       label.del_M2 = 'shape coef, LW  '; 
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 

%% And then add at the end of the pars_init_pMreconstruction:
%% reconstructed pM values for each tank
% this code automatically generates knot-coordinates for each knot abscissa
% collected in tpM.
% tpM is a 1-n vector which is passes as a global and which defined in run_reconstruction

[data, auxData, metaData, txtData, weights] = mydata_pMreconstruction();

for tank = 1:length(fieldnames(data))
    tanknames=fieldnames(data);
    tankname=char(tanknames(tank));
    for i = 1:length(tpM)
       fieldName     = ['pM_', tankname, '_', num2str(tpM(i))];
    % ------- MODIFY starting values here -------------------------------------   
        par.(fieldName)   = 363; % same starting value for all
        % -------------------------------------------------------------------------   
        free.(fieldName)  = 1;  % all the knot coordinates are estimated as free parameters
        units.(fieldName) = 'J/V/d';            
        label.(fieldName) = ['pM for ', 'pM_', tankname, '_', num2str(tpM(i))]; 
    end
end
    
%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
