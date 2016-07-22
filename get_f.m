clear all; clc; close all;

% metaData needed to automatically generate chemical parameters from
% pars_blank
metaData.phylum     = 'Chordata';   
metaData.class      = 'Actinopterygii'; 

[data, auxData, txtData, weights] = mydata_BPA0;

[par, txtPar] = pars_blank(metaData);

% initialise f value
par.f = 1; 
par.free.f = 1;

% parameters of the estimation
estim_options('default'); % runs estimation, uses nmregr method and filter
estim_options('max_step_number',18);  % set options for parameter estimation
estim_options('max_fun_evals',5e3);    % set options for parameter estimation
estim_options('lossfunction', 'E'); % there are three possibilities: 'E', 'F' or 'I'

% filter: 0 for hold, 1 for pass
% flag: indicator of reason for not passing the filter
filterNm = 'filter_abj';
[filter, flag] = feval(filterNm,par);
if ~filter
    fprintf('The seed parameter set is not realistic. \n');
      print_filterflag(flag);  
end

newData.tWw = data.tW_gw150A;
newWeights.tWw = weights.tW_gw150A;
[par, info, nsteps] = petregr_f('predict_tWw', par, newData, auxData, newWeights, filterNm); % WLS estimate parameters using overwrite

newPar = par;
prdData = predict_tWw(newPar, newData, auxData);

plot(newData.tWw(:,1), prdData.tWw, 'r', newData.tWw(:,1), newData.tWw(:,2),'mo')

