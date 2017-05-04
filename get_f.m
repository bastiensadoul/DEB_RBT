%% This script estimates functional response for each control dataset 
% the results are saved in .mat files

clear all; clc; close all;

% metaData needed to automatically generate chemical parameters from
% pars_init_Oncorhynchus_mykiss - which are the parameters from the blank
% taken from AmP
metaData.phylum     = 'Chordata';   
metaData.class      = 'Actinopterygii'; 
[data, auxData, txtData, weights] = mydata_BPA;
[par, metaPar, txtPar] = pars_init_Oncorhynchus_mykiss(metaData);

% initialise f value
par.f = 1; 
par.free.f = 1;

% parameters of the estimation
estim_options('default'); % runs estimation, uses nmregr method and filter
estim_options('max_step_number',100);  % set options for parameter estimation
estim_options('max_fun_evals',5e3);    % set options for parameter estimation
estim = 1;

% filter: 0 for hold, 1 for pass
% flag: indicator of reason for not passing the filter
filterNm = 'filter_abj';
[filter, flag] = feval(filterNm,par);
    if ~filter
        fprintf('The seed parameter set is not realistic. \n');
          print_filterflag(flag);  
    end

 
[nm, nst] = fieldnmnst_st(data); % cell array of string with fieldnames of data

newAuxData.temp = auxData.temp;
newAuxData.pMoA = 'control'; % choose physiological mode of action
       
    for j = 1:7 % replace with nst if you want to run with all data - otherwise I put 7 here because the first four fields are for the controls

        newData.tW = data.(nm{j});
        newWeights.tW = weights.(nm{j});

        if estim
        [newPar, info, nsteps] = petregr_f('predict_tW', par, newData, newAuxData, newWeights, filterNm); % WLS estimate parameters using overwrite
        else
        newPar = par;
        end

        % creates prdData for each tank: 1st col = time
        %                                2nd col = predict
        %                                3rd col = real data
        prdData.(nm{j}) = newData.tW(:,1);
        a=predict_tW(newPar, newData, newAuxData);
        prdData.(nm{j})(:,2) = a.tW;
        prdData.(nm{j})(:,3) = newData.tW(:,2);

        fprintf(['f for ',nm{j},' is: %2.3f \n'],newPar.f);

        figure()
        plot(newData.tW(:,1), prdData.(nm{j})(:,2), 'r', 'linewidth',2); hold on;
        plot(newData.tW(:,1), newData.tW(:,2),'bo','markerfacecolor','b')
        xlabel('time since fertilisation, d'); ylabel('wet weight, g')
        title(nm{j})
        set(gca,'Fontsize',12); 
        
        f.(nm{j})=newPar.f;

    end

        save('prdData', 'prdData');       
        save('f_prdData', 'f');
