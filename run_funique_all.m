%% run_funique_all
% Starrlight & Bastien May 289 2017
% This script estimates functional response for each dataset
% the results are saved in prdData_funique_all.mat and
% f_prdData_funique_all.mat
% The first seven datasets in mydata_BPA  are control
% we use petregr_f and estimate f for each data set separately taking all
% of the other parameters as given.
% see inside the code for the format in which the predictions and the f
% values are saved.
% prdData.mat can be used as input for the R script to study the way RE
% varies with time
% RE ïs computed as function of time and plotted 

clear all; clc; close all;

global study % the initial egg sized depends on whether it is the 150 or the 124 experiments

% metaData needed to automatically generate chemical parameters from
% pars_init_Oncorhynchus_mykiss - which are the parameters from the blank
% taken from AmP
metaData.phylum     = 'Chordata';   
metaData.class      = 'Actinopterygii'; 
[data, auxData, txtData, weights] = mydata_BPA;
[par, metaPar, txtPar] = pars_init_funiqueall(metaData);

% initialise f value for each tank
par.f = 0.7; 
par.free.f = 1;

% parameters of the estimation
estim_options('default'); % runs estimation, uses nmregr method and filter
estim_options('max_step_number',1e3);  % set options for parameter estimation
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

 % Get names and number of tanks
[nm, nst] = fieldnmnst_st(data); % cell array of string with fieldnames of data

newAuxData.temp = auxData.temp;
newAuxData.pMoA = 'control'; % choose physiological mode of action
       
    for j = 1:7 % replace with nst if you want to run with all data - otherwise put 7 to do controls only
        newData.tW = data.(nm{j});
        newWeights.tW = weights.(nm{j});
       
        if (strfind(nm{j}, '150') > 0)
            study = 'e150';             
        else
            study = 'e124';
        end
        
        if estim
        [newPar, info, nsteps] = petregr_f('predict_tW', par, newData, newAuxData, newWeights, filterNm); % WLS estimate parameters using overwrite
        [merr, rerr, prdInfo] = mre_st('predict_tW', newPar, newData, newAuxData, newWeights);
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

        save('prdData_funique_all', 'prdData');       
        save('f_prdData_funique_all', 'f');

        
 % Relative difference between data and model as function of time
 % red: 150, cyan 124, dashed lines are exposed individuals

figure()

for j = 1:7 % replace with nst if you want to run with all data - otherwise put 7 to do only controls
        if (strfind(nm{j}, '150') > 0)
            color = 'red';
             if (strfind(nm{j}, 'BPA') > 0)
                 linestyle = '--';
             else
                 linestyle = '-';
             end
        else
            color = 'cyan';
            if (strfind(nm{j}, 'BPA') > 0)
                 linestyle = '--';
                 else
                 linestyle = '-';
            end
        end
    EWw = prdData.(nm{j})(:,2);
    diff= (EWw-data.(nm{j})(:,2))./data.(nm{j})(:,2)*100;
    plot(data.(nm{j})(:,1), diff,'color',color,'linewidth',2,'linestyle',linestyle)
    hold on    
end
legend(nm{1:7})
xlabel('d, time since fertilization')
ylabel('-, relative difference between model and data')
set(gca,'Fontsize',12);

% Bargraph for f by condition?

 