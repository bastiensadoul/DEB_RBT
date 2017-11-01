%% run_funique_control
% This script estimates functional response for each control dataset
% the results are saved in prdData_f.mat files
% The first seven datasets in mydata_BPA  are control
% see inside the code for the format in which the predictions and the f
% values are saved.
% prdData.mat can be used as input for the R script to study the way RE
% varies with time

clear all; clc; close all;

global study % the initial egg sized depends on whether it is the 150 or the 124 experiments

% metaData needed to automatically generate chemical parameters from
% pars_init_Oncorhynchus_mykiss - which are the parameters from the blank
% taken from AmP
metaData.phylum     = 'Chordata';   
metaData.class      = 'Actinopterygii'; 
[data, auxData, txtData, weights] = mydata_BPA;
[par, metaPar, txtPar] = pars_init_funiqueall(metaData);
load f_prdData_funique_all % load f estimated for all of the tanks
f_150 = mean([f.tW_gw150A, f.tW_gw150B, f.tW_gw150C]); % mean value for the 150 experiments
f_124 = mean([f.tW_gw124iniA, f.tW_gw124iniB, f.tW_gw124iniC]); % mean value for the 150 experiments
f_124fin = f.tW_gw124fin; % value for the 124fin
 
   
 % Get names and number of tanks
[nm, nst] = fieldnmnst_st(data); % cell array of string with fieldnames of data
  
% add auxilliary data 
newAuxData.temp = auxData.temp;
newAuxData.pMoA = 'control'; % choose physiological mode of action
       
    for j = 1:7 % replace with nst if you want to run with all data - otherwise put 7 to only estimate for control
        if (strfind(nm{j}, '150') > 0)
            study = 'e150';   
            par.f = f_150;
        else
            study = 'e124';
            par.f = f_124;
        end 
        newData.tW = data.(nm{j});   
        % creates prdData for each tank: 1st col = time
        %                                2nd col = predict
        %                                3rd col = real data
        prdData.(nm{j}) = newData.tW(:,1);
        [pD, info] = predict_tW(par, newData, newAuxData);
        prdData.(nm{j})(:,2) = pD.tW;
        prdData.(nm{j})(:,3) = newData.tW(:,2);

%         figure()
%         plot(newData.tW(:,1), prdData.(nm{j})(:,2), 'r', 'linewidth',2); hold on;
%         plot(newData.tW(:,1), newData.tW(:,2),'bo','markerfacecolor','b')
%         xlabel('time since fertilisation, d'); ylabel('wet weight, g')
%         title(nm{j})
%         set(gca,'Fontsize',12); 
        
       % creates prdData for each tank: 1st col = time
        %                                2nd col = predict
        %                                3rd col = real data
        prdData.(nm{j}) = newData.tW(:,1);
        a = predict_tW(par, newData, newAuxData);
        prdData.(nm{j})(:,2) = a.tW;
        prdData.(nm{j})(:,3) = newData.tW(:,2);

    end
   save('prdData_funique_control', 'prdData');       
 
          
 % Relative difference between data and model as function of time
 % red: 150, cyan 124, dashed lines are exposed individuals
figure() 

for j = 1:nst % replace with nst if you want to run with all data - otherwise I put 7 here because the first four fields are for the controls
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
legend(nm)
xlabel('d, time since fertilization')
ylabel('-, relative difference between model and data using mean control f')
set(gca,'Fontsize',12);

% Bargraph for f by condition?

 