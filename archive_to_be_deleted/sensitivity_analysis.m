clear all; clc; close all;

% metaData needed to automatically generate chemical parameters from
% pars_blank
metaData.phylum     = 'Chordata';   
metaData.class      = 'Actinopterygii'; 

%[data, auxData, txtData, weights] = mydata_BPA0;
%[data, auxData, txtData, weights] = mydata_BPA03to300;

auxData.temp.tWw = C2K(12);
auxData.t0.tWw = 'dpf';

newData.tWw= [linspace(70,300, 200)
linspace(70,300, 200)]';

[par, txtPar] = pars_blank(metaData);

% initialise f value
par.f = 1; 
par.free.f = 1;

% parameters of the estimation
estim_options('default'); % runs estimation, uses nmregr method and filter
estim_options('max_step_number',100);  % set options for parameter estimation
estim_options('max_fun_evals',5e3);    % set options for parameter estimation
estim_options('lossfunction', 'I'); % there are three possibilities: 'E', 'F' or 'I'

estim = 0;

% filter: 0 for hold, 1 for pass
% flag: indicator of reason for not passing the filter
filterNm = 'filter_abj';
[filter, flag] = feval(filterNm,par);
if ~filter
    fprintf('The seed parameter set is not realistic. \n');
      print_filterflag(flag);  
end

   newAuxData.temp = auxData.temp;

   
n=5;
var='p_M';
vari=linspace(25,140,n);   
% 
% var='v';
% vari=linspace(0.01,0.06,n);   
% 
% var='E_G';
% vari=linspace(4500,13000,n); 
% 
% var='kap';
% vari=linspace(0.5,0.9,n); 

% var='f';
% vari=linspace(0.5,1,n); 

% var='E_Hb';
% vari=linspace(11,100,n); 

% var='E_Hj';
% vari=linspace(300,600,n); 



% Initialize the cell containing the text. For each "n" there is a cell.
LegendString = cell(1,n);

for i = 1:n
    
    eval(strcat('par.', var , '=vari(i)'));
    newAuxData.t0.tWw = auxData.t0.tWw;
    
    prdData = predict_tWw(par, newData, newAuxData);
    
    if i>1
        hold on
    end
    plot(newData.tWw(:,1), prdData.tWw);
    LegendString{i} = strcat(var,' = ', num2str(vari(i)));
    
%     f.(nm{j})=newPar.f;
   clear prdData
   
end

legend(LegendString, 'Location', 'northwest')
