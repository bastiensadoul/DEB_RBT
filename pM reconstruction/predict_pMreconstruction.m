function [prdData, info] = predict_pMreconstruction(par, data, auxData)
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
% notice that the knot-abscissa for the scaled functional response
% trajectory is passed as a global here
% the knot-abscissa are in a 1-n vector and is defined in the run file
global tpM tf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% gw150 

%% First get f over time back
% get tyf back
yf = zeros(length(tf),1); % create empty vector with knot-coordinates for the f values
for i = 1:length(tf)
eval(['yf(i) = f_',num2str(tf(i)),';'])
end
tyf = [tf yf];

%% Then create p_M vectors
% prepare input matrix to the ODE called typM - this is a 2-n matrix with
% the first colum the knot-abscissa and the second column the knot
% coordinates. The knot coordinates are free parameters and are passed into
% the predict via the 'par' structure which is defined in the pars_init
% file

[data, auxData, metaData, txtData, weights] = mydata_pMreconstruction();

for tank = 1:length(fieldnames(data))
    tanknames=fieldnames(data);
    tankname=char(tanknames(tank));
    
    ypM = zeros(length(tpM),1); % create empty vector with knot-coordinates for the f values
    for i = 1:length(tpM)
    filterChecks = ( eval(['pM_', tankname, '_', num2str(tpM(i)),' < 0']));   
      if filterChecks 
          info = 0;
          prdData = {};
          return;
      end
    eval(['ypM(i) = pM_',tankname, '_', num2str(tpM(i)),';'])
    end
    eval (['typM_', tankname, ' = [tpM ypM]', ';']);
end

% for tank = 1:length(fieldnames(data))
%     tanknames=fieldnames(data);
%     tankname=char(tanknames(tank));
%     eval(['typM', '_', tankname, '=', 'typM', ';']);
% end


%% compute temperature correction factors
  TC = tempcorr(C2K(8.5), T_ref, T_A);


  % initial conditions -
pars_UE0          = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
[U_E0, Lb, info]  = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve, f= 1
E0            = U_E0 * p_Am; % J, initial energy in the egg
LEH0  = [1e-5; E0; 0]; % 3-1 vector [J, cm, J] with initial conditions



  % predict for all tanks
  
for tank = 1:length(fieldnames(data))
    tanknames=fieldnames(data);
    tankname=char(tanknames(tank));  

    time = eval(['data.', tankname, '(:,1);']);    % gets the times from the data of the specific tank
    
    eval(['typM = typM', '_', tankname, ';']);          % gets the typM for the specific tank
    
    [~, LEH] =  ode23s(@dget_LEH_for_pMreconstr, time, [LEH0; 0; 0; 0],[],tyf, typM,  TC, par, cPar); 
    EW = LEH(:,1).^3 + w_E/ mu_E * LEH(:,2)/ d_E; % g, wet weight

    eval(['prdData.', tankname, ' = EW;']) ;       % saves the predictions in prdData for specific tank            
end

