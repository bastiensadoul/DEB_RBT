function [prdData, info] = predict_freconstruction(par, data, auxData)
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
% notice that the knot-abscissa for the scaled functional response
% trajectory is passed as a global here
% the knot-abscissa are in a 1-n vector and is defined in the run file
global tf 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% gw150 

% prepare input matrix to the ODE called tyf - this is a 2-n matrix with
% the first colum the knot-abscissa and the second column the knot
% coordinates. The knot coordinates are free parameters and are passed into
% the predict via the 'par' structure which is defined in the pars_init
% file
yf = zeros(length(tf),1); % create empty vector with knot-coordinates for the f values
for i = 1:length(tf)
filterChecks = ( eval(['f_',num2str(tf(i)),' < 0'])  || eval(['f_',num2str(tf(i)),' > 1']));   
  if filterChecks 
      info = 0;
      prdData = {};
      return;
  end
eval(['yf(i) = f_',num2str(tf(i)),';'])
end

tyf = [tf yf];


%% compute temperature correction factors
  TC = tempcorr(C2K(8.5), T_ref, T_A);


  % initial conditions -
pars_UE0          = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
[U_E0, Lb, info]  = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve, f= 1
E0            = U_E0 * p_Am; % J, initial energy in the egg
LEH0  = [1e-5; E0; 0]; % 3-1 vector [J, cm, J] with initial conditions

  % predict 
time = data.tW_gw150A(:,1);
[~, LEH] =  ode23s(@dget_LEH_for_freconstr, time, [LEH0; 0; 0; 0],[],tyf, TC, par, cPar);
EW = LEH(:,1).^3 + w_E/ mu_E * LEH(:,2)/ d_E; % g, wet weight


prdData.tW_gw150A = EW ;
prdData.tW_gw150B = EW ;
prdData.tW_gw150C = EW ;
