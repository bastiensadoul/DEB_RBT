function [prdData, info] = predict_reconstruction(par, data, auxData)
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
[~, LEH] =  ode23s(@dget_LEH_for_reconstr, time, [LEH0; 0; 0; 0],[],tyf, TC, par, cPar); 
EW = LEH(:,1).^3 + w_E/ mu_E * LEH(:,2)/ d_E; % g, wet weight

%   % predict 
% time = data.timescale(:,1);
% [~, LEH] =  ode23s(@dget_LEH, time, [LEH0; 0; 0; 0],[],tyf, TC, par, cPar); 

prdData.tW_gw150A = EW ;
prdData.tW_gw150B = EW ;
prdData.tW_gw150C = EW ;
%prdData.timescale = LEH(:, 1) ;
% 
% function dLEH = dget_LEH(t, LEH, tyf, TC, p, c)
% % pMoA: physiological mode of action
%   
%   L = LEH(1);     % cm, structural length
%   E = LEH(2);     % J, reserve E
%   E_H = LEH(3);   % J, maturity E_H
%   %E_R = LEH(4);    % #/d, cum reproductive output
%   Lb = LEH(5);    % cm, length
%   Lj = LEH(6);    % cm, length
%   
%   %f=0,6;
%   if t<64 
%       f=1;
%   else
%       f = interp1(tyf(:,1),tyf(:,2),t,'linear');
%       end
%  
%  %  Shape correction function:
%  % are numerical inaccuracies a problem ? (sta 16/04/17)
%   if E_H < p.E_Hb
%       s_M = 1;
%   elseif E_H < p.E_Hj && E_H >= p.E_Hb
%       s_M = L/ Lb;
%   else
%       s_M = Lj/ Lb;
%   end
%   
%   %  Shape correction function applies to surface-area specific
%   %  assimilation and energy conductance:
%   % pA = 0 before birth E_H = E_Hb
%   % all parameters with time in dimension need to be temp corrected
%   p.v    = p.v * s_M * TC; % cm/d, conductance
%   c.p_Am = c.p_Am * s_M * TC;
%   p.p_M  = p.p_M * TC;
%   p.k_J  = p.k_J * TC;
%   
% 
%   % Growth rate and mobilization rate:
%   pA  = f * c.p_Am * L^2 * (E_H >= p.E_Hb); % J/cm^2/d, maximum surface area-specific assimilation rate
%   L2  = L * L; L3 = L2 * L; L4 = L3 * L;
%   pC  = E/L3 * (p.E_G * p.v/ L  + p.p_M)/ (p.kap * E/ L3 + p.E_G);   % [p_C], J/cm^3 (2.12, Kooijman 2010)
%     if p.kap * pC < p.p_M  
%     r = (E * p.v/ L4  - p.p_M/ p.kap)/ (E/ L3 + p.E_G * c.kap_G/ p.kap); % 1/d, specific growth rate
%     else
%     r = (E * p.v/ L4 - p.p_M/ p.kap)/ (E/ L3 + p.E_G/ p.kap); % 1/d, specific growth rate
%     end
%   
%   % generate dH/da, dE/da, dL/da:  
%   dE_H =  ((1 - p.kap) * pC * L3 - p.k_J * E_H)  * (E_H < p.E_Hp);
%   dE_R =  p.kap_R * ((1 - p.kap) * pC * L3 - p.k_J * p.E_Hp) * (E_H >= p.E_Hp);
% %   dN   =  p.kap_R * dE_R/ E0; % #/d, cum reproductive output
%   dE   =  pA - pC * L3;
%   dL   =  L * r/ 3;
%   dLj  =  max(0,dL) * (E_H <= p.E_Hj); % extract Lj dynamically. max (0,dL) to be sure that Lj doesn't decrease even if dL is negative 
%   dLb  =  max(0,dL) * (E_H <= p.E_Hb); % to extract Lb dynamically with changing parameter
%   
%    
%   % Print in txt
% %   fileID = fopen(['test,'.txt'],'a');
% %   fprintf(fileID, [num2str(t),'\t',num2str(L),'\t',num2str(r),'\t',num2str(E),'\t', num2str(p.v), '\t', num2str(dL), '\t', num2str(p.E_G),'\t', num2str(p.kap),'\t', num2str(p.p_M),'\t',num2str(pC),'\t',num2str(c.kap_G),'\n']);
% %   fclose(fileID);
% %   
%   % Pack output 
%   dLEH = [dL; dE; dE_H; dE_R; dLb; dLj]; 
