function [prdData, info] = predict_tW(par, data, auxData)                    

global string 

cPar  = parscomp_st(par); 
vars_pull(par); vars_pull(cPar); vars_pull(data); vars_pull(auxData);

TC = tempcorr(T, T_ref, T_A);

age = data.tW(:,1);
if age(1) > 0
    ageIn = [0; age];
else
    ageIn = age;
end
[~, LEH] =  ode23s(@dget_LEH, ageIn, [LEH0; 0; 0; 0],[],f, TC, par, cPar, treatment); 
if age(1) > 0
    LEH(1,:) = [];
end

% the three zeros are for N, Lb, and Lj (see inside the ODE subfunction)

L = LEH(:,1);   % cm, structural length
E = LEH(:,2);   % J, energy in reserve
EW = L.^3 + w_E/ mu_E * E/ d_E; % g, wet weight

prdData.tW = EW;


function dLEH = dget_LEH(t, LEH, f, TC, p, c, treatment)

  
  L = LEH(1);     % cm, structural length
  E = LEH(2);     % J, reserve E
  E_H = LEH(3);   % J, maturity E_H
  %E_R = LEH(4);    % #/d, cum reproductive output
  Lb = LEH(5);    % cm, length
  Lj = LEH(6);    % cm, length

 
 %  Shape correction function:
 % are numberical inaccuracies a problem ? (sta 16/04/17)
  if E_H < p.E_Hb
      s_M = 1;
  elseif E_H < p.E_Hj && E_H >= p.E_Hb
      s_M = L/ Lb;
  else
      s_M = Lj/ Lb;
  end
  
  %  Shape correction function applies to surface-area specific
  %  assimilation and energy conductance:
  % pA = 0 before birth E_H = E_Hb
  % all parameters with time in dimension need to be temp corrected
  p.v    = p.v * s_M * TC; % cm/d, conductance
  c.p_Am = c.p_Am * s_M * TC;
  p.p_M  = p.p_M * TC;
  p.k_J  = p.k_J * TC;
  
  switch treatment
      
      case 'p_M'
        p_M_Q = p.p_M * p.delta; % J/d/cm^3, p_M at start at T
            if p_M_Q < p.p_M
            p_M_t = min(p.p_M, p_M_Q + (p.p_M - p_M_Q)/ p.t_f * t);
            else
            p_M_t = max(p.p_M, p_M_Q + (p.p_M - p_M_Q)/ p.t_f * t);
            end
       p.p_M = p_M_t;
  
      case 'E_G'         
      E_G_Q = p.E_G * p.delta; % J/d/cm^3, p_M of exposed organism
      E_G_t = max(p.E_G, E_G_Q + (p.E_G - E_G_Q)/ p.t_f * t);
      p.E_G = E_G_t; % overwrite control value
      fprintf('E_G_t = %3.1f \n',E_G_t);
  end

  % Growth rate and mobilization rate:
  pA  = f * c.p_Am * L^2 * (E_H >= p.E_Hb); % J/cm^2/d, maximum surface area-specific assimilation rate
  L2  = L * L; L3 = L2 * L; L4 = L3 * L;
  pC  = E/L3 * (p.E_G * p.v/ L  + p.p_M)/ (p.kap * E/ L3 + p.E_G);   % [p_C], J/cm^3 (2.12, Kooijman 2010)
    if p.kap * pC < p.p_M  
    r = (E * p.v/ L4  - p.p_M/ p.kap)/ (E/ L3 + p.E_G * c.kap_G/ p.kap); % 1/d, specific growth rate
    else
    r = (E * p.v/ L4 - p.p_M/ p.kap)/ (E/ L3 + p.E_G/ p.kap); % 1/d, specific growth rate
    end
  
  % generate dH/da, dE/da, dL/da:  
  dE_H =  ((1 - p.kap) * pC * L3 - p.k_J * E_H)  * (E_H < p.E_Hp);
  dE_R =  p.kap_R * ((1 - p.kap) * pC * L3 - p.k_J * p.E_Hp) * (E_H >= p.E_Hp);
%   dN   =  p.kap_R * dE_R/ E0; % #/d, cum reproductive output
  dE   =  pA - pC * L3;
  dL   =  L * r/ 3;
  dLj  =  dL * (E_H <= p.E_Hj); % extract Lj dynamically
  dLb  =  dL * (E_H <= p.E_Hb); % to extract Lb dynamically with changing parameter
  
  % Pack output 
  dLEH = [dL; dE; dE_H; dE_R; dLb; dLj]; 