function dLEH_pM = dget_LEH2(t, LEH_pM, f, TC, p, c, dp_M)
% pMoA: physiological mode of action
  
  L = LEH_pM(1);     % cm, structural length
  E = LEH_pM(2);     % J, reserve E
  E_H = LEH_pM(3);   % J, maturity E_H
  E_R = LEH_pM(4);    % #/d, cum reproductive output
  Lb = LEH_pM(5);    % cm, length
  Lj = LEH_pM(6);    % cm, length
  p_M = LEH_pM(7);    % p_M

 
 %  Shape correction function:
 % are numerical inaccuracies a problem ? (sta 16/04/17)
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
  
  if p_M > p_Mcontrol
      dp_M = -delta ;
  else
      dp_M = 0;
  end
   
  
  % Pack output 
  dLEH_pM = [dL; dE; dE_H; dE_R; dLb; dLj; dp_M]; 