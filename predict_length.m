function [EL_Data] = predict_length(f, TC_t, tL_data, dp='dp')     %tL_data (1st column is time in dpf, second is length)    
                                                              %the last parameter is the type of age provided in tL_data. : 
                                                                      % -
                                                                      % 'dpf' for days post fertilization    
                                                                      % 'dph' for days post hatch                  
                                                                      % 'dp'  for unknown (default)              
                                                              
% Parameters
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
a_h   = aUL(2,1);                 % d, age at hatch at f and T_ref
U_Eh  = aUL(2,2);                 % d cm^2, scaled reserve at hatch at T_ref
UT_Eh = U_Eh/ TC_t;               % d cm^2, scaled reserve at hatch at T
UT_Hh = E_Hh/ (p_Am * TC_t);      % d cm^2, scaled maturity at hatch at T
L_h   = aUL(2,3);                 % cm, structural length at hatch
aT_h  = a_h/ TC_t;                % d, age at hatch at f and T  
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
kT_M  = k_M * TC_t; 
rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate between first feeding ang end of V1-morph period
rT_B = rho_B * kT_M; % 1/d, von Bert, exponential growth rate after V1-morph period
L_b   = l_b * L_m;   % cm, length at first feeding
L_j = l_j * L_m;     % cm, length at metamorphosis
L_i = l_i * L_m;     % cm, ultimate length
aT_b  = t_b/ kT_M;   % d, age at birth
aT_j = t_j/ kT_M;    % d, age at end V1-morph period

% tL_data
tL = tL_data(:,1);
if dp=='dpf'
  t_emb = tL(tL < aT_b & tL > aT_h,1);  
  %t_emb = t_emb - aT_h;
end
if dp=='dph'
  tL=tL+aT_h
  t_emb = tL(tL < aT_b & tL > aT_h,1);  
end
if dp=='dp'
  t_emb = [];  
end
 
if isempty(t_emb) == 0          % if t_emb is not empty
pars = [v * TC_t; g; L_m; k_J * TC_t; kap];
 if t_emb(1) > aT_h                % if first value is after hatch
 t_emb = [0;t_emb];    
 [AA, ULH] = ode23s(@dget_ulh_modified, t_emb, [UT_Eh; L_h; UT_Hh],[],pars); 
 ULH(1,:) = [];
 else 
 [AA, ULH] = ode23s(@dget_ulh_modified, t_emb, [UT_Eh; L_h; UT_Hh],[],pars);
 end 
 if length(t_emb) == 2
 ULH = ULH([1 end],:);
 end
 L_emb = ULH(:,2);   % cm, embryo structural length
else
L_emb = [];
end


t_bj = tL(tL >= aT_b & tL < aT_j,1);        % selects times during V1-morph period
t_ji = tL(tL >= aT_j,1);                    % selects times after metamorphosis
% time-length 
L_bj = L_b * exp(t_bj * rT_j/ 3);                                % cm, expected length during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j - aT_h));  % cm, expected length after V1-morph period

EL_Data = [L_emb; L_bj; L_jm]/ del_M;                            % cm, structural length



function dULH = dget_ulh_modified(t,ULH,p)
  % change in state variables during embryo stage
  % called by fnget_lnpars_r, get_pars

  % parameter vector : [v ;g; Lm ;kJ; kap]
  v = p(1); g = p(2); Lm = p(3); kJ = p(4); kap = p(5);
  
  % unpack state variables
  U = ULH(1); % U = M_E/{J_{EAm}}
  L = ULH(2); % structural length
  H = ULH(3); % H = M_H/{J_{EAm}}

  eL3 = U * v; % eL3 = L^3 * m_E/ m_Em
  gL3 = g * L^3;
  SC = L^2 * (1 + L/(g * Lm)) * g * eL3/ (gL3 + eL3); % J_EC/{J_EAm}
  dU = - SC;
  dL = v * (eL3 - L^4/ Lm)/ (3 * (eL3 + gL3));
  dH = (1 - kap) * SC - kJ * H;

  % pack derivatives
  dULH = [dU; dL; dH];


