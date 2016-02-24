function [prdData, info] = predict_Oncorhynchus_mykiss_BPA3and30(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
   
  
%    customized filters for allowable parameters of the standard DEB model (std)
%   for other models consult the appropriate filter function.
  filterChecks = k * v_Hp >= f_tW_gw150A_BPA3^3 ||  ~reach_birth(g, k, v_Hb, f_tW_gw150A_BPA3) || ...
             k * v_Hp >= f_tW_gw150B_BPA3^3 ||  ~reach_birth(g, k, v_Hb, f_tW_gw150B_BPA3) || ...
            k * v_Hp >= f_tW_gw150C_BPA3^3 ||  ~reach_birth(g, k, v_Hb, f_tW_gw150C_BPA3) || ...
            k * v_Hp >= f_tW_gw150A_BPA30^3 ||  ~reach_birth(g, k, v_Hb, f_tW_gw150A_BPA30) || ...
            k * v_Hp >= f_tW_gw150B_BPA30^3 ||  ~reach_birth(g, k, v_Hb, f_tW_gw150B_BPA30) || ...
            k * v_Hp >= f_tW_gw150C_BPA30^3 ||  ~reach_birth(g, k, v_Hb, f_tW_gw150C_BPA30);    
          
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  end  
%   

%% compute temperature correction factors
  TC_gw150 = tempcorr(temp.tW_gw150A_BPA3, T_ref, T_A);

  
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(1, pars_UE0); % d.cm^2, initial scaled reserve  
  E_0 = U_E0 * p_Am;     % J, initial reserve
  ELH_init = [E_0; 1e-4; 0];
  

% Study gw150_BPA3
p = [g, k, l_T, v_Hb, v_Hj, v_Hp];

[lj, lp, lb, info] = get_lj(p, f_tW_gw150A_BPA3); L_b = lb * L_m; L_j = lj * L_m;
[a, ELH] = ode45(@dget_ELH_pj, [0; tW_gw150A_BPA3(:,1)], ELH_init, [], L_b, L_j, L_m, p_Am, v, g, k_J, kap, f_tW_gw150A_BPA3, E_Hb, E_Hj, TC_gw150);
ELH(1,:) = []; E = ELH(:,1); L = ELH(:,2); EW_tW_gw150A_BPA3 = L.^3 + E * w_E/ mu_E/ d_E;  % g, wet weight

[lj, lp, lb, info] = get_lj(p, f_tW_gw150B_BPA3); L_b = lb * L_m; L_j = lj * L_m;
[a, ELH] = ode45(@dget_ELH_pj, [0; tW_gw150B_BPA3(:,1)], ELH_init, [], L_b, L_j, L_m, p_Am, v, g, k_J, kap, f_tW_gw150B_BPA3, E_Hb, E_Hj, TC_gw150);
ELH(1,:) = []; E = ELH(:,1); L = ELH(:,2); EW_tW_gw150B_BPA3 = L.^3 + E * w_E/ mu_E/ d_E;  % g, wet weight

[lj, lp, lb, info] = get_lj(p, f_tW_gw150C_BPA3); L_b = lb * L_m; L_j = lj * L_m;
[a, ELH] = ode45(@dget_ELH_pj, [0; tW_gw150C_BPA3(:,1)], ELH_init, [], L_b, L_j, L_m, p_Am, v, g, k_J, kap, f_tW_gw150C_BPA3, E_Hb, E_Hj, TC_gw150);
ELH(1,:) = []; E = ELH(:,1); L = ELH(:,2); EW_tW_gw150C_BPA3 = L.^3 + E * w_E/ mu_E/ d_E;  % g, wet weight

%-------------------------------------
% Study gw150_BPA30
p = [g, k, l_T, v_Hb, v_Hj, v_Hp];

[lj, lp, lb, info] = get_lj(p, f_tW_gw150A_BPA30); L_b = lb * L_m; L_j = lj * L_m;
[a, ELH] = ode45(@dget_ELH_pj, [0; tW_gw150A_BPA30(:,1)], ELH_init, [], L_b, L_j, L_m, p_Am, v, g, k_J, kap, f_tW_gw150A_BPA30, E_Hb, E_Hj, TC_gw150);
ELH(1,:) = []; E = ELH(:,1); L = ELH(:,2); EW_tW_gw150A_BPA30 = L.^3 + E * w_E/ mu_E/ d_E;  % g, wet weight

[lj, lp, lb, info] = get_lj(p, f_tW_gw150B_BPA30); L_b = lb * L_m; L_j = lj * L_m;
[a, ELH] = ode45(@dget_ELH_pj, [0; tW_gw150B_BPA30(:,1)], ELH_init, [], L_b, L_j, L_m, p_Am, v, g, k_J, kap, f_tW_gw150B_BPA30, E_Hb, E_Hj, TC_gw150);
ELH(1,:) = []; E = ELH(:,1); L = ELH(:,2); EW_tW_gw150B_BPA30 = L.^3 + E * w_E/ mu_E/ d_E;  % g, wet weight

[lj, lp, lb, info] = get_lj(p, f_tW_gw150C_BPA30); L_b = lb * L_m; L_j = lj * L_m;
[a, ELH] = ode45(@dget_ELH_pj, [0; tW_gw150C_BPA30(:,1)], ELH_init, [], L_b, L_j, L_m, p_Am, v, g, k_J, kap, f_tW_gw150C_BPA30, E_Hb, E_Hj, TC_gw150);
ELH(1,:) = []; E = ELH(:,1); L = ELH(:,2); EW_tW_gw150C_BPA30 = L.^3 + E * w_E/ mu_E/ d_E;  % g, wet weight


%--------------------------------------
% pack to output
prdData.tW_gw150A_BPA3 = EW_tW_gw150A_BPA3 ;
prdData.tW_gw150B_BPA3 = EW_tW_gw150B_BPA3 ;
prdData.tW_gw150C_BPA3 = EW_tW_gw150C_BPA3 ;
prdData.tW_gw150A_BPA30 = EW_tW_gw150A_BPA30 ;
prdData.tW_gw150B_BPA30 = EW_tW_gw150B_BPA30 ;
prdData.tW_gw150C_BPA30 = EW_tW_gw150C_BPA30 ;


% info = 1;

%% sub subfuctions

function dELH = dget_ELH_pj(t, ELH, Lb, Lj, Lm, p_Am, v, g, kJ, kap, f, Hb, Hj, TC)
%     dELH = dget_ELH_pj(t, ELH, Lb, Lj, Lm, p_Am, v, g, kJ, kap, f, Hb, Hj, tT, T_A, T_ref)
  %  change in state variables during juvenile stage
  %  dELH = dget_ELH_p_pj(t, ELH)
  %  ELH: 3-vector
  %   E: reserve E
  %   L: structural length
  %   H: maturity E_H
  %  dELH: change in reserve, length, scaled maturity
  
 
  %  unpack variables
  E = ELH(1); L = ELH(2); H = ELH(3);
  
%   TC = tempcorr(C2K(spline1(t, tT)), T_ref, T_A);
  vT = v * TC; kT_J = kJ * TC; pT_Am = p_Am * TC;
 
  if H < Hb 
    s = 1; % -, multiplication factor for v and {p_Am}
  elseif H < Hj
    s = L/ Lb;
  else
    s = Lj/ Lb;
  end
  e = vT * E/ L^3/ pT_Am; % -, scaled reserve density; 
  rT = s * vT * (e/ L - 1/ Lm/ s)/ (e + g); % 1/d, spec growth rate
  pT_C = E * (s * vT/ L - rT); % cm^2, scaled mobilisation
  
  % generate dH/dt, dE/dt, dL/dt
  dH = (1 - kap) * pT_C - kT_J * H;
  dE = (L > Lb) * s * pT_Am * f * L^2 * (H>= Hb) - pT_C;
  dL = rT * L/3;

  % pack derivatives
  dELH = [dE; dL; dH];
















% function [EW, EL] = predict_WL(f, TC, timeSinceFertilization, p, c)                    
% % [EW, EL] = predict_WL(f, TC, timeSinceHatch, p, c)
% % Inputs:
% % f, scalar, scaled func response
% % TC, scalar, temperature correction
% % timeSinceHatch, n-vector, time since hatch
% % p, structure with parameters
% % c, structure with compound parameters
% 
% % Life cycle parameters
% pars_tj = [c.g c.k c.l_T c.v_Hb c.v_Hj c.v_Hp];   
% [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
% 
% % Initial parameters
% pars_UE0 = [c.V_Hb; c.g; p.k_J; c.k_M; p.v]; % compose parameter vector
% U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
% UT_E0 = U_E0/ TC; % cm * d , scaled initial reserve at T
% 
% % Hatch
% [U_H, aUL] = ode45(@dget_aul, [0; c.U_Hh; c.U_Hb], [0 U_E0 1e-10], [], p.kap, p.v, p.k_J, c.g, c.L_m);
% a_h   = aUL(2,1);                 % d, age at hatch at f and T_ref
% aT_h  = a_h/ TC;                % d, age at hatch at f and T  
% 
% % Somatic maintenance coefficient corrected by T
% kT_M  = c.k_M * TC; 
% 
% % Birth
% aT_b  = t_b/ kT_M;   % d, age at birth
% L_b = l_b * c.L_m;   % cm, length at birth
% 
% % Metamorphosis
% aT_j = t_j/ kT_M;    % d, age at metamorphosis
% L_j = l_j * c.L_m;   % cm, length at metamorphosis
% 
% % Ultimate
%                      % no age, age = infinite
% L_i = l_i * c.L_m;   % cm, length at ultimate
% 
% % Von bert coefficients for difference periods
% rT_j = rho_j * kT_M; % 1/d, between first feeding and end of V1-morph period
% rT_B = rho_B * kT_M; % 1/d, after V1-morph period
% 
% %   % reproduction
% %   pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector
% %   RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % ultimate reproduction rate
% % 
% %   % life span
% %   pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
% %   t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
% %   aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
% 
% % Predictions during time
% 
% t = timeSinceFertilization;            % d, age since fertilization
% t_0b = t(t < aT_b,1);                 % ages during the embryo period
% t_bj = t(t >= aT_b & t < aT_j,1);     % selects times during V1-morph period
% t_ji = t(t >= aT_j,1);                % selects times after metamorphosis
% 
% if isempty(t_0b) == 0     % if t_emb is not empty    
% t_0b = [0;t_0b];   
% [a, LUH] = ode45(@dget_LUH, t_0b, [1e-4 UT_E0 0], [], p.kap, p.v * TC, p.k_J * TC, c.g, c.L_m);
%     if length(t_0b) == 2
%     LUH = LUH(end,:);
%     else
%     LUH = LUH(2:end,:);    
%     end
%     L_emb = LUH(:,1);                                  % cm, embryo structural length
%     E_emb = LUH(:,2) * c.p_Am * TC;                    % J, embryo energy in reserve
%     Ww_emb = p.d_V * L_emb.^3 + c.w_E/ p.mu_E * E_emb; % g, embryo wet weight
% else
% L_emb = []; Ww_emb = [];
% end
% 
% % time-length 
% L_bj = L_b * exp((t_bj - aT_b) * rT_j/ 3); Ww_bj = L_bj.^3 * (1 + c.w * f);   % cm,g, length and weight during V1-morph period
% L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j));                      % cm, length after V1-morph period
% Ww_jm = L_jm.^3 * (1 + c.w * f);                                              % g, weight after V1-morph period
% 
% EL = [L_emb; L_bj; L_jm]; % cm, structural length
% EW = [Ww_emb; Ww_bj; Ww_jm]; % g, wet weight
% 
% 
% 
