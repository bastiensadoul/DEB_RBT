function [prdData, info] = predict_Oncorhynchus_mykiss_SSAF(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  
  % customized filters to contrain a parameter - optional - delete if not
%   needed
  filterChecks = f_tW>1 || f_tW <0 || ...         % f contrained to not be larger than 1
                 f_tWL>1 || f_tWL <0;   % 
  
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  end  
  
  
%% compute temperature correction factors
  TC_ah   = tempcorr(temp.ah, T_ref, T_A);  
  TC_ap   = tempcorr(temp.ap, T_ref, T_A);
  TC_am   = tempcorr(temp.am, T_ref, T_A);
  TC_Ri   = tempcorr(temp.Ri, T_ref, T_A);
  TC_tW   = tempcorr(temp.tW, T_ref, T_A);
   
% parameter vector for DEBtool:
  pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];
  pars_JO = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose pars
  p_ref = p_Am * L_m^2; % J/d, max assimilation power at max size

%% zero -variate data  

  % life cycle
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

  % initial 
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve  
  E_0 = U_E0 * p_Am;     % J, initial reserve
  Wd0 = E_0 * w_E/ mu_E ; % g, egg dry weight 
    
  % hatch
  [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  aT_h = aUL(2,1)/ TC_ah;              % d, age at hatch at f and T
  L_h = aUL(2,3); % cm, strucural length at hatch
  E_h = aUL(2,2) * p_Am; % J, energy in reserves at hatch
  Wdh = (d_V * L_h^3 + w_E/ mu_E * E_h); % g, dry weight at hatch
  
  % birth
  L_b   = L_m * l_b;                  % cm, structural length at birth of foetus  at f = 1
  aT_b  = t_b/ k_M/ TC_ah;            % d, age at birth of foetus at f and T
  Wdb   = d_V * L_b^3 * (1 + f * w);  % g, dry weight at birth at f 
     
  % puberty: this is moved to after uni-variate data, since is works with f_tLW

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector
  RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % ultimate reproduction rate

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T

  % pack to output
  prdData.Wd0 = Wd0;
  prdData.ah = aT_h;
  prdData.Wdh = Wdh;
  prdData.ab = aT_b;  
  prdData.Wdb = Wdb;
  prdData.am = aT_m;
  prdData.Li = Lw_i;
  prdData.Wi = Ww_i;
  prdData.Ri = RT_i;

  %% uni-variate data with f = 1
  
 
   
  % McKenPed2007
  
  % length - weight - respiration at f and T
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f); 
  % SSAF - small size at age family
  % Ww-JO
  L = (WwJO_1(:,1)/ (1 + f * w)) .^ (1/3);  % cm, structural length
  pACSJGRD = p_ref * scaled_power_j(L, f, pars_JO, l_b, l_j, l_p);  % J/d, powers
  J_M = - (n_M\n_O) * eta_O * pACSJGRD(:, [1 7 5])';  % mol/d: J_C, J_H, J_O, J_N in rows
  EJO = - J_M(3,:)' * TC_WJO * 1e3;         % mmol O2/d, O2 consumption 
  % t-Ww , and t-L
  rT_B = TC_WwJO * rho_B * k_M; rT_j = TC_WwJO * rho_j * k_M; % 1/d, von Bert, exponential growth rate
  L_j = l_j * L_m; L_i = l_i * L_m;
  L_0 = (76.5/ (1 + f * w))^(1/3); % cm, structural length at t = 0 
  if L_0 < L_j 
    tj = log(L_j/ L_0) * 3/ rT_j; % d, time since beginning of experiment that metamorphosis occurs
    t_bj = tWw(tWw(:,1) < tj,1); % select times between birth & metamorphosis
    L_bj = L_0 * exp(t_bj * rT_j/3); % exponential growth as V1-morph
    t_ji = tWw(tWw(:,1) >= tj,1); % selects times after metamorphosis
    L_ji = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - tj)); % cm, expected length at time
    L = [L_bj; L_ji]; % catenate lengths
  else 
    L = L_i - (L_i - L_0) * exp( - rT_B * tWw(:,1)); % cm, expected length at time
  end
  EW = L.^3 * (1 + f * w);
  EL =   spline(tL(:,1), [tWw(:,1) L] )./ del_M;
  
  
  % pack to output
  prdData.tW       = EW;
  prdData.tL      = EL;
  prdData.WJO   = EJO; 
 

