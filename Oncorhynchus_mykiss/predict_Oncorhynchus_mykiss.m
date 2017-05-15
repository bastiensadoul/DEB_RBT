function [prdData, info] = predict_Oncorhynchus_mykiss(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  
  % customized filters to contrain a parameter - optional - delete if not
%   needed
  filterChecks = f_tW3>1 || f_tW3 <0 || ...         % f contrained to not be larger than 1
                 f_tW1>1 || f_tW1 <0 || ...
                 f_tW2>1 || f_tW2 <0  || ...
                 E_Hh > E_Hb || E_Hh <= 0; % 
  
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  end  
  
  
% temperature correction factors
  TC_ah5  = tempcorr(temp.ah_5, T_ref, T_A);
  TC_ah   = tempcorr(temp.ah, T_ref, T_A);  
  TC_ap   = tempcorr(temp.ap, T_ref, T_A);
  TC_am   = tempcorr(temp.am, T_ref, T_A);
  TC_Ri   = tempcorr(temp.Ri, T_ref, T_A);
  TC_tW  = tempcorr(temp.tW1, T_ref, T_A); % DaviKlemm2014
  TC_tW3   = tempcorr(temp.tW3, T_ref, T_A); % YaniHisa2002
  TC_Tah  = tempcorr(Tah(:,1), T_ref, T_A);
  TC_tWde = tempcorr(temp.tWde, T_ref, T_A);
  TC_WJO = tempcorr(temp.WJO, T_ref, T_A);
  TC_Wie1985 = tempcorr(temp.Wie1985, T_ref, T_A);
  
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
  aT_h5 = aUL(2,1)/ TC_ah5;              % d, age at hatch at f and T
  aT_h = aUL(2,1)/ TC_ah;              % d, age at hatch at f and T
  L_h = aUL(2,3); % cm, strucural length at hatch
  E_h = aUL(2,2) * p_Am; % J, energy in reserves at hatch
  Wdh = (d_V * L_h^3 + w_E/ mu_E * E_h); % g, dry weight at hatch
  
  % birth
  L_b   = L_m * l_b;                  % cm, structural length at birth of foetus  at f = 1
  aT_b5 = t_b/ k_M/ TC_ah5;           % d, age at birth of foetus at f and T
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
  prdData.ah_5 = aT_h5;
  prdData.Wdh = Wdh;
  prdData.ab = aT_b;  
  prdData.ab_5 = aT_b5;
  prdData.Wdb = Wdb;
  prdData.am = aT_m;
  prdData.Li = Lw_i;
  prdData.Wi = Ww_i;
  prdData.Ri = RT_i;

  %% uni-variate data with f = 1
  
  % T-ah - Vels1987 - at f and T
  EaT_h =  aUL(2,1) ./ TC_Tah;   % d, age at hatch at f and T
  
  % tWde and tWde_E - NinnStev2006 - at f and T
  vT = v * TC_tWde; kT_J = TC_tWde * k_J; kT_M = k_M * TC_tWde; pT_Am = p_Am * TC_tWde;
  UT_E0 = U_E0/ TC_tWde; aT_b = t_b/ kT_M; aT_j = t_j/ kT_M; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; 
  L_j = l_j * L_m;
  if tWde_E(1,1) > 0
     time = [0; tWde_E(:,1)];
  else
     time = tWde_E(:,1);
  end
  % embryo yolk
  [a, LUH] = ode45(@dget_LUH, time, [1e-10 UT_E0 0], [], kap, vT, kT_J, g, L_m);
  EWde_E = 1e3 * max(0, LUH(:,2) * pT_Am * w_E/ mu_E - f * m_Em * d_V * LUH(:,1) .^ 3 );
  % embryo body mass
  t_bj = time(time >= aT_b & time < aT_j); t_ji = time(time >=  aT_j); 
  L_0b = LUH(time < aT_b,1);
  L_bj = L_b * exp((t_bj - aT_b) * rT_j/ 3);
  L_ji = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j)); % cm, expected length at time
  Wd_0b = (1 + f * m_Em) * d_V * L_0b.^3;
  Wd_bi = [L_bj;L_ji].^3 * d_V * (1 + f * w); % predicted dry structure in mg
  EWde  = 1e3 * [Wd_0b; Wd_bi];  % mg, yolk-free embryo dry mass
  % remove the first values in case we had appended a zero in front
  if tWde_E(1,1) >0
  EWde_E =  EWde_E(2:end);
  EWde = EWde(2:end);     
  end
  
  % Oxygen consumption   
  
  % KieAls1998
  L = forkLength.WJO * del_M ; % structural length
  pACSJGRD = p_ref * scaled_power_j(L, f, pars_JO, l_b, l_j, l_p);  % J/d, powers
  J_M = - (n_M\n_O) * eta_O * pACSJGRD(:, [1 7 5])';  % mol/d: J_C, J_H, J_O, J_N in rows
  EJO = - J_M(3,:)' .* TC_WJO * 1e3;         % mmol O2/d, O2 consumption 

  % Wie1985   
  L = (Wie1985(:,1) /(1 + f * w) ) .^ (1/3) ;  % estimated structural length from weights
  pACSJGRD = p_ref * scaled_power_j(L, f, pars_JO, l_b, l_j, l_p);  % J/d, powers
  J_M = - (n_M\n_O) * eta_O * pACSJGRD(:, [1 7 5])';  % mol/d: J_C, J_H, J_O, J_N in rows
  EJO_Wie1985 = - J_M(3,:)' .* TC_Wie1985 * 1e3;      % mmol O2/d, O2 consumption 

  % L-Ww,  
  EL3 = (LW(:,1) * del_M2).^3 * (1 + f * w); % g, wet mass
  
  
    
  % t-L and t-Ww, a-p, L_p Ww_p DaviKenn2014
  TC = TC_tW;
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tW1);
  % puberty 
  L_p  = L_m * l_p;                 % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, total length at puberty at f
  aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at f and T
  Ww_p = L_p^3 * (1 + f_tW1 * w);   % g, wet weight at puberty at f

  [EWw, EL1]  = predict_WL(f_tW1, TC, tL1(:,1), par, cPar);
  EL1 = EL1/ del_M; % cm, structural length
  EW1 = predict_WL (f_tW1, TC, tW1(:,1), par, cPar); % g, wet weight
  
  [EWw, EL2]  = predict_WL(f_tW2, TC, tL2(:,1), par, cPar);
  EL2 = EL2/ del_M; % cm, structural length
  EW2 = predict_WL (f_tW2, TC, tW2(:,1), par, cPar); % g, wet weight

  % t-Ww , YaniHisa2002
  TC = TC_tW3; % -, temp correction factor
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tW3);
  rT_B = TC * rho_B * k_M; rT_j = TC * rho_j * k_M; % 1/d, von Bert, exponential growth rate
  L_j  = l_j * L_m; L_i = l_i * L_m;
  L_0  = (W0.tW3/ (1 + f_tW3 * w))^(1/3); % cm, structural length at t = 0 
  if L_0 < L_j
    tj = log(L_j/ L_0) * 3/ rT_j; % d, time since beginning of experiment that metamorphosis occurs
    t_bj = tW3(tW3(:,1) < tj,1); % select times between birth & metamorphosis
    L_bj = L_0 * exp(t_bj * rT_j/3); % exponential growth as V1-morph
    t_ji = tW3(tW3(:,1) >= tj,1); % selects times after metamorphosis
    L_ji = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - tj)); % cm, expected length at time
    L = [L_bj; L_ji]; % catenate lengths
  else 
    L = L_i - (L_i - L_0) * exp( - rT_B * tW3(:,1)); % cm, expected length at time
  end
  EW3 = L.^3 * (1 + f_tW3 * w);
  
  % pack to output
   % pack to output:
  prdData.Tah      = EaT_h;
  prdData.tWde_E   = EWde_E;
  prdData.tWde     = EWde; 
  prdData.WJO      = EJO;
  prdData.Wie1985  = EJO_Wie1985; 
  prdData.LW      = EL3;
  prdData.tL1       = EL1; % DaviKlemm2014
  prdData.tL2       = EL2; % DaviKlemm2014
  prdData.tW1       = EW1; % DaviKlemm2014
  prdData.tW2       = EW2; % DaviKlemm2014
  prdData.tW3       = EW3; % YaniHisa2002
   
  prdData.ap       = aT_p; % DaviKlemm2014
  prdData.Lp       = Lw_p; % DaviKlemm2014
  prdData.Wp       = Ww_p; % DaviKlemm2014

%% Subfunctions :
  
function [EW, EL] = predict_WL(f, TC, timeSinceHatch, p, c)                    
% [EW, EL] = predict_WL(f, TC, timeSinceHatch, p, c)
% Inputs:
% f, scalar, scaled func response
% TC, scalar, temperature correction
% timeSinceHatch, n-vector, time since hatch
% p, structure with parameters
% c, structure with compound parameters
  
pars_tj = [c.g c.k c.l_T c.v_Hb c.v_Hj c.v_Hp];   
pars_UE0 = [c.V_Hb; c.g; p.k_J; c.k_M; p.v]; % compose parameter vector
      
% Calculate Parameters
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
UT_E0 = U_E0/ TC; % cm * d , scaled initial reserve at T
[U_H, aUL] = ode45(@dget_aul, [0; c.U_Hh; c.U_Hb], [0 U_E0 1e-10], [], p.kap, p.v, p.k_J, c.g, c.L_m);
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);

a_h   = aUL(2,1);                 % d, age at hatch at f and T_ref
aT_h  = a_h/ TC;                % d, age at hatch at f and T  
kT_M  = c.k_M * TC; 
rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate between first feeding and end of V1-morph period
rT_B = rho_B * kT_M; % 1/d, von Bert, exponential growth rate after V1-morph period
L_b = l_b * c.L_m; L_j = l_j * c.L_m; L_i = l_i * c.L_m;     % cm, length at birth, metamorphosis, ultimate
aT_b  = t_b/ kT_M; aT_j = t_j/ kT_M;    % d, age at birth, metamorphosis at T

t = timeSinceHatch + aT_h; % d, age since fertilization
t_0b = t(t < aT_b,1);    % ages during the embryo period
t_bj = t(t >= aT_b & t < aT_j,1);    % selects times during V1-morph period
t_ji = t(t >= aT_j,1);                % selects times after metamorphosis

if isempty(t_0b) == 0     % if t_emb is not empty    
t_0b = [0;t_0b];   
[a, LUH] = ode45(@dget_LUH, t_0b, [1e-4 UT_E0 0], [], p.kap, p.v * TC, p.k_J * TC, c.g, c.L_m);
    if length(t_0b) == 2
    LUH = LUH(end,:);
    else
    LUH = LUH(2:end,:);    
    end
    L_emb = LUH(:,1);   % cm, embryo structural length
    E_emb = LUH(:,2) * c.p_Am * TC;   % J, embryo energy in reserve
    Ww_emb = p.d_V * L_emb.^3 + c.w_E/ p.mu_E * E_emb; % g, embryo wet weight
else
L_emb = []; Ww_emb = [];
end

% time-length 
L_bj = L_b * exp((t_bj - aT_b) * rT_j/ 3); % cm length and weight during V1-morph period
Ww_bj = L_bj.^3 * (1 + c.w * f);   % g weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j));   % cm, length after V1-morph period
Ww_jm = L_jm.^3 * (1 + c.w * f); % g, weight after V1-morph period

EL = [L_emb; L_bj; L_jm]; % cm, structural length
EW = [Ww_emb; Ww_bj; Ww_jm]; % g, wet weight


