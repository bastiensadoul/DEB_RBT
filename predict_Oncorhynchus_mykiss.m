function [prdData, info] = predict_Oncorhynchus_mykiss(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
   
  if E_Hh >= E_Hb
   info = 0; prdData = {}; return
  end

%% compute temperature correction factors
  TC_ah_8_5 = tempcorr(temp.ah_8_5, T_ref, T_A);
  TC_ab_8_5 = tempcorr(temp.ab_8_5, T_ref, T_A);
  TC_Tah = tempcorr(Tah(:,1), T_ref, T_A);
  TC_Tab = tempcorr(Tab(:,1), T_ref, T_A);
  TC_Tah_Velsen = tempcorr(Tah_Velsen(:,1), T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  TC_tW = tempcorr(temp.tW, T_ref, T_A);
  TC_tT_Davidson2014 = tempcorr(mean(temp.tT_Davidson2014(:,2)), T_ref, T_A);
  TC_gw150meancontrol = tempcorr(temp.tW_gw150meancontrol, T_ref, T_A);
  TC_gw124bvarmeancontrol = tempcorr(temp.tW_gw124bvarmeancontrol, T_ref, T_A);


  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

  % initial 
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve

  % hatch
  [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  a_h = aUL(2,1);                 % d, age at hatch at f and 8.5C
  aT_h_8_5 = a_h/ TC_ah_8_5;      % d, age at hatch at f and 8.5C

  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth of foetus  at f = 1
  Lw_b = L_b/ del_M;                % cm, total length at birth
  aT_b_8_5 = t_b/ k_M/ TC_ab_8_5;   % d, age at birth at f and 8.5C

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, total length at puberty at f
  aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at f and T

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
  prdData.ah_8_5 = aT_h_8_5;
  prdData.ab_8_5 = aT_b_8_5;
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lb = Lw_b;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wi = Ww_i;
  prdData.Ri = RT_i;

  %% uni-variate data
  
  % t-Ww-data from YaniHisa2002
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tW);
  rT_B = TC_tW * rho_B * k_M; rT_j = TC_tW * rho_j * k_M; % 1/d, von Bert, exponential growth rate
  aT_b = t_b/ k_M/ TC_tW; aT_j = t_j/ k_M/ TC_tW;
  L_b = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;
  W_0 = initWeight.tW;
  L_0 = (W_0/ (1 + f_tW * w))^(1/3); % cm, structural length at t = 0 
  if L_0 < L_j
    tj = log(L_j/ L_0) * 3/ rT_j;
    t_bj = tW(tW(:,1) < tj,1); % select times between birth & metamorphosis
    L_bj = L_0 * exp(t_bj * rT_j/3); % exponential growth as V1-morph
    t_ji = tW(tW(:,1) >= tj,1); % selects times after metamorphosis
    L_ji = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - tj)); % cm, expected length at time
    L = [L_bj; L_ji]; % catenate lengths
  else 
    L = L_i - (L_i - L_0) * exp( - rT_B * tW(:,1)); % cm, expected length at time
  end
  EWw = L.^3 * (1 + f_tW * w);
  
  % T-ah and T_ab from From1991:
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  a_h = aUL(2,1);                 % d, age at hatch at f and T_ref
  a_b = aUL(3,1);                 % d, age at hatch at f and T_ref
  EaT_h =  a_h ./ TC_Tah;              % d, age at hatch at f and T
  EaT_b =  a_b ./ TC_Tab;              % d, age at birth at f and T
  
  % T-ah from Velsen:
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  a_h = aUL(2,1);                 % d, age at hatch at f and T_ref
  EaT_h_Velsen =  a_h ./ TC_Tah_Velsen;              % d, age at hatch at f and T
  
% Davidson2014
U_E0 = initial_scaled_reserve(f_tWL_Davidson2014, pars_UE0); % d.cm^2, initial scaled reserve
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
a_h   = aUL(2,1);                 % d, age at hatch at f and T_ref
U_Eh  = aUL(2,2);                 % d cm^2, scaled reserve at hatch at T_ref
UT_Eh = U_Eh/ TC_tT_Davidson2014; % d cm^2, scaled reserve at hatch at T
UT_Hh = E_Hh/ (p_Am * TC_tT_Davidson2014); % d cm^2, scaled maturity at hatch at T
L_h   = aUL(2,3);                 % cm, structural length at hatch
aT_h  = a_h/ TC_tT_Davidson2014;  % d, age at hatch at f and T  
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tWL_Davidson2014);
kT_M  = k_M * TC_tT_Davidson2014; rT_B = rho_B * kT_M; rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate
L_b   = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;
aT_b  = t_b/ kT_M; aT_j = t_j/ kT_M;


% tL_Davidson2014
tL = tL_Davidson2014(:,1);  
t_emb = tL(tL < aT_b - aT_h,1); % selects times after hatch and before birth    
if isempty(t_emb) == 0
pars = [v * TC; g; L_m; k_J * TC; kap];
 if t_emb(1) > 0
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
t_bj = tL(tL >= aT_b - aT_h & tL < aT_j - aT_h,1); % selects times after hatch and before birth
t_ji = tL(tL >= aT_j - aT_h,1); % selects times after metamorphosis
% time-length 
L_bj = L_b * exp(t_bj * rT_j/ 3);
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j - aT_h)); % cm, expected length at time
EL_Davidson2014 = [L_emb; L_bj; L_jm]/ del_M;                   % cm, structural length
  
% T-Ww davidson2014
tWw = tW_Davidson2014(:,1);
t_emb = tWw((tWw <= aT_b - aT_h),1); % selects times after hatch and before birth    
if isempty(t_emb) == 0
pars = [v * TC; g; L_m; k_J * TC; kap];
 if t_emb(1) > 0
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
E_emb = ULH(:,1) * p_Am * TC;   % J, embryo energy in reserve
Ww_emb = d_V * L_emb.^3 + w_E/ mu_E * E_emb; % g, embryo wet weight
else
Ww_emb = [];
end
t_bj = tWw(tWw >= aT_b - aT_h & tWw < aT_j - aT_h,1); % selects times after hatch and before birth
t_ji = tWw(tWw >= aT_j - aT_h,1); % selects times after metamorphosis
L_bj = L_b * exp(t_bj * rT_j/ 3); Ww_bj = L_bj.^3 * (1 + w * f);
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j - aT_h)); Ww_jm = L_jm.^3 * (1 + w * f);
EW_Davidson2014 = [Ww_emb; Ww_bj; Ww_jm];  % g, wet weight
   
% time -weight tW_gw150meancontrol
tWw = tW_gw150meancontrol(:,1);
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tW_gw150meancontrol);
kT_M = k_M * TC_gw124bvarmeancontrol; 
rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; L_b = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;
aT_b = t_b/ kT_M; aT_j = t_j/ kT_M;
tT_j = aT_j - aT_b; % d, time since birth at metamorphosis
L_bj = L_b * exp(tWw((tWw <= tT_j),1) * rT_j/ 3); L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWw((tWw > tT_j),1) - tT_j)); 
Ww_bj = L_bj.^3 * (1 + w * f);   Ww_jm = L_jm.^3 * (1 + w * f); 
EW_gw150meancontrol = [Ww_bj; Ww_jm]; % g, wet weight
  
% tW_gw124bvarmeancontrol-data
tWw = tW_gw124bvarmeancontrol(:,1);
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tW_gw124bvarmeancontrol);
kT_M = k_M * TC_gw124bvarmeancontrol; 
rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
L_b = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m; aT_b = t_b/ kT_M; aT_j = t_j/ kT_M;
tT_j = aT_j - aT_b; % d, time since birth at metamorphosis
L_bj = L_b * exp(tWw((tWw(:,1) <= tT_j),1) * rT_j/ 3); L_jm = L_i - (L_i - L_j) * exp( - rT_B * (tWw((tWw(:,1) > tT_j),1) - tT_j)); 
Ww_bj = L_bj.^3 * (1 + w * f);   Ww_jm = L_jm.^3 * (1 + w * f); 
EW_gw124bvarmeancontrol = [Ww_bj; Ww_jm]; % g, wet weight
 
% pack to output
prdData.tW = EWw;
prdData.Tah = EaT_h ;
prdData.Tab = EaT_b ;
prdData.Tah_Velsen = EaT_h_Velsen ;
prdData.tL_Davidson2014 = EL_Davidson2014;
prdData.tW_Davidson2014 = EW_Davidson2014 ;
prdData.tW_gw150meancontrol = EW_gw150meancontrol ;
prdData.tW_gw124bvarmeancontrol = EW_gw124bvarmeancontrol ;

%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
