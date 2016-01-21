function [prdData, info] = predict_Oncorhynchus_mykiss(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
   
  if E_Hh >= E_Hb
   info = 0; prdData = {}; return
  end
  
   % customized filters for allowable parameters of the standard DEB model (std)
  % for other models consult the appropriate filter function.
  filterChecks = k * v_Hp >= f_tW^3 || ...         % constraint required for reaching puberty with f_tL
                 ~reach_birth(g, k, v_Hb, f_tW);   % constraint required for reaching birth with f_tL
  
  if filterChecks  
    info = 0;
    prdData = {};
    return;
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
  TC_Davidson2014 = tempcorr(mean(temp.tT_Davidson2014(:,2)), T_ref, T_A);
  TC_gw150 = tempcorr(temp.tW_gw150A, T_ref, T_A);
  TC_gw124 = tempcorr(temp.tW_gw124iniA, T_ref, T_A);
  TC_ind = tempcorr(temp.tL_ind, T_ref, T_A);
 
  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

  % initial 
  pars_UE0 = [V_Hb g k_J k_M v]; % compose parameter vector
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

  
%-------------------------------------
% Build in pred_pars, all the parameters needed in predict_length and predict_weight 
%(parameters common to all studies, independent from f and T

pred_pars=[pars_UE0 U_Hh U_Hb kap v k_J g L_m pars_tj k_M E_Hh p_Am w];
 
  
  %% uni-variate data
  
  
%-------------------------------------
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
  
%-------------------------------------
% T-ah and T_ab from From1991:
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
a_h = aUL(2,1);                 % d, age at hatch at f and T_ref
a_b = aUL(3,1);                 % d, age at hatch at f and T_ref
EaT_h =  a_h ./ TC_Tah;              % d, age at hatch at f and T
EaT_b =  a_b ./ TC_Tab;              % d, age at birth at f and T
  
%-------------------------------------
% T-ah from Velsen:
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
a_h = aUL(2,1);                 % d, age at hatch at f and T_ref
EaT_h_Velsen =  a_h ./ TC_Tah_Velsen;              % d, age at hatch at f and T

%-------------------------------------
% Davidson2014

% tL_Davidson2014
EL_Davidson2014 = predict_length(f_tWL_Davidson2014, TC_Davidson2014, tL_Davidson2014, 'dph', pred_pars);
EL_Davidson2014 = EL_Davidson2014/del_M;

% T-Ww davidson2014
EW_Davidson2014=predict_weight(f_tWL_Davidson2014, TC_Davidson2014, tW_Davidson2014, 'dph', pred_pars);

%-------------------------------------
% Study gw150

EW_tW_gw150A=predict_weight(f_tW_gw150, TC_gw150, tW_gw150A, 'dpf', pred_pars);
EW_tW_gw150B=predict_weight(f_tW_gw150, TC_gw150, tW_gw150B, 'dpf', pred_pars);
EW_tW_gw150C=predict_weight(f_tW_gw150, TC_gw150, tW_gw150C, 'dpf', pred_pars);

%-------------------------------------
% Study gw124ini

% let's try with adding 64 to the date like in 150 and say it is dpf
tW_gw124iniA(:,1)=tW_gw124iniA(:,1)+64;
tW_gw124iniB(:,1)=tW_gw124iniB(:,1)+64;
tW_gw124iniC(:,1)=tW_gw124iniC(:,1)+64;

EW_tW_gw124iniA=predict_weight(f_tW_gw124, TC_gw124, tW_gw124iniA, 'dpf', pred_pars);
EW_tW_gw124iniB=predict_weight(f_tW_gw124, TC_gw124, tW_gw124iniB, 'dpf', pred_pars);
EW_tW_gw124iniC=predict_weight(f_tW_gw124, TC_gw124, tW_gw124iniC, 'dpf', pred_pars);

%-------------------------------------
% Study gw124fin

% Same here
tW_gw124fin(:,1)=tW_gw124fin(:,1)+64;

EW_tW_gw124fin=predict_weight(f_tW_gw124, TC_gw124, tW_gw124fin, 'dpf', pred_pars);


%-------------------------------------
% Study ind

% tL_ind-data
EL_ind=predict_length(f_tLW_ind, TC_ind, tL_ind, 'dpf', pred_pars);
EL_tL_ind=EL_ind/del_M2;

% tW_ind-data
EW_tW_ind=predict_weight(f_tLW_ind, TC_ind, tL_ind, 'dpf', pred_pars);




%--------------------------------------
% pack to output
prdData.tW = EWw;
prdData.Tah = EaT_h ;
prdData.Tab = EaT_b ;
prdData.Tah_Velsen = EaT_h_Velsen ;
prdData.tL_Davidson2014 = EL_Davidson2014;
prdData.tW_Davidson2014 = EW_Davidson2014 ;
prdData.tW_gw150A = EW_tW_gw150A ;
prdData.tW_gw150B = EW_tW_gw150B ;
prdData.tW_gw150C = EW_tW_gw150C ;
prdData.tW_gw124iniA = EW_tW_gw124iniA ;
prdData.tW_gw124iniB = EW_tW_gw124iniB ;
prdData.tW_gw124iniC = EW_tW_gw124iniC ;
prdData.tW_gw124fin = EW_tW_gw124fin ;
prdData.tW_ind = EW_tW_ind;
prdData.tL_ind=  EL_tL_ind;

%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  
  
  
  %%%%%%%%%%------------------- get length   ------------------------------------------------------
  
  function [EL_Data] = predict_length(f, TC_t, tL_data, dp, pred_pars)                      %tL_data (1st column is time in dpf, second is length)    
                                                                               % 'dpf' for days post fertilization    
                                                                               % 'dph' for days post hatch 
% Extract common parameters
pars_UE0 = pred_pars(1:5);
U_Hh=pred_pars(6);
U_Hb=pred_pars(7);
kap=pred_pars(8);
v=pred_pars(9);
k_J=pred_pars(10);
g=pred_pars(11);
L_m=pred_pars(12);
pars_tj=pred_pars(13:18);
k_M=pred_pars(19);
E_Hh=pred_pars(20);
p_Am=pred_pars(21);
w=pred_pars(22);

                                                                        
% Calculate Parameters
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
rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate between first feeding and end of V1-morph period
rT_B = rho_B * kT_M; % 1/d, von Bert, exponential growth rate after V1-morph period
L_b = l_b * L_m;   % cm, length at birth
L_j = l_j * L_m;     % cm, length at metamorphosis
L_i = l_i * L_m;     % cm, ultimate length
aT_b  = t_b/ kT_M;   % d, age at birth
aT_j = t_j/ kT_M;    % d, age at end V1-morph period

% tL_data
tL = tL_data(:,1);

if strcmp(dp, 'dph') == 1                            % if data expressed in dph, need to add the predicted value for age at hatch
  tL=tL+aT_h;
  t_emb = tL(aT_h < tL & tL < aT_b,1);
end
t_emb = tL(aT_h < tL & tL < aT_b,1);    % available ages between hatch and birth in the dataset

if isempty(t_emb) == 0                  % if t_emb is not empty
pars = [v * TC_t; g; L_m; k_J * TC_t; kap];
 if t_emb(1) > aT_h                     % if first age value is after hatch
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

EL_Data = [L_emb; L_bj; L_jm];                            % cm, structural length




 %%%%%%%%%%------------------- get weight   ------------------------------------------------------


function [EW_Data] = predict_weight(f, TC_t, tW_data, dp, pred_pars)        %tL_data (1st column is time in dpf, second is length)    
                                                                               % 'dpf' for days post fertilization    
                                                                               % 'dph' for days post hatch
% Extract common parameters
pars_UE0 = pred_pars(1:5);
U_Hh=pred_pars(6);
U_Hb=pred_pars(7);
kap=pred_pars(8);
v=pred_pars(9);
k_J=pred_pars(10);
g=pred_pars(11);
L_m=pred_pars(12);
pars_tj=pred_pars(13:18);
k_M=pred_pars(19);  
E_Hh=pred_pars(20);
p_Am=pred_pars(21);
w=pred_pars(22);
                                                              
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
rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate between first feeding and end of V1-morph period
rT_B = rho_B * kT_M; % 1/d, von Bert, exponential growth rate after V1-morph period
L_b = l_b * L_m;   % cm, length at birth
L_j = l_j * L_m;     % cm, length at metamorphosis
L_i = l_i * L_m;     % cm, ultimate length
aT_b  = t_b/ kT_M;   % d, age at birth
aT_j = t_j/ kT_M;    % d, age at end V1-morph period

% tW_data
tW = tW_data(:,1);

if strcmp(dp, 'dph') == 1               % if data expressed in dph, need to add the predicted value for age at hatch
  tW=tW+aT_h;
  t_emb = tW(aT_h < tW & tW < aT_b,1);
end
t_emb = tW(aT_h < tW & tW < aT_b,1);    % available ages between hatch and birth in the dataset

if isempty(t_emb) == 0                  % if t_emb is not empty
pars = [v * TC_t; g; L_m; k_J * TC_t; kap];
 if t_emb(1) > aT_h                     % if first age value is after hatch
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

t_bj = tW(tW >= aT_b - aT_h & tW < aT_j - aT_h,1); % selects times after hatch and before birth
t_ji = tW(tW >= aT_j - aT_h,1); % selects times after metamorphosis
L_bj = L_b * exp(t_bj * rT_j/ 3); Ww_bj = L_bj.^3 * (1 + w * f);
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j - aT_h)); Ww_jm = L_jm.^3 * (1 + w * f);
EW_Data = [Ww_emb; Ww_bj; Ww_jm];                           % cm, structural length




