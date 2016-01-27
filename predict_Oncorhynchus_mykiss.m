function [prdData, info] = predict_Oncorhynchus_mykiss(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
   
  
  % customized filters for allowable parameters of the standard DEB model (std)
  % for other models consult the appropriate filter function.
  filterChecks = k * v_Hp >= f_tW^3 || ...         % constraint required for reaching puberty with f_tWw
                 ~reach_birth(g, k, v_Hb, f_tW) || ...   % constraint required for reaching birth with f_tWw
                 k * v_Hp >= f_LW^3 ||    ~reach_birth(g, k, v_Hb, f_LW) || ...   
                 k * v_Hp >= f_tWL^3 ||   ~reach_birth(g, k, v_Hb, f_tWL) || ...   
                 E_Hh >= E_Hb || ...
                 f_tWL>1 || f_tW > 1 ;   
             
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  end  
  

%% compute temperature correction factors
  TC_ah5 = tempcorr(temp.ah_5, T_ref, T_A);
  TC_ah = tempcorr(temp.ah, T_ref, T_A);  
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  TC_tW = tempcorr(temp.tW, T_ref, T_A);
  TC_tWw = tempcorr(temp.tWw, T_ref, T_A);
  TC_Tah = tempcorr(Tah(:,1), T_ref, T_A);
  
  

  % life cycle
  pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

  % initial 
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve  
  E_0 = U_E0 * p_Am;     % J, initial reserve
  Wd0 = 1e3 * E_0 * w_E/ mu_E ; % mg, egg dry weight 
    
  % hatch
  [U_H aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  aT_h5 = aUL(2,1)/ TC_ah5;              % d, age at hatch at f and T
  aT_h = aUL(2,1)/ TC_ah;              % d, age at hatch at f and T
  L_h = aUL(2,3); % cm, strucural length at hatch
  E_h = aUL(2,2) * p_Am; % J, energy in reserves at hatch
  Wdh = 1e3 * (d_V * L_h^3 + w_E/ mu_E * E_h); % mg, dry weight at hatch
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth of foetus  at f = 1
  aT_b5 = t_b/ k_M/ TC_ah5;           % d, age at birth of foetus at f and T
  aT_b = t_b/ k_M/ TC_ah;           % d, age at birth of foetus at f and T
  Wdb = 1e3 * d_V * L_b^3 * (1 + f * w); % mg, dry weight at birth
     
  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, total length at puberty at f
  aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at f and T
  Ww_p = 1e3 * d_V * L_p^3 * (1 + f * w); % g, wet weight at puberty at f

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
  prdData.ah_5 = aT_h5;
  prdData.ah = aT_h;
  prdData.ab_5 = aT_b5;
  prdData.ab = aT_b;  
  prdData.ap = aT_p;
  prdData.am = aT_m;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wi = Ww_i;
  prdData.Ri = RT_i;
  prdData.Wd0 = Wd0;
  prdData.Wdh = Wdh;
  prdData.Wdb = Wdb;
  prdData.Wp =  Ww_p;

  %% uni-variate data
  
  % t-Ww-data
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tW);
  rT_B = TC_tW * rho_B * k_M; rT_j = TC_tW * rho_j * k_M; % 1/d, von Bert, exponential growth rate
  aT_b = t_b/ k_M/ TC_tW; aT_j = t_j/ k_M/ TC_tW;
  L_b = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;
  L_0 = (W0.tW/ (1 + f_tW * w))^(1/3); % cm, structural length at t = 0 
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
  EW = L.^3 * (1 + f_tW * w);
    
  % L-Ww,  
  LWw = (LWw(:,1) * del_M).^3 * (1 + f_LW * w); % g, wet mass
  
  %-------------------------------------
% Build in pred_pars, all the parameters needed in predict_length and predict_weight 
%(parameters common to all studies, independent from f and T

  % t-L and t-Ww, DaviKenn2014
  [EWw, EL]  = predict_WL(f_tWL, TC_tWw, tL, 'dph', par, cPar);
  EL  = EL/ del_M; % cm, structural length
  EWw = predict_WL (f_tWL, TC_tWw, tWw, 'dph', par, cPar); 

  % T-ah, Vels1987:
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
  [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  a_h = aUL(2,1);                 % d, age at hatch at f and T_ref
  EaT_h =  a_h ./ TC_Tah;              % d, age at hatch at f and T
  
  
  
  
  % pack to output
  prdData.tW  = EW;
  prdData.LWw = LWw;
  prdData.tL  = EL;
  prdData.tWw = EWw;
  prdData.Tah = EaT_h;
     
  
  
%% %%%%%%%%------------------- get length   ------------------------------------------------------
  
  function [EW, EL] = predict_WL(f, TC, tL_data, dp, p,c)                    
      %tL_data (1st column is time in dpf, second is length)                                                                                 
      % dp, string 'dpf' for days post fertilization     or 'dph' for days post hatch 
  
pars_tj = [c.g c.k c.l_T c.v_Hb c.v_Hj c.v_Hp];   
pars_UE0 = [c.V_Hb; c.g; p.k_J; c.k_M; p.v]; % compose parameter vector
      
% Calculate Parameters
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
[U_H, aUL] = ode45(@dget_aul, [0; c.U_Hh; c.U_Hb], [0 U_E0 1e-10], [], p.kap, p.v, p.k_J, c.g, c.L_m);
a_h   = aUL(2,1);                 % d, age at hatch at f and T_ref
U_Eh  = aUL(2,2);                 % d cm^2, scaled reserve at hatch at T_ref
UT_Eh = U_Eh/ TC;               % d cm^2, scaled reserve at hatch at T
UT_Hh = p.E_Hh/ (c.p_Am * TC);      % d cm^2, scaled maturity at hatch at T
L_h   = aUL(2,3);                 % cm, structural length at hatch
aT_h  = a_h/ TC;                % d, age at hatch at f and T  
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
kT_M  = c.k_M * TC; 
rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate between first feeding and end of V1-morph period
rT_B = rho_B * kT_M; % 1/d, von Bert, exponential growth rate after V1-morph period
L_b = l_b * c.L_m;   % cm, length at birth
L_j = l_j * c.L_m;     % cm, length at metamorphosis
L_i = l_i * c.L_m;     % cm, ultimate length
aT_b  = t_b/ kT_M;   % d, age at birth
aT_j = t_j/ kT_M;    % d, age at end V1-morph period

% tL_data
tL = tL_data(:,1);

if strcmp(dp, 'dph') == 1                            % if data expressed in dph, need to add the predicted value for age at hatch
  tL=tL+aT_h;
end
t_emb = tL(aT_h < tL & tL < aT_b,1);    % available ages between hatch and birth in the dataset

if isempty(t_emb) == 0     % if t_emb is not empty
pars = [p.v * TC; c.g; c.L_m; p.k_J * TC; p.kap];
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
 E_emb = ULH(:,1) * c.p_Am * TC;   % J, embryo energy in reserve
 Ww_emb = d_V * L_emb.^3 + w_E/ mu_E * E_emb; % g, embryo wet weight
else
L_emb = []; Ww_emb = [];
end

t_bj = tL(tL >= aT_b & tL < aT_j,1);        % selects times during V1-morph period
t_ji = tL(tL >= aT_j,1);                    % selects times after metamorphosis
% time-length 
L_bj = L_b * exp((t_bj - aT_b) * rT_j/ 3); Ww_bj = L_bj.^3 * (1 + c.w * f);   % cm,g, length and weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j));   % cm, length after V1-morph period
Ww_jm = L_jm.^3 * (1 + c.w * f); % g, weight after V1-morph period

EL = [L_emb; L_bj; L_jm]; % cm, structural length
EW = [Ww_emb; Ww_bj; Ww_jm]; % g, wet weight



 
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
  


  