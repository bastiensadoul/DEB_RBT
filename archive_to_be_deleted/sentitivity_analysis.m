
%clear all; close all; clc

load('results_Oncorhynchus_mykiss.mat', 'par');


timevector=linspace(0, 900, 100)';
f=1;

c = parscomp_st(par); 
vars_pull(par);                     % no need to call a param using "par."
TC_test=tempcorr(C2K(20), T_ref, T_A);
vars_pull(c);

% Initial conditions
pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector, 
                                   % V_Hb: scaled maturity at birth
                                   % g: energy investment ratio (specific cost for structure)
                                   % k_J : maturity maintenance rate coeff
                                   % k_M : somatic maintenance rate coefficient
                                   % v : conductance
      
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
UT_E0 = U_E0/ TC_test; % cm * d , scaled initial reserve at T

% Age, energy in reserves and structural length at hatch and birth using U_H, scaled energy at hatch and birth (E divided by p_Am)
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);

% Get scaled age and scaled length at metamorphosis, puberty and birth
pars_tj = [g k l_T v_Hb v_Hj v_Hp];
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);

% Transform scaled age in scaled time (since t=0 : fertilization)
a_h   = aUL(2,1);                 % d, age at hatch at f and T_ref
aT_h  = a_h/ TC_test;             % d, age at hatch at f and T  
t = timevector + aT_h; % d, age since fertilization

% Correct von Bert for Temperature
kT_M  = k_M * TC_test;  % somatic maintenance rate coefficient
rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate between first feeding and end of V1-morph period
rT_B = rho_B * kT_M; % 1/d, von Bert, exponential growth rate after V1-morph period

% Transform structural length in length
L_b = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;     % cm, length at birth, metamorphosis, ultimate

% Transform scaled age to age
aT_b  = t_b/ kT_M; aT_j = t_j/ kT_M;    % d, age at birth, metamorphosis at T

% Select time points between t=0 to birth and metamorphosis
t_0b = t(t < aT_b,1);              % ages during the embryo period
t_bj = t(t >= aT_b & t < aT_j,1);  % selects times during V1-morph period
t_ji = t(t >= aT_j,1);             % selects times after metamorphosis

% Get length, scaled reserve and scaled maturity between time 0 and birth
if isempty(t_0b) == 0     % if t_emb is not empty    
t_0b = [0;t_0b];   
[a, LUH] = ode45(@dget_LUH, t_0b, [1e-4 UT_E0 0], [], kap, v * TC_test, k_J * TC_test, g, L_m);
    if length(t_0b) == 2
    LUH = LUH(end,:);
    else
    LUH = LUH(2:end,:);    
    end
    L_emb = LUH(:,1);   % cm, embryo structural length
    E_emb = LUH(:,2) * p_Am * TC_test;   % J, embryo energy in reserve
    Ww_emb = d_V * L_emb.^3 + w_E/ mu_E * E_emb; % g, embryo wet weight
else
L_emb = []; Ww_emb = [];
end

% time-length 
L_bj = L_b * exp((t_bj - aT_b) * rT_j/ 3); % cm length and weight during V1-morph period
Ww_bj = L_bj.^3 * (1 + w * f);   % g weight during V1-morph period
L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j));   % cm, length after V1-morph period
Ww_jm = L_jm.^3 * (1 + w * f); % g, weight after V1-morph period

EL = [L_emb; L_bj; L_jm]; % cm, structural length
EW = [Ww_emb; Ww_bj; Ww_jm]; % g, wet weight



% Plot 
figure(1)
plot(t,EW,'g','linewidth',2)

