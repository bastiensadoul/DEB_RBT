%% This script checks that predict_tW works !
% the output of predict_tW is compared with the results of using von Bert.
% growth at constant food -

% three circles on the plot represent: hatch, birth and metamorphosis
% respectively

% this scrip can be useful to look at how 'f' impacts the timing of
% metamorphosis - the period during metabolic acceleration is rather
% crucial and will have 

clear all; close all; clc


load('results_Oncorhynchus_mykiss.mat') % load parameter values of the blank taken from AmP
% re-initialise the value of f:
f     = 0.7;
par.f = f;
auxData.temp.tW = C2K(8.5); % set temperature 
c = parscomp_st(par);  vars_pull(c); vars_pull(par)

% initial conditions -
pars_UE0          = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
[U_E0, Lb, info]  = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve, f= 1


% output of predict_tW:
auxData.pMoA = 'control'; % no parameter modification

TC = tempcorr(auxData.temp.tW, T_ref, T_A);

d.tW = linspace(0,150,100)'; % = data.tW_gw150A(:,1);
prdData = predict_tW(par, d, auxData);

figure(1)
plot(d.tW, prdData.tW, 'b','linewidth',3);
hold on;

%% growth with von Bert:
% checks the exactitude of the dyanmic ODE in run_tW - notice here that
% the ab etc are computed using the maternel effect rule

pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];   
      
% Calculate Parameters
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh/TC; U_Hb/TC], [0; U_E0/TC; 1e-10], [], kap, v * TC, k_J * TC, g, L_m);
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
aT_h = aUL(2,1); % d, age at hatch at f and T_ref
E_h  = aUL(2,2) * p_Am * TC; % E, energy in reserve at f 
L_h  =  aUL(2,3); % E, energy in reserve at f 
Ww_h = L_h^3 + w_E/ mu_E/ d_E * E_h; % g, wet weight at hatch
kT_M = k_M * TC; 
rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate between first feeding and end of V1-morph period
rT_B = rho_B * kT_M; % 1/d, von Bert, exponential growth rate after V1-morph period
L_b  = l_b * L_m; L_j = l_j * L_m; L_i = l_i * c.L_m;     % cm, length at birth, metamorphosis, ultimate
L_p  = l_p * L_m;
aT_b = t_b/ kT_M; aT_j = t_j/ kT_M;  aT_p = t_p/ kT_M;     % d, age at birth, metamorphosis, puberty at T
Ww_b = L_b^3 * (1 + w * f); Ww_j = L_j^3 * (1 + w * f); Ww_p = L_p^3 * (1 + w * f);  % g, wet weight at birth, metamorphosis, puberty

t = d.tW(:,1); % d, age since fertilization
t_0b = t(t < aT_b,1);    % ages during the embryo period
t_bj = t(t >= aT_b & t < aT_j,1);    % selects times during V1-morph period
t_ji = t(t >= aT_j,1);                % selects times after metamorphosis

if isempty(t_0b) == 0     % if t_emb has values
    if t_0b(1) > 0        
    timeIn = [0;t_0b];   % append 0, if there is no zero
    else 
    timeIn = t_0b;       % otherwise don't append a zero
    end
    [a, LUH] = ode45(@dget_LUH, timeIn, [1e-4 U_E0/TC 0], [], kap, v * TC, k_J * TC, c.g, L_m);
       if  t_0b(1) > 0
           LUH = LUH(2:end,:);  % remove the appended valeu
       end
            if length(timeIn) == 2 
            LUH = LUH(end,:);
            end
    L_emb = LUH(:,1);   % cm, embryo structural length
    E_emb = LUH(:,2) * p_Am * TC;   % J, embryo energy in reserve
    Ww_emb = L_emb.^3 + w_E/ mu_E * E_emb/ d_E; % g, embryo wet weight
else
L_emb = []; Ww_emb = []; % emryo lengths and weights are empty if no embryo ages
end
% time-length 
L_bj  = L_b * exp((t_bj - aT_b) * rT_j/ 3); % cm length and weight during V1-morph period
Ww_bj = L_bj.^3 * (1 + w * f);   % g weight during V1-morph period
L_jm  = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - aT_j));   % cm, length after V1-morph period
Ww_jm = L_jm.^3 * (1 + w * f); % g, weight after V1-morph period
EL = [L_emb; L_bj; L_jm]; % cm, structural length
EW = [Ww_emb; Ww_bj; Ww_jm]; % g, wet weight

% continue plotting on the figure
figure(1)
plot(t,EW,'g','linewidth',2)
plot(aT_b, Ww_b, 'ro')
plot(aT_h, Ww_h, 'ro')
if t(end) > aT_j
plot(aT_j, Ww_j, 'ro')
end
if t(end) > aT_p
plot(aT_p, Ww_p, 'ro')
end
xlim([0 t(end)+3]); ylim([0 max(EW(end),prdData.tW(end))])
xlabel('age, dpf'); ylabel('wet weight, g');
set(gca,'Fontsize',12); 


