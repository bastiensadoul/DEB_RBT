%% Effects of a parameter change on growth and the energy budget

clear all; close all; clc


[data, auxData, txtData, weights] = mydata_BPA; % load the data matrixes of the control
load('results_Oncorhynchus_mykiss.mat') % load parameter values of the control -
c = parscomp_st(par);  vars_pull(c); vars_pull(par)

par.f        = 0.76;
pars_UE0      = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0          = initial_scaled_reserve(par.f, pars_UE0); % d.cm^2, initial scaled reserve, f= 1
E0            = U_E0 * p_Am; % J, initial energy in the egg
auxData.LEH0  = [1e-5; E0; 0]; % 3-1 vector [J, cm, J] with initial conditions
auxData.T     = C2K(8.5); % K, kelvin


%% contols
auxData.treatment = 'control';
data.tW = linspace(0,600,100)'; % = data.tW_gw150A(:,1);
prdData = predict_tW(par, data, auxData);

figure(1)
plot(data.tW, prdData.tW, 'b','linewidth',3);
hold on;
plot(data.tW_gw150A(:,1), data.tW_gw150A(:,2), 'ro','markersize',8,'markerfacecolor','r');


%% growth with von Bert:
% checks the exactitude of the dyanmic ODE in run_tW - notice here that
% the ab etc are computed using the maternel effect rule
TC = tempcorr(auxData.T, T_ref, T_A);
f = par.f;
pars_tj = [g; k; l_T; v_Hb; v_Hj; v_Hp];   
      
% Calculate Parameters
[U_H, aUL] = ode45(@dget_aul, [0; U_Hh/TC; U_Hb/TC], [0; U_E0/TC; 1e-10], [], kap, v * TC, k_J * TC, g, L_m);
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
aT_h   = aUL(2,1); % d, age at hatch at f and T_ref
E_h = aUL(2,2) * p_Am * TC; % E, energy in reserve at f 
L_h  =  aUL(2,3); % E, energy in reserve at f 
Ww_h = L_h^3 + w_E/ mu_E/ d_E * E_h; % g, wet weight at hatch
kT_M  = k_M * TC; 
rT_j = rho_j * kT_M; % 1/d, von Bert, exponential growth rate between first feeding and end of V1-morph period
rT_B = rho_B * kT_M; % 1/d, von Bert, exponential growth rate after V1-morph period
L_b = l_b * L_m; L_j = l_j * L_m; L_i = l_i * c.L_m;     % cm, length at birth, metamorphosis, ultimate
L_p = l_p * L_m;
aT_b  = t_b/ kT_M; aT_j = t_j/ kT_M;  aT_p = t_p/ kT_M;     % d, age at birth, metamorphosis, puberty at T
Ww_b = L_b^3 * (1 + w * f); Ww_j = L_j^3 * (1 + w * f); Ww_p = L_p^3 * (1 + w * f);  % g, wet weight at birth, metamorphosis, puberty

t = data.tW(:,1); % d, age since fertilization
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


%% increase or decrease p_M values
auxData.treatment = 'p_M'; % p_M is modified
par.t_f      = 200; % day, dpf when when parameter reaches normal value again
par.delta    = 50; % factor by wich the parameter is modified at start from BPA

pD = predict_tW(par, data, auxData); % linear change in parameter value from 0 to t_f
diff = (pD.tW - prdData.tW)./pD.tW * 100;

figure(1)
plot(t,pD.tW,'r','linewidth',2,'linestyle','--' )

pT_M = p_M * TC; % J/cm^3/d, vol linked som maintenance at T      
pT_M_Q = pT_M * par.delta; % J/d/cm^3, p_M at start at T
    if pT_M_Q < pT_M
    pT_M_t = min(pT_M, pT_M_Q +(pT_M - pT_M_Q)/ par.t_f * t);
    else
    pT_M_t = max(pT_M, pT_M_Q +(pT_M - pT_M_Q)/ par.t_f * t);
    end 
figure(2); hold on;
plot(data.tW, diff,'r');
xlim([64 t(end)])
xlabel('age, dpf'); ylabel('diff to control curve, [p_M]'); set(gca,'Fontsize',12); 
figure(3); hold on;
plot(t, pT_M_t,'r')
xlabel('age, dpf'); ylabel('change in p_M'); set(gca,'Fontsize',12); 
  
 %% increase costs for growth:
 % this can only increase
auxData.treatment = 'E_G';
par.t_f      = 42; % day, dpf when when parameter reaches normal value again
par.delta    = 5; % factor by wich the parameter is modified at start from BPA

pD = predict_tW(par, data, auxData); % linear change in parameter value from 0 to t_f
diff = (pD.tW - prdData.tW)./pD.tW * 100;

figure(1)
plot(t,pD.tW,'g','linewidth',2,'linestyle','--' )

% ylim([0 0.5])
% xlim([0 75])

E_G_Q = E_G * par.delta; % J/d/cm^3, p_M at start at T
E_G_t = max(E_G, E_G_Q +(E_G - E_G_Q)/ par.t_f * t); % how E_G changes

figure(4)
plot(data.tW, diff,'g');
xlabel('age, dpf'); ylabel('diff to control curve, [E_G]'); set(gca,'Fontsize',12); 
% xlim([0 75])

figure(5)
plot(t, E_G_t,'g')
xlabel('age, dpf'); ylabel('change in [E_G]'); set(gca,'Fontsize',12); 
  
  