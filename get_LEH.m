%% script calles dget_LEH
% makes lots of plots and test out the ODE function

% load parameter values
load('results_Oncorhynchus_mykiss.mat')

cPar  = parscomp_st(par); 
vars_pull(par); vars_pull(cPar);  

TC = tempcorr(C2K(8.5), T_ref, T_A); % -, TC temperature correction

pMoA = 'control';

f = 0.7;

% initial conditions -
pars_UE0          = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
[U_E0, Lb, info]  = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve, f= 1
E0            = U_E0 * p_Am; % J, initial energy in the egg
LEH0  = [1e-5; E0; 0]; % 3-1 vector [J, cm, J] with initial conditions

age = linspace(0,400,100)';
[~, LEH] =  ode23s(@dget_LEH, age, [LEH0; 0; 0; 0],[],f, TC, par, cPar, pMoA); 

% the three zeros are for N, Lb, and Lj (see inside the ODE subfunction)

L = LEH(:,1);   % cm, structural length
E = LEH(:,2);   % J, energy in reserve
EH = LEH(:,3);   % J, energy in reserve
ER = LEH(:,4);   % J, energy in reserve
Lb = LEH(:,5);   % J, energy in reserve
Lj = LEH(:,6);   % J, energy in reserve
EW = L.^3 + w_E/ mu_E * E/ d_E; % g, wet weight


figure(1)
subplot(221)
plot(age, L/del_M,'g')
ylabel('stuct length, cm')
subplot(222)
plot(age, E,'r')
ylabel('energy in reserve, J')
subplot(223)
plot(age, E./(E_m .* L.^3),'r')
ylabel('scaled energy in reserve, -')
ylim([0 1.2])

xlim([75 400])
subplot(224)
plot(age, EW,'r')
ylabel('wet weight, g')


LEH0 = LEH(end,:); % get initial conditions
fnew = 0;
age = linspace(400,1000,100)';
[~, LEH] =  ode23s(@dget_LEH, age, LEH0,[],fnew, TC, par, cPar, pMoA); 

L = LEH(:,1);   % cm, structural length
E = LEH(:,2);   % J, energy in reserve
EH = LEH(:,3);   % J, energy in reserve
ER = LEH(:,4);   % J, energy in reserve
Lb = LEH(:,5);   % J, energy in reserve
Lj = LEH(:,6);   % J, energy in reserve
EW = L.^3 + w_E/ mu_E * E/ d_E; % g, wet weight


figure(2)
subplot(221)
plot(age, L/del_M,'g')
ylabel('stuct length, cm')
subplot(222)
plot(age, E,'r')
ylabel('energy in reserve, J')
subplot(223)
plot(age, E./(E_m .* L.^3),'r')
ylabel('scaled energy in reserve, -')
ylim([0 1.2])
subplot(224)
plot(age, EW,'r')
ylabel('wet weight, g')

