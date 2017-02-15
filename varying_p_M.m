clear all; close all; clc
[~, ~, metaDataAux, ~, ~] = feval('mydata_Oncorhynchus_mykiss');  
[par, metaPar, txtPar] = feval('pars_init_Oncorhynchus_mykiss', metaDataAux);

% Data
data=[1:1000]';
f=1;
TC=8.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%% p_M following Spring and Damper model from Sadoul et al. 2015 (Plos One).

ks=1; % empirical
cs=50; % empirical
Fpert=50; % Force of the perturbation at the beginning of the perturbation
param=[ks;cs];
tminmax=[0;40]; % begining of the perturbation and end of the perturbation
ic = [0;0]; % initial conditions of p_M (vector of speed and accceleration)
tspan = [0:1000]; % solve between the two values
M = 0; % there is a mass on the damper (makes oscillations) / problem with M>0


% Call an ODE solver and plot the results.
[tvec, yvec] = ode45(@spring_damper_model, tspan, ic , [], param, Fpert, M, tminmax);

p_M_BPA = yvec(:,1) + 133 ; % 133 is the p_M of the control fish

figure(1);
plot(tvec,p_M_BPA), title('p_M for BPA condition over Time');

%%%%%%%%%%%%%%%%%%%%%%%%% Predict using p_M varying over time.
% 
% out=predict_forvarying_p_M(f, TC, data, par, p_M_BPA);
% 
% figure(2);
% plot(data,out), title('predict over time');




