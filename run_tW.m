%% Effects of a parameter change on growth and the energy budget

clear all; close all; clc


[data, auxData, txtData, weights] = mydata_BPA; % load the data matrixes of the control
load('results_Oncorhynchus_mykiss.mat') % load parameter values of the control -

c = parscomp_st(par);  vars_pull(c); vars_pull(par)

par.f        = 0.76;

TC = tempcorr(auxData.temp.tW, T_ref, T_A); % -, TC temperature correction


%% contols
auxData.pMoA = 'control';
d.tW = linspace(0,600,100)'; % = data.tW_gw150A(:,1);
t = d.tW;
controlData = predict_tW(par, d, auxData);

figure(1)
plot(t, controlData.tW, 'b','linewidth',3);
hold on;
plot(data.tW_gw150A(:,1), data.tW_gw150A(:,2), 'ro','markersize',8,'markerfacecolor','r');

%% increase or decrease p_M values
auxData.pMoA = 'p_M'; % p_M is modified
par.t_f      = 200; % day, dpf when when parameter reaches normal value again
par.delta    = 50; % factor by wich the parameter is modified at start from BPA

pD = predict_tW(par, d, auxData); % linear change in parameter value from 0 to t_f
diff = (controlData.tW - pD.tW)./controlData.tW;

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
plot(d.tW, diff,'r');
xlim([64 t(end)])
xlabel('age, dpf'); ylabel('diff to control curve, [p_M]'); set(gca,'Fontsize',12); 
figure(3); hold on;
plot(t, pT_M_t,'r')
xlabel('age, dpf'); ylabel('change in p_M'); set(gca,'Fontsize',12); 
  
 %% increase costs for growth:
 % this can only increase
auxData.pMoA = 'E_G';
par.t_f      = 42; % day, dpf when when parameter reaches normal value again
par.delta    = 5; % factor by wich the parameter is modified at start from BPA

pD = predict_tW(par, d, auxData); % linear change in parameter value from 0 to t_f
diff = (controlData.tW - pD.tW)./controlData.tW;

figure(1)
plot(t,pD.tW,'g','linewidth',2,'linestyle','--' )

% ylim([0 0.5])
% xlim([0 75])

E_G_Q = E_G * par.delta; % J/d/cm^3, p_M at start at T
E_G_t = min(E_G, E_G_Q +(E_G - E_G_Q)/ par.t_f * t); % how E_G changes

figure(4)
plot(d.tW, diff,'g');
xlabel('age, dpf'); ylabel('diff to control curve, [E_G]'); set(gca,'Fontsize',12); 
% xlim([0 75])

figure(5)
plot(t, E_G_t,'g')
xlabel('age, dpf'); ylabel('change in [E_G]'); set(gca,'Fontsize',12); 
  
  