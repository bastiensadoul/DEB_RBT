% close all;  
% clear all;  clc
% 
% % Calls AmPtool - pars has parameter values in columns and each species in
% % lines
% pars = read_allStat('p_Am','v', 'E_Hb', 'E_Hj', 'E_Hp','s_M', 'L_i','p_M',  'kap', 'E_m', ...
%     'h_a', 'a_m', 'L_m', 'T_A', 'c_T', 'J_Oi', 'Wd_i', 'R_i', 'N_i', 'Wd_0');
% 
% % unpack parameters & other statistics
% % keep in mind that parameters are at T_ref 
% % statistics (with time in their dimension) are at T_typical
% p_Am = pars(:,1);
% v    = pars(:,2);
% E_Hb = pars(:,3);
% E_Hj = pars(:,4);
% E_Hp = pars(:,5);
% s_M  = pars(:,6);
% L_i  = pars(:,7);
% p_M  = pars(:,8);
% kap  = pars(:,9);
% E_m  = pars(:,10);
% h_a  = pars(:,11);
% a_m  = pars(:,12);
% L_m  = pars(:,13);
% T_A = pars(:,14);
% c_T = pars(:,15);
% 
% J_Oi = pars(:,16); JT_Oi = J_Oi./ c_T; % correct to reference temp
% Wd_i = pars(:,17);
% R_i = pars(:,18);
% N_i = pars(:,19);
% Wd_0 = pars(:,20);
% 
% [selM, srcM] = select_01('Animalia','Actinopterygii'); % indices for Bony fish
% [selB, srcB] = select_01('Animalia','Salmoniformes'); % indices for salmoniformes
% % [selA, srcA] = select_01('Animalia','Salmoniformes');  % indices for Arctica Islandica
% [selT, srcT] = select_01('Animalia','Oncorhynchus_mykiss'); % indices for Tridacna gigas
% 

%% Max respiration (at T_ref) as function of Wd_i (ultimate dry weight) 

figure(); hold on
plot(log10(Wd_i(selM)), log10(JT_Oi(selM)), 'k.', 'markersize',8)
plot(log10(Wd_i(selB)), log10(JT_Oi(selB)),   'bo', 'markersize',6,'markerfacecolor','b')
plot(log10(Wd_i(selT)), log10(JT_Oi(selT)),'bs', 'markersize',12,'markerfacecolor','b')
xlabel('Max ultimate dry weight, _{10}log Wd_i, g')
ylabel('Max respiration, _{10}log J_{Oi} mol/d')
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Wdi_JTOi.png')


%% Max respiration (at T_ref) as function of Wd_i (ultimate dry weight) 
figure(); hold on
plot(log10(Wd_i(selM)), log10(JT_Oi(selM)./ Wd_i(selM)), 'k.', 'markersize',8)
plot(log10(Wd_i(selB)), log10(JT_Oi(selB)./ Wd_i(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
plot(log10(Wd_i(selT)), log10(JT_Oi(selT)./ Wd_i(selT)),'bs', 'markersize',12,'markerfacecolor','b')
xlabel('ultimate dry weight, _{10}log Wd_i, g')
ylabel('max spec. O_2, _{10}log J_{Oi} mol/g.d')
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Wdi_MJTOi.png')

%% Mean life span (at species-specific typical temp) as function of p_M (given at ref. temp.)
figure(); hold on
plot(log10(p_M(selM)), log10(a_m(selM)), 'k.', 'markersize',8)
plot(log10(p_M(selB)), log10(a_m(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
plot(log10(p_M(selT)), log10(a_m(selT)),'bs', 'markersize',12,'markerfacecolor','b')
xlabel('Max somat. maint. rate, _{10}log p_M, cm')
ylabel('Mean life span, _{10}log a_m')
xlim([min(log10(p_M(selM))) max(log10(p_M(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/pM_am.png')


%% mean life span at T_ref as function of p_M (at reference temperature)
figure(); hold on
plot(log10(p_M(selM)), log10(a_m(selM) .* c_T(selM)), 'k.', 'markersize',8)
plot(log10(p_M(selB)), log10(a_m(selB).* c_T(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
plot(log10(p_M(selT)), log10(a_m(selT).* c_T(selT)),'bs', 'markersize',12,'markerfacecolor','b')
xlabel('Max somat. maint. rate, _{10}log p_M, cm')
ylabel('Mean life span, _{10}log a_m')
xlim([min(log10(p_M(selM))) max(log10(p_M(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/pM_am_Tref.png')

%% Weibull ageing acceleration as function of ultimate structural length
figure(); hold on
plot(log10(L_i(selM)), log10(h_a(selM)), 'k.', 'markersize',8)
plot(log10(L_i(selB)), log10(h_a(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
plot(log10(L_i(selT)), log10(h_a(selT)),'bs', 'markersize',12,'markerfacecolor','b')
xlabel('Ultimate struct. length, _{10}log L_i, cm')
ylabel('Weibull ageing acceleration, _{10}log h_a')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Li_ha.png')

%%
figure(); hold on
plot(log10(L_i(selM)),log10(p_M(selM)), 'k.', 'markersize',8)
plot(log10(L_i(selB)),log10(p_M(selB)),'bo', 'markersize',6,'markerfacecolor','b')
plot(log10(L_i(selT)), log10(p_M(selT)),'bs', 'markersize',12,'markerfacecolor','b')
plot([min(log10(L_i(selM))) max(log10(L_i(selM)))],[log10(18) log10(18)], 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
xlabel('Ultimate struct. length, _{10}log L_i, cm')
ylabel('Somatic maint. _{10}log [p_M], J/cm^3/d')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Li_pM.png')

%%
figure(); hold on
plot(log10(L_i(selM)), log10(v(selM).* s_M(selM)), 'k.', 'markersize',8)
plot(log10(L_i(selB)), log10(v(selB).* s_M(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
% plot(log10(L_i(selBivalvia)), log10(v(selBivalvia)), 'bo', 'markersize',6)
% plot(log10(L_i(selA)), log10(v(selA).* s_M(selA)),'gs', 'markersize',12,'markerfacecolor','g')
% plot(log10(L_i(selA)), log10(v(selA)),'gs', 'markersize',12)
plot(log10(L_i(selT)), log10(v(selT).* s_M(selT)),'bs', 'markersize',12,'markerfacecolor','b')
plot(log10(L_i(selT)), log10(v(selT)),'bs', 'markersize',12)
plot([min(log10(L_i(selM))) max(log10(L_i(selM)))],[log10(0.02) log10(0.02)], 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
xlabel('Ultimate struct. length, _{10}log L_i, cm')
ylabel('Energy conductance _{10}log v, cm/d')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Li_v.png')

%% Assimilation as function of maximum size
figure(); hold on
plot(log10(L_i(selM)), log10(p_Am(selM).* s_M(selM)), 'k.', 'markersize',8)
plot(log10(L_i(selB)), log10(p_Am(selB).* s_M(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
% plot(log10(L_i(selA)), log10(p_Am(selA).* s_M(selA)), 'gs', 'markersize',12,'markerfacecolor','g')
% plot(log10(L_i(selA)), log10(p_Am(selA)), 'gs', 'markersize',12)
plot(log10(L_i(selT)), log10(p_Am(selT).* s_M(selT)), 'bs', 'markersize',12,'markerfacecolor','b')
plot(log10(L_i(selT)), log10(p_Am(selT)), 'bs', 'markersize',12)
% plot(log10(L_m(selMollusca)), log10(22.5 * L_m(selMollusca)), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
plot(log10(L_i(selM)), log10(22.5 * L_i(selM)), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
xlabel('Ultimate struct. length, _{10}log L_i, cm')
ylabel('Max assim rate _{10}log {p_{Am}}, J/cm^2/d')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Li_pAm.png')

%% E_H^b as function of L_inf
figure(); hold on
plot(log10(L_i(selM)), log10(E_Hb(selM)), 'k.', 'markersize',8)
plot(log10(L_i(selB)), log10(E_Hb(selB)),'bo', 'markersize',6,'markerfacecolor','b')
% plot(log10(L_i(selA)), log10(E_Hb(selA).* s_M(selA)), 'gs', 'markersize',12,'markerfacecolor','g')
plot(log10(L_i(selT)), log10(E_Hb(selT).* s_M(selT)), 'bs', 'markersize',12,'markerfacecolor','b')
plot(log10(L_i(selM)),log10(0.275 * L_i(selM).^3), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
xlabel('Ultimate struct. length, _{10}log L_i, cm')
ylabel('Mat. level at birth _{10}log E_H^b, J')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Li_EHb.png')


%% E_H^j as function of L_inf
figure(); hold on
plot(log10(L_i(selM)), log10(E_Hj(selM)), 'k.', 'markersize',8)
plot(log10(L_i(selB)),log10(E_Hj(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
% plot(log10(L_i(selA)),log10(E_Hj(selA)), 'gs', 'markersize',12,'markerfacecolor','g')
plot(log10(L_i(selT)),log10(E_Hj(selT)), 'bs', 'markersize',12,'markerfacecolor','b')
plot(log10(L_i(selM)),log10(0.275 * L_i(selM).^3), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
xlabel('Ultimate struct. length, _{10}log L_i, cm')
ylabel('Mat. level at metam. _{10}log E_H^j, J')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Li_EHj.png')

%% E_H^p as function of L_inf
figure(); hold on
plot(log10(L_i(selM)), log10(E_Hp(selM)), 'k.', 'markersize', 8)
plot(log10(L_i(selB)),log10(E_Hp(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
% plot(log10(L_i(selA)),log10(E_Hp(selA)), 'gs', 'markersize',12,'markerfacecolor','g')
plot(log10(L_i(selT)),log10(E_Hp(selT)), 'bs', 'markersize',12,'markerfacecolor','b')
plot(log10(L_i(selM)),log10(166 * L_i(selM).^3), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
xlabel('Ultimate struct. length, _{10}log L_i, cm')
ylabel('Mat. level at puberty _{10}log E_H^j, J')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Li_EHp.png')



%% E_m as function of L_inf
figure(); hold on
plot(log10(L_i(selM)), log10(E_m(selM)), 'k.', 'markersize', 8)
plot(log10(L_i(selB)), log10(E_m(selB)), 'bo', 'markersize',6,'markerfacecolor','b')
% plot(log10(L_i(selA)), log10(E_m(selA)), 'gs', 'markersize',12,'markerfacecolor','g')
plot(log10(L_i(selT)), log10(E_m(selT)), 'bs', 'markersize',12,'markerfacecolor','b')
plot(log10(L_i(selM)), log10(22.5 * L_i(selM)/ 0.02), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
xlabel('Ultimate struct. length, _{10}log L_i, cm')
ylabel('Max. res. dens. _{10}log [E_m], J/cm^3')
xlim([min(log10(L_i(selM))) max(log10(L_i(selM)))]  )
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/Li_Em.png')


%%
pMsurv = surv(p_M);
xsurv = pMsurv(:,1); ysurv = pMsurv(:,2);
pMsurv_B = surv(p_M(selB));
% pMsurv_A = surv(p_M(selA));

figure();hold on
plot(log10(xsurv),ysurv, 'b','linewidth',2)
plot(log10(pMsurv_B(:,1)),pMsurv_B(:,2), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
% plot([log10(p_M(selA)) log10(p_M(selA))],[0 1], 'k--')
plot([-0.5 2 ],[0.5 0.5], 'k')
xlabel('log_{10} [p_M]')
ylabel('survivor function, -')
set(gca, 'FontSize', 15, 'Box', 'on')


%%
vsurv = surv(v(selM)); % survival and embryo values for mollusca
xsurv = vsurv(:,1); ysurv = vsurv(:,2);
vsurv_B = surv(v(selB)); % survival and embryo values for bivalvia
% vsurv_A = surv(v(selA));

vsurv_adM = surv(v(selM) .* s_M(selM)); % survival and post-metam values for mollusca
vsurv_adB = surv(v(selB) .* s_M(selB)); % survival and post-metam values for bivalvia

figure();hold on
plot(log10(xsurv),ysurv, 'b','linewidth',2)
plot(log10(vsurv_adM(:,1)),vsurv_adM(:,2), 'b--','linewidth',2) % post-metam values
plot(log10(vsurv_B(:,1)),vsurv_B(:,2), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
plot(log10(vsurv_adB(:,1)),vsurv_adB(:,2), 'r--','linewidth',2) % post-metam values

% plot([log10(v(selA)) log10(v(selA))],[0 1], 'k--')
% plot([-4 -1 ],[0.5 0.5], 'k')
xlabel('log_{10} v ')
ylabel('survivor function, -')
set(gca, 'FontSize', 15, 'Box', 'on')
saveas(gca,'figs/survV.png')
% 
%%
kapsurv = surv(kap);
xsurv = kapsurv(:,1); ysurv = kapsurv(:,2);
kapsurv_B = surv(kap(selB));
% kapsurv_A = surv(kap(selA));

figure();hold on
plot(xsurv,ysurv, 'b','linewidth',2)
plot(kapsurv_B(:,1),kapsurv_B(:,2), 'color',[0.7 0.7 0.7], 'linewidth',1, 'linestyle', '--')
% plot([kap(selA) kap(selA)],[0 1], 'k--')
plot([0 1 ],[0.5 0.5], 'k')
xlabel('kappa, -')
ylabel('survivor function, -')
set(gca, 'FontSize', 15, 'Box', 'on')

