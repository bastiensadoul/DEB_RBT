load('results_Oncorhynchus_mykiss.mat')


[stat1 txtStat] = statistics_st('abj', par, C2K(8.5), 0.7) ;


par.v = par.v * 0.6;

[stat2 txtStat] = statistics_st('abj', par, C2K(8.5), 0.7) ;

p = par; c = parscomp_st(p);
UT_E0 = stat2.U_E0/ stat2.c_T; % cm^2 * d, scaled energy in egg at T and f
vT = p.v * stat2.c_T; kT_J = p.k_J * stat2.c_T; pT_Am = c.p_Am * stat2.c_T;
LUH = [1e-9; UT_E0; 0]; % initial parameters
[t, LUH] = ode23s(@dget_LUH, [0;66], LUH, [], p.kap, vT, kT_J, c.g, c.L_m);

L = LUH(:,1); 
U = LUH(:,2); % cm and J, structural length and energy in reserve
E = U * pT_Am;


stat2.Wd_b 
Wd_66 = p.d_V * L(end)^3 + c.w_E/ p.mu_E * E(end)


