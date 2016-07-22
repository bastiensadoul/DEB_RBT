% predicts Weigth as funciton of age 

function [prdData, info] = predict_tWw(p, data, auxData)
  
 % unpack par, data, auxData

 c = parscomp_st(p);          

 
 filterChecks = c.k * c.v_Hp >= p.f ||  ~reach_birth(c.g, c.k, c.v_Hb, p.f);    
         
  if filterChecks  
    info = 0;
    prdData = {};
    return;
  else
    info = 1;
  end   
 
 
%% compute temperature correction factors
 TC  = tempcorr(auxData.temp.tWw, p.T_ref, p.T_A);
   
params = [c.g, c.k, c.l_T, c.v_Hb, c.v_Hj, c.v_Hp];

% set the initial conditions depending on whether the experiment is days
% post fertilization (dpf) or days since first feeding which assume is
% birth 'dpb'

if strcmp(auxData.t0.tWw, 'dpf') == 1

pars_UE0 = [c.V_Hb; c.g; p.k_J; c.k_M; p.v]; % compose parameter vector
U_E0 = initial_scaled_reserve(1, pars_UE0); % d.cm^2, initial scaled reserve  
E_0 = U_E0 * c.p_Am;     % J, initial reserve
ELH_init = [E_0; 1e-4; 0];
time = [0; data.tWw(:,1)];
[lj, lp, lb, info] = get_lj(params, p.f); L_b = lb * c.L_m; L_j = lj * c.L_m;
[a, ELH] = ode45(@dget_ELH_pj, time, ELH_init, [], L_b, L_j, c.L_m, c.p_Am, p.v, c.g, p.k_J, p.kap, p.f, p.E_Hb, p.E_Hj, TC);
ELH(1,:) = []; E = ELH(:,1); L = ELH(:,2); EWw = L.^3 + E * c.w_E/ p.mu_E/ p.d_E;  % g, wet weight

elseif strcmp(auxData.t0.tWw, 'dpb') == 1
[lb, info] = get_lb([c.g, c.k, c.v_Hb], 1); Lb = lb * c.L_m; 
ELH_init = [c.E_m * Lb^3; Lb; p.E_Hb]; 
if data.tWw(1,1) == 0
time = data.tWw(:,1);
else
 time = [0; data.tWw(:,1)];
end
[lj, lp, lb, info] = get_lj(params, p.f); L_b = lb * c.L_m; L_j = lj * c.L_m;
[a, ELH] = ode45(@dget_ELH_pj, time, ELH_init, [], L_b, L_j, c.L_m, c.p_Am, p.v, c.g, p.k_J, p.kap, p.f, p.E_Hb, p.E_Hj, TC);

if data.tWw(1,1) ~= 0
ELH(1,:) = [];
end

E   = ELH(:,1); L = ELH(:,2); 
EWw = L.^3 + E * c.w_E/ p.mu_E/ p.d_E;  % g, wet weight

end

   
prdData.tWw = EWw; 

function dELH = dget_ELH_pj(t, ELH, Lb, Lj, Lm, p_Am, v, g, kJ, kap, f, Hb, Hj, TC)
%     dELH = dget_ELH_pj(t, ELH, Lb, Lj, Lm, p_Am, v, g, kJ, kap, f, Hb, Hj, tT, T_A, T_ref)
  %  change in state variables during juvenile stage
  %  dELH = dget_ELH_p_pj(t, ELH)
  %  ELH: 3-vector
  %   E: reserve E
  %   L: structural length
  %   H: maturity E_H
  %  dELH: change in reserve, length, scaled maturity
  
 
  %  unpack variables
  E = ELH(1); L = ELH(2); H = ELH(3);
  
%   TC = tempcorr(C2K(spline1(t, tT)), T_ref, T_A);
  vT = v * TC; kT_J = kJ * TC; pT_Am = p_Am * TC;
 
  if H < Hb 
    s = 1; % -, multiplication factor for v and {p_Am}
  elseif H < Hj
    s = L/ Lb;
  else
    s = Lj/ Lb;
  end
  e = vT * E/ L^3/ pT_Am; % -, scaled reserve density; 
  rT = s * vT * (e/ L - 1/ Lm/ s)/ (e + g); % 1/d, spec growth rate
  pT_C = E * (s * vT/ L - rT); % cm^2, scaled mobilisation
  
  % generate dH/dt, dE/dt, dL/dt
  dH = (1 - kap) * pT_C - kT_J * H;
  dE = (L > Lb) * s * pT_Am * f * L^2 * (H>= Hb) - pT_C;
  dL = rT * L/3;

  % pack derivatives
  dELH = [dE; dL; dH];

  
