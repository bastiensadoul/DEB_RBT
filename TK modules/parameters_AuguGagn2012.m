% ---------------- parameters_AuguGagn2011b ------------------------------
% parameters for zebrafish Danio rerio established in Augustine et al 2011 
% Comp Biochem Physiol A
%  Executed in functions:
%  - fit_adult_cotrol;m: [W Lobs R] = fit_adult_control(par, aW, aL, aR)

% created october 21 2011
% last modified: 

% ------------------------------------------------------------------------
%-------------------------------------------------------------------------
% parameters (AuguGagn2011a):
T_A     = 3000;        % K, Arhennius temperature
del_Y   = .1054*1.07;       % -, shape coefficient to convert vol-length to physical length
del_M   = .1325;       % -, shape coefficient to convert vol-length to physical length
v       = 0.0278;      % cm/d, energy conductance
kap     = 0.437;       % -, alloaction fraction to soma = growth + somatic maintenance
kap_R   = 0.95;        % -, reproduction efficiency
p_M     = 500.9;       % J/d.cm^3, [p_M] vol-specific somatic maintenance
k_J     = 0.0166;      % 1/d, < k_M = p_M/E_G, maturity maint rate coefficient
E_G     = 4652;        % J/cm^3, [E_G], spec cost for structure
p_Am    = 246.3;       % J/cm^3/d, maximum surface area-specific assimilation 
E_Hb    = 0.54;        % J, E_H^b cumulated energy invested in maturation at birth
E_Hp    = 2062;        % J, E_H^p cumulated energy invested in maturation at puberty  
E_Hj    = 19.66;       % J, E_H^j cumulated energy invest in maturation at metamorphosis

% parameters that link moles to grams (wet weight), volumes and energy
d_V   = 0.2; % g/cm^3, specific densities of structure
mu_E  = 500000;                   % J/mol, chemical potential of reserve
mu_V  = 500000;                   % J/ mol, chemical potential of strucure
w_E   = 23.9;                     % g/mol, mol-weights for reserve
w_V   = 23.9;                     % g/ mol, mol-weights for structure
M_V   = d_V/ w_V;                 % mol/ cm^3, specific structural mass
E_m   = p_Am/ v ;                 % J/ cm^3 [E_m], maximum reserve density
M_Em  = E_m/ mu_E;                % mol/ cm^3 [M_Em], maximum reserve desity in c-moles

% compound parameters:
% used for DEBtool routines (initial_scaled_reserve)
g   = E_G * v/ (kap * p_Am);   % -, energy investment ratio
k_M = p_M/ E_G;                % 1/d, somatic maintenance rate coefficient
k   = k_J/ k_M;                % -, maintenance ratio
L_m = kap * p_Am/ p_M;         % cm, ultimate structural length
kap_G = mu_V * d_V * 1/E_G * 1/ w_V; % -, growth conversion efficiency

U_Hb = E_Hb/ p_Am;                             % cm^2 d, scaled maturity at birth
V_Hb = U_Hb/ (1 - kap);                        % cm^2 d, scaled maturity at birth
v_Hb = V_Hb * g^2 * k_M^3/ v^2;                %  -, scaled maturity density at birth
U_Hj = E_Hj/ p_Am;                             % cm^2 d, scaled maturity at metamorphosis
V_Hj = U_Hj/ (1 - kap);                        % cm^2 d, scaled maturity at metamorphosis
v_Hj = V_Hj * g^2 * k_M^3/ v^2;                % -, scaled maturity density at metamorphosis 


% -------------------------------------------------------------------------
