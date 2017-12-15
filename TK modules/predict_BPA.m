function prdData = predict_BPA(par, data, auxData)


vars_pull(par); vars_pull(data); vars_pull(auxData);

% L:  length 0-9 dpf, c_d = 0, 84 and 1050.4 nmol/l, f=0, 25°C 
% Wd: dry mass  0-9 dpf, c_d = 0, 84 and 1050.4 nmol/l, f=0, 25°C
% MQ: internal concentration c_d = 0, 84 and 1050.4 nmol/l, f=0, 25°C
% global param

% ------------------------------------------------------------------------
%-------------------------------------------------------------------------

% parameters_AuguGagn2012

% % ---------------------- unpack par:
% NEC     = par(1)  ;       % mol U cm^-3 structure, no effect internal concentration
% M_QT    = par(2);         % mol U cm^-3 structure, tolerance concentration
% P_Vd    = par(3);         % liter water cm^-3 V, partition coefficient water/structure
% P_EV    = par(4) ;        % -, partition coefficient reserve/structure
% k_e     = par(5) ;        % d^-1, elimination rate
% f_SBemb = par(6);         %-6- -, scaled functional response, BourSimo2008



% --------------------- ELIMINATION ---------------------------------------
c_d = 0;                              % nmol/l, environmental concentration
Wd0 = w_E/ mu_E * E0;

% ----------------------------------------------- 
MQ  = Bour09F1(1,1);                     % nmol/g, amount in egg (observation)
CI = [1e-4,E0, 0,0,0,MQ * Wd0];
[~, LEHRU] = ode45(@ode_LEHRU_bj_absolute,a_SB,CI,[],p,p_tox,f_SBemb,c_d, param,0); 
L  = LEHRU(:,1);                        % cm, structural length
E  = LEHRU(:,2);                        % J, energy in reserve and reproduction buffer
Wd = d_V * L.^3 + w_E/ mu_E * E;        % g, dry weight 
eMQ = LEHRU(:,6);                        % nmol, amount of uranium per unit structure
eMQ1  = eMQ./ Wd ;                       % nmol/ g, observed internal concentration 


% ----------------------------------------------- 
MQ  = Bour09F1(2,1);                     % nmol/g, amount in egg (observation)
CI = [1e-4,E0, 0,0,0,MQ * Wd0];
[~, LEHRU] = ode45(@ode_LEHRU_bj_absolute,a_SB,CI,[],p,p_tox,f_SBemb,c_d, param,0); 
L  = LEHRU(:,1);                        % cm, structural length
E  = LEHRU(:,2);                        % J, energy in reserve and reproduction buffer
Wd = d_V * L.^3 + w_E/ mu_E * E;        % g, dry weight 
eMQ = LEHRU(:,6);                        % nmol, amount of uranium per unit structure
eMQ2  = eMQ./ Wd ;                        % nmol/ g, observed internal concentration 

% ----------------------------------------------- 
MQ  = Bour09F1(3,1);                     % nmol/g, amount in egg (observation)
CI = [1e-4,E0, 0,0,0,MQ * Wd0 ];
[~, LEHRU] = ode45(@ode_LEHRU_bj_absolute,a_SB,CI,[],p,p_tox,f_SBemb,c_d, param,0); 
L  = LEHRU(:,1);                        % cm, structural length
E  = LEHRU(:,2);                        % J, energy in reserve and reproduction buffer
Wd = d_V * L.^3 + w_E/ mu_E * E;        % g, dry weight 
eMQ = LEHRU(:,6);                        % nmol, amount of uranium per unit structure
eMQ3  = eMQ./ Wd ;               % nmol/ g, observed internal concentration 

% -------- output eBour09F1: 
eBour09F1 = [...
    eMQ1(end) % nmol/g at last day of observation
    eMQ2(end) % nmol/g at last day of observation
    eMQ3(end) % nmol/g at last day of observatoin
    ];

 


%% Subfunction:

function dLEHRU = ode_LEHRU_bj_absolute(t, LEHRU, par, f, c_d)
% dLEHRU = ode_LEHRU(~, LEHRU,p,p_tox,f, c_d,Lbj,lab, param)
% f: scaled functional response, c_d: mol/l environmental concentration
% p     = [v * TC; k_J * TC;  kap; E_Hb; E_Hj; E_Hp; p_Am * TC ; E_G * TC; p_M * TC; mu_E; d_V; w_V; w_E; del_M; del_Y; E_ERs; kap_R; kap_G];
% p_tox = [P_EV; P_Vd; k_e; M_QT; NEC] 
% Lbj = [L_b; L_j] if starting after a=0  or Lbj = [] if starting with a=0;
% param: str 'p_M', 'E_G' or 'p_Am' or 'control'

% !!!!!!!!!!!!!! does not have batch preparation !!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!! works in absolute amount !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    persistent L_b L_j 

% v    = p(1) ; kJ   = p(2);  kap  = p(3); EHb  = p(4); EHj  = p(5); EHp  = p(6); p_Am = p(7);
% E_G  = p(8); p_M  = p(9); mu_E = p(10); d_V = p(11); w_V = p(12);  kapR = p(17);kapG = p(18);
% 
% P_EV = p_tox(1); P_Vd = p_tox(2); k_e = p_tox(3); M_QT = p_tox(4); NEC = p_tox(5); 

% --- COMPOUND PARAMETERS -----------------------------------------------
% !!! make sur scaling is done before the parameter modification !!!!!!!!!
L_m  = kap * p_Am/ p_M;  % cm, ultimate structural length
E_m  = p_Am/ v ;         % J/ cm^3 [E_m], maximum reserve density
M_Em = E_m/ mu_E;        % mol/ cm^3 [M_Em], maximum reserve desity in c-moles
M_V  = d_V/ w_V;         % mol/ cm^3 [M_V], specific structural mass
% -----------------------------------------------------------------------

% --------------- unpack LEHRU ------------------------------------------
L   =  LEHRU(1); % cm, volumetric structural length
E   =  LEHRU(2); % J,   energy in reserve 
EH  =  LEHRU(3); % J, E_H maturity
ER  =  LEHRU(4); % J, E_R reproduction buffer
EB  =  LEHRU(5); % J, E_B energy in batch material
MQ  =  LEHRU(6); % nmolU, M_Q in body

% -------------- scaled state variables ----------------------------------
l     = L/ L_m;                           % -, scaled length
e     = E/ L^3/ E_m; eR = (ER + EB)/ L^3/ E_m;  % -, scaled energy in reserve and reproduction buffer

% ------------------- stress function and parameter modification --------
        
%s   = (MQ/ L^3 - NEC)* (MQ>NEC) / M_QT; % parameter modification is linear in the internal concentration

% for debugging purposes: sometimes intitial [M_Q] is egg is very high
% (structure close to zero) which makes for very important stress function
s = 0;

switch param
    
    case 'control'
E_G = E_G * 1;       
    
    case 'E_G'
E_G = E_G * (1 + s) * (E_G* (1 + s) >= 0);   % J/cm^3, [E_G] cost for synthesis of a unit structure
if E_G == 0
    disp('E_G is zero')
end
    case 'p_M'
p_M = p_M * (1 + s) ;   % J/d/cm^3, [p_M] volume specific somatic maintenance 

    case 'p_Am'
p_Am = p_Am * (1 - s) * (p_Am * (1 - s) >= 0); % J/d/cm^2 {p_Am}, surface area-specific assimilation
if p_Am == 0
    disp('p_Am is zero')
end

end
% -----------------------------------------------------------------------

if EH <= EHb % embryo stage
ML = 1;
pA   = 0;    
r = (E * v * ML/ L - p_M * L^3/ kap)/ (E + E_G * L^3/ kap) * (E >= p_M * L^4/ (kap * v * ML)) + ...
        (E * v * ML/ L - p_M * L^3/ kap)/ (E + kapG * E_G * L^3/ kap) * (E < p_M * L^4/ (kap * v * ML));pC   = E * (v/ L - r);                               % J d^-1, mobilization flux
L_b  = L;  
dER0  = 0;
dEB  = 0;


elseif EH <= EHj % uptil metamorphosis
ML = L/ L_b; 
pA   = f * p_Am * ML * L^2;    
r = (E * v * ML/ L - p_M * L^3/ kap)/ (E + E_G * L^3/ kap) * (E >= p_M * L^4/ (kap * v * ML)) + ...
        (E * v * ML/ L - p_M * L^3/ kap)/ (E + kapG * E_G * L^3/ kap) * (E < p_M * L^4/ (kap * v * ML));pC   = E * (v * ML/ L - r);                               % J d^-1, mobilization flux
L_j  = L;
dER0  = 0;
dEB   = 0;

elseif EH < EHp % juvenile II stage
ML = L_j/ L_b; 
pA   = f * p_Am * ML * L^2;   
r = (E * v * ML/ L - p_M * L^3/ kap)/ (E + E_G * L^3/ kap) * (E >= p_M * L^4/ (kap * v * ML)) + ...
        (E * v * ML/ L - p_M * L^3/ kap)/ (E + kapG * E_G * L^3/ kap) * (E < p_M * L^4/ (kap * v * ML));pC   = E * (v * ML/ L - r);                               % J d^-1, mobilization fluxb
dER0  = 0;
dEB  = 0;

else % adult: allocation to reproduction starts
    
ML   = L_j/ L_b;

    if male == 1 % batch preparation starts under stimulus of partner (assumed to be male here)
pCm  = E_m * (E_G * v * ML * L^2 + p_M * L^3)/ (kap * E_m + E_G);
dEB = kapR * ((1 - kap) * pCm - kJ * EHp)  * (ER0 > 0);
    else
        dEB = 0;
    end
    
if E >= p_M * L^4/ (kap * v * ML) % enough energy in reserve to cover somatic maintenance

    r    = (E * v * ML/ L^4 - p_M/ kap)/ (E/ L^3 + E_G/ kap); % d^-1, specific growth rate  
    pC   = E * (v * ML/ L - r);
    

dER0  = max(0, (1 - kap) * pC - kJ * EHp) - dEB; 

else %if E < p_M * L^4/ (kap * v * ML) % not enough energy in reserve to cover somatic maintenance
    
    if ER0 > 0% energy is taken from reproduction buffer to cover somatic maintenance
    r = 0;
    else       % nothing is left in the reproducion buffer and shrinking occurs
    r    = (E * v * ML/ L - p_M * L^3/ kap)/ (E + kapG * E_G  * L^3/ kap); % d^-1, specific growth rate  
    end
   
pC    = E * (v * ML/ L - r);
 
dER0  = max(0, (1 - kap) * pC - kJ * EHp) - (p_M * L^3 - kap * pC) * (ER0 > 0) * (ER0 > 0)  - dEB * (ER0 > 0 ) ;

end 

end
 
dEH  = max(0, (1 - kap) * pC - kJ * EH) * (EH < EHp);
dE    = pA - pC;                                          % J, energy in reserve  
dL    = r/ 3 * L;                                         % cm, structural length
P_WV  = 1 + M_Em/ M_V * P_EV * (e + eR);                  % mol, mass of uranium in entire body
dMQ   = k_e/l * P_Vd * c_d - MQ * k_e/ l/ P_WV ;          % nmol, M_Q, concentration of uranium 
  
% pack dLEHRU
dLEHRU = [dL; dE; dEH; dER0; dEB; dMQ];