function [L_out, Wd_out, MQ_out, LEQ ] = growth(f, p, p_tox, CI, c_d, param, time)
% function [L_out, Wd_out, MQ_out ] = growth(f, p, p_tox, CI, c_d, param, time)
% f: scaled functional response, c_d: mol/l environmental concentration
% p     = [v * TC; k_J * TC;  kap; E_Hb; E_Hj; E_Hp; p_Am * TC ; E_G * TC; p_M * TC; mu_E; d_V; w_V; w_E; del_M; del_Y; E_ERs; kap_R];
% p_tox = [P_EV; P_Vd; k_e; M_QT; NEC] 
% Lbj = [L_b; L_j] if starting after a=0  or Lbj = [] if starting with a=0;
% calls ode_LEHRU
% outputs:
%- L_out:  mm observed physical length, 
%- Wdout:  µg dry mass, 
%- MQ_out: nmol/ g observed internal concentration
%- LEQ n-3 matrix with state variables: L(cm), E (J) and [M_Q] (nmol/l)

    % no batch preparation
male = 0;

 mu_E = p(10); d_V = p(11); w_E = p(13); del_Y = p(15); 

     % -, shape coefficient to convert vol-length to physical length

[~, LEHRU] = ode45(@ode_LEHRU_bj,time,CI,[],p,p_tox,f,c_d, param, male);
L  = LEHRU(:,1);                        % cm, structural length
E  = LEHRU(:,2);                        % J, energy in reserve and reproduction buffer
MQ = LEHRU(:,6);                        % mol/cm^3, amount of uranium per unit structure
Q  = MQ .* L.^3;                        % mol, mass of uranium in entire body

LEQ = [L E MQ]; % pack state variables

Wd = d_V * L.^3 + w_E/ mu_E * E;        % g, dry weight 
L_out   = L * 10/ del_Y;                % mm, observed physical length
Wd_out  = Wd * 1e6;                     % µg, dry mass
MQ_out  = Q./ Wd ;                       % nmol/ g, observed internal concentration 

