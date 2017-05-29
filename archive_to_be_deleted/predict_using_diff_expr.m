
clear all; close all; clc
load('results_Oncorhynchus_mykiss.mat') % load parameter values of the blank taken from AmP
parcomp=parscomp_st(par);

% f and T
f=1;
T=C2K(8.5);
TC= exp(par.T_A/ par.T_ref - par.T_A/ T);

% Maximum structural length is z from pars_init
Lm = par.z ;           % cm

% calculate p_Am
p_Am = Lm * par.p_M/ par.kap ;       % J/d/cm2

% Initial V, Lb and Lj
V0=0.001;
Lbini = 0 ;                         % Lb is used only after Eh>Ehb in ode
Ljini = 0 ;                         % Lj is used only after Eh>Ehj in ode

% Calculate E0
f_mom = 1 ;                          % we assume that the mother had f=1
pars_UE0          = [parcomp.V_Hb; parcomp.g; par.k_J; parcomp.k_M; par.v]; % compose parameter vector
[U_E0, Lb, info]  = initial_scaled_reserve(f_mom, pars_UE0); % d.cm^2, initial scaled reserve, f= 1
E0                = U_E0 * p_Am; % J, initial energy in the egg

EVEhErLbLj0=[E0 ; V0 ; 0 ; 0 ; Lbini ; Ljini ];

% Run ode on given time period
timeframe=linspace(0, 900, 1000);
calcpar.p_Am=p_Am;

[~, E_V_Eh_Er_Lb_Lj] = ode23s(@dget_dE_dV_dEh_dEr_dLb_dLj, timeframe, EVEhErLbLj0, [], f, par, calcpar, TC);

E=E_V_Eh_Er_Lb_Lj(:,1);
V=E_V_Eh_Er_Lb_Lj(:,2);
Eh=E_V_Eh_Er_Lb_Lj(:,3);
Er=E_V_Eh_Er_Lb_Lj(:,4);
Lb=E_V_Eh_Er_Lb_Lj(:,5);
Lj=E_V_Eh_Er_Lb_Lj(:,6);

L=V.^(1/3)/par.del_M;



% Get length from L and Lm.
l=L*Lm ;

% Get dry weight Structure + Reserve + Reproduction
Wd = par.d_V * V + (parcomp.w_E/par.mu_E) * E /par.d_E + (parcomp.w_E/par.mu_E) * Er / par.d_E;   % w_E and mu_E are in add_chem. w_E/mu_E is in g/J

% Get wet weight
Ww = Wd / par.d_V  ;  %%%%%%%?????????????????????????????????????????????????

% Plot

figure(1)
plot(timeframe, V, 'b','linewidth',3);
ylabel('Structure, cm^3');

figure(2)
plot(timeframe, l, 'b','linewidth',3);
ylabel('length, cm');

figure(3)
plot(timeframe, E, 'b','linewidth',3);
ylabel('reserve, J');

figure(4)
plot(timeframe, Er, 'b','linewidth',3);
ylabel('Er, J');

figure(5)
plot(timeframe, Eh, 'b','linewidth',3);
ylabel('Eh, J');

figure(6)
plot(timeframe, Ww, 'b','linewidth',3);
ylabel('Wet weight, g');

%%%%%%%%% Add data points
% tW_gw150=[...
%     0.0        0.1301370      1.0000000        0.1280822      1.0000000        0.1264407      1.0000000
%    14.0        0.2308772      0.9760274        0.2238596      0.9760274        0.2456349      0.8542373
%    28.0        0.4258993      0.9520548        0.4257143      0.9589041        0.4395062      0.8237288
%    42.0        0.7931408      0.9486301        0.7700000      0.9589041        0.7590909      0.8203390
%    56.0        1.3334545      0.9417808        1.2082143      0.9589041        1.2599174      0.8203390
%    70.0        2.1609489      0.9383562        1.9246429      0.9589041        2.0731405      0.8203390
%    84.0        2.9242424      0.9349315        2.6357843      0.9554795        2.7353293      0.8169492
%    98.0        4.1257576      0.9349315        3.8529412      0.9554795        4.0053892      0.8169492
%   112.0        6.2828283      0.9349315        5.5490196      0.9554795        5.7710843      0.8120572
%  126.0        8.0085859      0.9349315        7.7745098      0.9554795        7.8078313      0.8120572
%  140.0       10.4131980      0.9302096       10.2745098      0.9554795        9.7969880      0.8120572
%  154.0       14.6192893      0.9302096       13.6764706      0.9554795       13.3734940      0.8120572
%  175.0       22.0558376      0.9302096       21.8137255      0.9554795       21.3253012      0.8120572
%  196.0       31.8622449      0.9254878       31.8137255      0.9554795       31.4457831      0.8120572
%  217.0       42.5765306      0.9254878       42.0833333      0.9554795       41.9578313      0.8120572
%  245.0       59.9744898      0.9254878       58.0392157      0.9554795       59.5481928      0.8120572
%  273.0       76.7692308      0.9207659       76.5686275      0.9554795       76.7575758      0.8120572
%  357.0       114.1857143      0.9113221      120.3918033      0.9414283      116.1598639      0.8120572
% %  357.5      119.2000000      0.9113221      124.8000000      0.9414283     118.0000000      0.8120572                 %  because of cull effect
% %  385.0      163.6111111      0.9113221      163.5000000      0.9414283      167.5675676      0.7958161
% %  412.0      187.9166667      0.9113221      189.1250000      0.9414283      191.7567568      0.7958161
% ];
% 
% tW_gw150(:,1)=tW_gw150(:,1)+64;         % to put in dpf
% data.tW_gw150A = tW_gw150(:,[1 2]);



%%%%%%%%% DIFFERENTIAL EXPRESSIONS for Reserve (E), structural volume,
%%%%%%%%% maturity and reproduction

function dE_dV_dEh_dEr_dLb_dLj = dget_dE_dV_dEh_dEr_dLb_dLj(t, EVEhErdLbdLj, f, par, calcpar, TC)

    %%%%% Initial variables
    E=EVEhErdLbdLj(1);
    V=EVEhErdLbdLj(2);
    Eh=EVEhErdLbdLj(3);
    Er=EVEhErdLbdLj(4);
    Lb=EVEhErdLbdLj(5);
    Lj=EVEhErdLbdLj(6);


    %%%%% To take into acount the V1-morph phase    ??????????????????????????
    if Eh < par.E_Hb
      s_M = 1;         % It is Lb/Lb
    elseif Eh < par.E_Hj && Eh >= par.E_Hb
      s_M = V^(1/3)/Lb;                 
    else
      s_M = Lj/Lb;       
    end


    %%%%% Correction by temperature
    p_TAm = calcpar.p_Am * TC * s_M;
    p_TM = par.p_M * TC * s_M;
    vT = par.v * TC * s_M;


    %%%%% What is entering: nothing before birth
    if (Eh >= par.E_Hb)
        pA = f*p_TAm*V^(2/3);
    else
        pA=0 ;
    end
    print(pA);

    %%%%% What is going to somatic maintenance
    pS = p_TM*V;
    %pS = p_TM*V + p_TT*V^(2/3);   % from the book ???  don't understand!

    %%%%% What is going out of reserve : [p_C], J/d (2.12, Kooijman 2010)
    pC = E*(par.E_G*vT*V^(2/3)+pS)/(par.kap*E+par.E_G*V);      

    %%%%% What is going towards structure
    pG = par.kap*pC - pS;

    %%%%% Increase of Volume
    dV = pG / par.E_G;    % or called r 

    %%%%% What is going to maturity maintenance
    pJ = par.k_J * Eh;

    %%%%% What is going to maturity / repro
    pR = (1-par.kap) * pC - pJ;

    %%%%% Increase of maturity till E_Hp after in reproduction
    if Eh < par.E_Hp
        dEh = pR;
        dEr=0;
    else 
        dEh=0;
        dEr = pR;
    end

    %%%%% Increase of Lb and Lj
    dL=dV.^(1/3)/par.del_M;

    if Eh < par.E_Hb
      dLb = dL ;
    else
      dLb=0;
    end
    
    if Eh < par.E_Hj
      dLj = dL;                 
    else
      dLj = 0;
    end

    %%%%% dE
    dE=pA-pC ;   % J/d

    %%%%%%%%%%%%%%% Pack output 
    dE_dV_dEh_dEr_dLb_dLj = [dE; dV; dEh; dEr; dLb; dLj]; 

end





