% ---------------- AuguGagn2011 b ----------------------------------------
% estimates toxico-kinetic parameters for zebrafish
% accumulaation in embryo + juvenile: BourSimo08,
% depuratoin from egg and early juv: Bou09
% accumulation in adult: BariBuet2005, BariAdam2007
% created: beginning october 2011
% last modified: october 21 2011

%clear all; 
clc; close all


% define parameter which is modified: 'p_M', 'E_G', or 'p_Am'
global param
% param = 'p_Am';
param = 'E_G'; % 'E_G', 'p_Am', 'p_M'


          
% ----- TIME VECTOR for embryo/juv observations (25+273 K):
t = (0:64)'; % age vector

          

  
  % --------- PARAMETER VECTOR:
  f_SBemb = 0.8;
  
 % ------- tox parameters specific to each mode of action:
switch param
    
    case 'E_G'

NEC      = 0;%0.003258 * 0.5;    % nmol U cm^-3 structure, no effect internal concentration 
M_QT     = 0.2109;            %*3000;%*0.8;  % nmol U cm^-3 structure, tolerance concentration
P_Vd = 0.00389*1.3;%*0.4; 
P_EV     = 0.222*0.5*3*100;         %*0.5;%1.038;  % -, partition coefficient reserve/structure
k_e      = 0.0006 * 10 * 0.5;%4e-005;  % d^-1, elimination rate


    case 'p_M'

NEC  = 0.003258;    % nmol U cm^-3 structure, no effect internal concentration 
M_QT = 0.2109 * 600;%*0.8;  % nmol U cm^-3 structure, tolerance concentration
P_Vd =  0.0191*800;     % 0.07004*5;  % liter water cm^-3 V, partition coefficient water/structure
P_EV = 0.222*10;%1.038;  % -, partition coefficient reserve/structure
k_e  = 0.0006;%4e-005;  % d^-1, elimination rate
P_Vd_emb = 0.00389*10; 

    case 'p_Am'
        
NEC  = 0.03591;       % nmol U cm^-3 structure, no effect internal concentration 
M_QT = 0.6938;        % nmol U cm^-3 structure, tolerance concentration
P_Vd = 0.0009282;      % liter water cm^-3 V, partition coefficient water/structure
P_EV = 0.04316;       % -, partition coefficient reserve/structure
k_e  = 0.021;         % d^-1, elimination rate

end

par_txt = {...
    'NEC, nmol/cm^3, no effect internal concentration'
    'M_QT, nmol/ cm^3, tolerance concentration'
    'P_Vd, liter water cm^-3 V, partition coefficient water/structure'
    'P_EV, -, partition coefficient reserve/structure'
    'k_e, 1/d, elimination rate'
    ' f_SBemb -, scaled functional response, BourSimo2008'
    };

par = [...
    NEC     0        % -1- nmol U cm^-3 structure, no effect internal concentration
    M_QT    0        % -2- nmol U cm^-3 structure, tolerance concentration
    P_Vd    1        % -3- liter water cm^-3 V, partition coefficient water/structure
    P_EV    0        % -4- -, partition coefficient reserve/structure
    k_e     0        % -5- d^-1, elimination rate
    f_SBemb 0        % -6- -, scaled functional response, BourSimo2009
    ];     
 

if 0
par_new = nmregr('estimation_AuguGagn2011b_EmbJuv',par, BourSimo08,  Bour09F1); 
else 
    par_new = par(:,1);
end

[eBourSimo08, eBour09F1 ] = estimation_AuguGagn2011b_EmbJuv(par_new, BourSimo08,  Bour09F1 );

fprintf('Mode of action %s .\n', param)

printpar(par_txt,par_new,par(:,1),'predicted par & old par') 

% ------------------ PLOT RESULTS ----------------------------------------
figure(1)
subplot(221)
plot(t([3 10]), BourSimo08(10:11,2), 'm*', t([3 10]), eBourSimo08(10:11,1), 'mo')
hold on
plot(t([3 10]), BourSimo08(12:13,2), 'r*', t([3 10]), eBourSimo08(12:13,1), 'ro')
box off
xlabel('age, dpf')
ylabel('nmol U/g W_d')
xlim([0 9])

% -------------------------------------------
subplot(222) 
plot(t([3 10]), Bour09F1(1,1:2),'b*',  t(10), eBour09F1(1), 'bo',...
     t([3 10]), Bour09F1(2,1:2),'m*',  t(10), eBour09F1(2),'mo', ...
     t([3 10]), Bour09F1(3,1:2),'r*',  t(10), eBour09F1(3),'ro')
xlabel('age, dpf')
ylabel('nmol U/ g')
% -------------------------------------------
subplot(223)
plot(t([3 10]), BourSimo08(4:5,2), 'b*', t([3 10]), eBourSimo08(4:5,1), 'bo',...
     t([3 10]), BourSimo08(6:7,2),'m*',  t([3 10]), eBourSimo08(6:7,1),'mo', ...
     t([3 10]), BourSimo08(8:9,2),'r*',  t([3 10]), eBourSimo08(8:9,1),'ro')
xlabel('age, dpf')
ylabel('dry mass, \mu g')
xlim([0 9])
box off

% -------------------------------------------
subplot(224)
plot(t(10), BourSimo08(1,2), 'b*', t(10), eBourSimo08(1,1), 'bo', ...
     t(10), BourSimo08(2,2),'m*',  t(10), eBourSimo08(2,1),'mo',...
     t(10), BourSimo08(3,2),'r*',  t(10), eBourSimo08(3,1),'ro')
xlabel('age, dpf')
ylabel('total length, mm')
xlim([0 9])
box off










