%% predict_tW
% starrlight sugustine & bastien sadoul
% last modified 2017/05/04
% computes growth as function of age and allows parameter to be modified
% according to physiological mode of action pMoA -
% the pMoA is a string that is a field in the auxData structure 

function [prdData, info] = predict_tW(par, data, auxData)                    

global study

cPar  = parscomp_st(par); 
vars_pull(par); vars_pull(cPar); vars_pull(data); vars_pull(auxData);

TC = tempcorr(temp.tW, T_ref, T_A); % -, TC temperature correction

% Creates a pMod vector if doesnt exist
if exist('pMod','var')== 0
    pMod = [];
end



% initial conditions -
% if study==150
%     e=e150;
% else
%     e=e124;
% end
pars_UE0          = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
[U_E0, Lb, info]  = initial_scaled_reserve(par.(study), pars_UE0); % d.cm^2, initial scaled reserve, f= 1
E0            = U_E0 * p_Am; % J, initial energy in the egg
LEH0  = [1e-5; E0; 0]; % 3-1 vector [J, cm, J] with initial conditions

age = data.tW(:,1);

% Creates vector of time to predict
if age(1) > 0
   ageIn = linspace(0,age(end),age(end)+1)';
   ageIn = union(ageIn, age);
 %   ageIn = [0; age];
else
   ageIn = linspace(age(1),age(end),age(end))';
   ageIn = union(ageIn, age);
%     ageIn = age;
end

% Predict
[~, LEH] =  ode23s(@dget_LEH, ageIn, [LEH0; 0; 0; 0],[],f, TC, par, cPar, pMod, pMoA); 
% if age(1) > 0
%     LEH(1,:) = [];
% end

% the three zeros are for N, Lb, and Lj (see inside the ODE subfunction)
% indice = ind2sub(ageIn,age);
[tf, loc] = ismember(age, ageIn);  
L = LEH(loc,1);   % cm, structural length  
E = LEH(loc,2);   % J, energy in reserve
EW = L.^3 + w_E/ mu_E * E/ d_E; % g, wet weight

prdData.tW = EW;
