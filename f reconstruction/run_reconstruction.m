%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trajectory reconstruction for BPA0
%-------------------------------------------------------------------------

clear all; close all; clc

global pets tf
pets = {'reconstruction'};

% 1-n vector with the knot-abscissa for the functional response trajectory
% that we are trying to reconstruct on the basis of growth data

tf = round(linspace(64,600,6)'); % choose the interval between knots which are estimated


%   'pars_init_method':
%     0 - get initial estimates from automatized computation (default)
%     1 - read initial estimates from .mat file (for continuation)
%     2 - read initial estimates from pars_init file
%   'results_output':
%     0 - prints results to screen (default)
%     1 - prints results to screen, saves to .mat file
%     2 - saves data to .mat file and graphs to .png files
estim_options('default'); 
estim_options('max_step_number',5e2); % change to 5e3 for the last round
estim_options('max_fun_evals',5e3);   % don't change this

estim_options('pars_init_method', 2); 
estim_options('results_output', 1); % 1, saves parameters into results_reconstruction
estim_options('method', 'nm'); % set no if you don't want to estimate

% % KEEP THESE COMMENTED IF YOU ARE NOT DOING CONTINUATIONS FROM THE
% RESULT.MAT FILE
% estim_pars; close all;
% estim_pars; close all;
% estim_pars; close all;
% estim_pars; close all;
% estim_pars; close all;
estim_pars; 

%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make additional plots here :
%-------------------------------------------------------------------------

clear par ;
%%%%%%% Start with parameters from the results
load('results_reconstruction.mat');
[data, auxData, metaData, txtData, weights] = mydata_reconstruction

%%%%%%% New predictions

cPar = parscomp_st(par); vars_pull(par); 
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
% compute temperature correction factors
TC = tempcorr(C2K(8.5), T_ref, T_A);

% initial conditions -
pars_UE0          = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
[U_E0, Lb, info]  = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve, f= 1
E0            = U_E0 * p_Am; % J, initial energy in the egg
LEH0  = [1e-5; E0; 0]; % 3-1 vector [J, cm, J] with initial conditions

% get tyf back
yf = zeros(length(tf),1); % create empty vector with knot-coordinates for the f values
for i = 1:length(tf)
eval(['yf(i) = f_',num2str(tf(i)),';'])
end
tyf = [tf yf];

% predict 
time = data.tW_gw150A(:,1);
[~, LEH] =  ode23s(@dget_LEH_for_reconstr, time, [LEH0; 0; 0; 0],[],tyf, TC, par, cPar); 
EW = LEH(:,1).^3 + w_E/ mu_E * LEH(:,2)/ d_E; % g, wet weight

diff= (EW-data.tW_gw150A(:,2))./data.tW_gw150A(:,2)*100;

figure()
plot(time, EW, 'g', 'linewidth',2); hold on
plot(time, data.tW_gw150A(:,2),'bs','markersize',10,'markerfacecolor','b');

figure()
plot(time, diff, 'g', 'linewidth',2);
ylabel('Diff');


%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the reconstructed trajectory:
%-------------------------------------------------------------------------

% prepare 2-n matrix with the first colum the knot-abscissa and the second
% column the knot-coordinates. The knot coordinates are the parameters
% which were estimated. 
% Use the results which were saved in results_reconstruction.mat
load('results_reconstruction.mat') ; 
vars_pull(par); 
tyfs = [tf, spline(tf, tyf)]; % smoothed trajectory - this version gives problems
figure()
plot(tyf(:,1), tyf(:,2),'bs','markersize',10,'markerfacecolor','b'); hold on
plot(tyfs(:,1), tyfs(:,2),'g', 'linewidth',2);
xlabel('time, day'); ylabel('scaled func resp., -')
legend('estimated knots','smoothed trajectory')
ylim([0 1]); 
set(gca,'fontsize',14)



% figure(6)
% 
% [prdData, info] = predict_reconstruction(par, data, auxData);
% plot(data.timescale, prdData.timescale,'b')

