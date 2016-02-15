close all
global pets

pets = {'Oncorhynchus_mykiss'};
check_my_pet(pets); 

estim_options('default'); % runs estimation, uses nmregr method and filter
estim_options('max_step_number',5e3);  % set options for parameter estimation
estim_options('max_fun_evals',5e3);    % set options for parameter estimation

% 'pars_init_method': 0 - get initial estimates from automatized computation (default)
%                     1 - read initial estimates from .mat file 
%                     2 - read initial estimates from pars_init file 
% 'results_output':   0 - prints results to screen; (default)
%                     1 - prints results to screen, saves to .mat file
%                     2 - saves data to .mat file and graphs to .png files
%                     (prints results to screen using a customized results file when there is one)
% 'method':           'nm' - use Nelder-Mead method (default); 'no' - do not estimate;

estim_options('pars_init_method', 2);
estim_options('results_output',1);
estim_options('method', 'no');

estim_pars;