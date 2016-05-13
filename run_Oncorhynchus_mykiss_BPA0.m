clear all; close all; clc

global pets

% pets = {'Oncorhynchus_mykiss'};
% pets = {'Oncorhynchus_mykiss_BPA0'};
 pets = {'Oncorhynchus_mykiss_BPA03to300'};


estim_options('default'); % runs estimation, uses nmregr method and filter

estim_options('max_step_number',500);  % set options for parameter estimation
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
estim_options('results_output', 2);
estim_options('method', 'nm');
% estim_options('filter', 0); % we no longer have all of the model parameters from the 'abj' model, so it is best to put this to zero

estim_pars;          % run estimation