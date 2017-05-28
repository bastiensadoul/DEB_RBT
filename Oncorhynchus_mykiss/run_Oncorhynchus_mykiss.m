close all; 

global pets

pets = {'Oncorhynchus_mykiss'};
check_my_pet(pets)
estim_options('default'); % runs estimation, uses nmregr method and filter
estim_options('max_step_number',5e3);  % set options for parameter estimation
estim_options('max_fun_evals',5e4);    % set options for parameter estimation

estim_options('pars_init_method',2);  % 1 from mat files (continuation)
estim_options('results_output', 0);    %
estim_options('method', 'no');
estim_pars;          % run estimation

