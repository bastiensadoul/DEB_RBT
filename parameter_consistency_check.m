clear all; clc

load('results_Oncorhynchus_mykiss.mat');

% food is 42% protein, 16% fat
% fat  is roughly 17.2 kJ/g , protein is roughly 38.9 kJ/g for protein
% energy content of food in kJ/g can be computed like this:
% JX = 17.2 * 0.16 + 38.9 * 0.42; % kJ/g

JX = 21.9 ; % kJ/d, pers. comm. Steve Sommerfelt 18/02/16

% we don't know the digestion efficiency 
kap_X = par.kap_X;

% use this to get p_Am needed lower down
c = parscomp_st(par);


% the data:
tLW = [...
30	NaN	1	12  NaN ;
61	NaN	2	13.9	0.061612903 ;
91	NaN	12	14	0.53 ;
122	NaN	37	13.5	0.588709677 ;
152	NaN	68	13.4	1.105666667 ;
183	20.14568672	130	13.3	2.58 ;
213	21.5443469	182	12.5	1.82 ;
244	25.04752863	297	12.7	4.340322581 ;
274	30.37304642	566	13.5	12.374 ;
305	32.15692038	685	13.6	6.602580645 ;
335	35.37616057	943	14.8	13.244 ;
366	38.472504	1230	15.2	15.18322581 ;
396	41.5666287	1501	13	17.43433333 ;
427	43.30037284	1713	13	12.99354839 ;
457	45.39591421	2002	13.1	11.46366667 ;
488	47.81766661	2307	13	19.08709677 ;
518	49.18429576	2570	12.7	16.13066667 ;
549	51.38700584	2836	12.6	10.98322581 ;
579	53.4819702	3136	12.5	11.1 ;
610	55.32407387	3556	12.5	11.92258065 ;
640	57.17912711	4038	12.5	18.47666667 ;
671	59.42611275	4533	12.7	19.16129032 ;
701	60.72011595	4858	12.9	24.15833333 ;
732	61.56382501	4970	13.4	22.32774194 ;
762	62.32757629	5012	13.6	30.24 ;
793	61.98695269	5097	13.4	1.946774194];

% correct parameter for experimental temperature :
Texp = mean(C2K(tLW(:, 4)));
TC = tempcorr(Texp, par.T_ref, par.T_A);
% calculate the structural length
L = tLW(6:end,2) .* par.del_M; % cm, structural length
% Calculate the food ingested in Joules per day assuming maximum feeding
% (f=1) :
pX_pred = c.p_Am * TC .* L.^2/ kap_X; % J, predicted maximum food intake for each size

% "Real" energetic value of the food
Mx = tLW(2:end,5); % g of food
pX = JX .* Mx * 1e3; % J in the food

plot(tLW(2:end,1), pX, 'ro');
hold on
plot(tLW(6:end,1), pX_pred, 'bo');
xlabel('time since hatch')
ylabel('J given')
legend('observed','predicted')





