function [data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss

%% set metadata
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Salmoniformes'; 
metaData.family     = 'Salmonidae';
metaData.species    = 'Oncorhynchus_mykiss'; 
metaData.species_en = 'Rainbow trout'; 
metaData.T_typical  = C2K(15.5); % K, body temp
metaData.data_0     = {'ah_T'; 'ab_T'; 'ap'; 'am'; 'Lb'; 'Lp'; 'Li'; 'Wd0'; 'Wdh'; 'Wdb'; 'Wwi'; 'Ri'};  % tags for different types of zero-variate data
metaData.data_1     = {'L-Ww'; 't-Ww'; 't-L'; 'tWde'; 'tWde_E'; 'T-ah'; 'Ww-JO'; 'WLO'; 'Wie1985'; 'tW150'; 'tW124ini'; 'tW124fin'}; % tags for different types of uni-variate data

metaData.COMPLETE = 2.4; % using criteria of LikaKear2011

metaData.author   = {'Bas Kooijman'};        
metaData.date_subm = [2014 09 26];                           
metaData.email    = {'bas.kooijman@vu.nl'};                 
metaData.address  = {'VU University Amsterdam'}; 

metaData.author_mod_1   = {'Starrlight Augustine'};        
metaData.date_mod_1 = [2016 01 26];                           
metaData.email_mod_1    = {'starrlight.augustine@akvaplan.niva.no'};                 
metaData.address_mod_1  = {'akvaplan-niva'};

metaData.author_mod_2   = {'Bastien Sadoul';'Starrlight Augustine' };        
metaData.date_mod_2 = [2016 02 01];                           
metaData.email_mod_2    = {'bastien.sadoul@hotmail.fr';'starrlight.augustine@akvaplan.niva.no'};                 
metaData.address_mod_2  = {'University of Calgary';'akvaplan-niva'};

metaData.author_mod_3   = {'Bas Kooijman'};        
metaData.date_mod_3 = [2016 04 07];                           
metaData.email_mod_3    = {'bas.kooijman@vu.nl'};                 
metaData.address_mod_3  = {'VU University Amsterdam'};

metaData.author_mod_4   = {'Bastien Sadoul';'Starrlight Augustine' };       
metaData.date_mod_4 = [2017 04 15];                           
metaData.email_mod_4    = {'bastien.sadoul@hotmail.fr';'sta@akvaplan.niva.no'};                 
metaData.address_mod_4  = {'University of Calgary';'Akvaplan-niva'};
% we added respiration data and re-did parameter estimation

metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
metaData.date_acc    = [2017 04 16]; 

%% set data
% zero-variate data
    data.Wd0 = 0.1007 * 0.4111;  units.Wd0 = 'g'; label.Wd0 = 'dry egg weight'; bibkey.Wd0 = 'FromRasm1991';   
comment.Wd0 = 'wet weight multiplied by percent dry matter of large eggs, see Table 1'; 
    data.ah = 19 + 13;  units.ah = 'd'; label.ah = 'age at hatching at 10 degrees'; bibkey.ah = 'FromRasm1991';   
temp.ah = C2K(10); units.temp.ah = 'K'; label.temp.ah = 'temperature';
    data.ah_5 = 33 + 34;  units.ah_5 = 'd'; label.ah_5 = 'age at hatching at 5 degrees'; bibkey.ah_5 = 'FromRasm1991';   
temp.ah_5 = C2K(5); units.temp.ah_5 = 'K'; label.temp.ah_5 = 'temperature';
    data.Wdh = 0.1101 * 0.3523;  units.Wdh = 'g'; label.Wdh = 'weight at hatch'; bibkey.Wdh = 'FromRasm1991';   
comment.Wdh = 'large eggs, 10 deg C, Table 2, wet weight times percent dry matter';
    data.ab = 19 + 13 + 22; units.ab = 'd'; label.ab = 'age at birth at 10 degrees'; bibkey.ab = 'FromRasm1991';   
temp.ab = C2K(10); units.temp.ab = 'K'; label.temp.ab = 'temperature';
    data.ab_5 = 33 + 34 + 52; units.ab_5 = 'd'; label.ab_5 = 'age at birth at 5 degrees'; bibkey.ab_5 = 'FromRasm1991';   
temp.ab_5 = C2K(5); units.temp.ab_5 = 'K'; label.temp.ab_5 = 'temperature';
    data.Wdb = 0.1665 * 0.1929;  units.Wdb = 'g'; label.Wdb = 'weight at birth'; bibkey.Wdb = 'FromRasm1991';   
comment.Wdb = 'large eggs, 10 deg C, Table 2, wet weight times percent dry matter';
    data.ap =28 + (20 * 30); units.ap = 'd';    label.ap = 'age at puberty';         bibkey.ap = 'DaviKenn2014';
 temp.ap = C2K(13); units.temp.ap = 'K'; label.temp.ap = 'temperature';
 comment.ap = 'GSI assessment indicated that rapid egg growth started at 20 Mo post hatch, and we assume that it takes about 28 d to hatch according to Vels1987';
data.am = 11*365;  units.am = 'd';    label.am = 'life span';              bibkey.am = {'fishbase'};   
  temp.am = C2K(5); units.temp.am = 'K'; label.temp.am = 'temperature';
data.Lp = 54;      units.Lp = 'cm';   label.Lp = 'forked length at puberty'; bibkey.Lp = 'DaviKenn2014';
data.Wp = 3.5 * 1e3;   units.Wp = 'g';    label.Wp = 'wet weight at puberty';    bibkey.Wp = 'DaviKenn2014';
    comment.Wp = 'weigth at 20 montsh post hatch, when rapid egg growth occurred'; 
data.Li = 120;     units.Li = 'cm';   label.Li = 'ultimate total length';  bibkey.Li = 'fishbase';
data.Wi = 25400;   units.Wi = 'g';    label.Wi = 'ultimate wet weight';    bibkey.Wi = 'fishbase';
data.Ri = data.Wi * 2.5/ 365; units.Ri = '#/d'; label.Ri = 'maximum reprod rate'; bibkey.Ri = 'Wiki';   
  temp.Ri = C2K(5); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
  comment.Ri = '2000 till 3000 eggs per kg';
  
% uni-variate data

% L_Ww data
data.LW = [... % length (cm), wet weight (g)
5.859	2.310; 9.567	9.342; 12.167	14.040; 12.552	43.337
12.889	7.340; 13.178	29.876; 13.756	32.171; 14.334	38.968
16.597	61.652; 24.928	170.345; 25.217	111.830; 27.095	195.273
27.769	251.609; 27.817	213.339; 28.780	224.668; 28.925	179.651
29.888	260.774; 30.610	323.867; 30.658	449.950; 31.332	368.949
32.055	308.216; 32.681	317.268; 33.355	362.347; 35.136	484.057
37.207	504.475; 37.400	738.636; 39.326	612.701; 40.722	716.371
42.119	824.543; 42.648	705.258; 43.274	813.373; 43.467	887.684
43.660	1009.274; 45.345	1173.753; 46.116	1009.459; 48.138	1095.164; 49.920	1376.723];
units.LW = {'cm', 'g'};  label.LW = {'length', 'wet weight'};  bibkey.LW = {'ChenSnow2015', 'BudyThie2002',  'StraStut1997',  'WeatGill1981'};
comment.LW = 'Compiled from 4 publications';

% T-ab data from Velsen 1987
% given as the 50% value
Tah = [...
2	115; 2.5	106; 2.8	93; 3	111; 3.2	101;4.5	72.9; 4.8	75; 5	72
5	68; 5	64; 5.8	63; 6.1	61; 6.2	61; 6.5	57.5;7	56
7	60; 7.2	45; 7.5	43; 7.5	43; 7.7	48; 7.7	44; 7.7	46.5; 7.8	48; 7.8	49
7.8	44; 7.9	42; 7.9	43; 7.9	48; 7.9	46; 7.9	46; 8	41; 8.7	40.3
9.2	41; 9.3	35; 9.5	36; 9.5	36; 9.5	36; 9.5	36; 10	38; 10	34; 10	33; 10	35.5; 10.3	29.6
10.3	28; 10.4	32.1; 10.7	27; 10.7	29.2; 10.8	29.5; 11.3	30; 11.5	30.3; 11.5	27
11.5	28; 11.7	24; 12	25; 12.2	26; 12.2	23; 12.4	24; 12.5	27
12.8	24.5; 12.9	18; 13	28; 14.5	21; 15	26; 15	22; 15.5	18; 17.5	18];
Tah(:,1) = C2K(Tah(:,1));
data.Tah = Tah; units.Tah = {'K','d'}; label.Tah = {'Temperature', 'age at hatch'}; bibkey.Tah = {'Vels1987'};
comment.Tah = 'age is given as the 50% value';

% t-Wd-Wdyolk from NinnStev2006
% age dpf, dry W mg of embryo, dry W mg of yolk
tWdY=[...
24	1.01	32.74;
30	3.08	31.8;
45	12.6	15.55;
50	22.91	8.24;
60	46.55	2.53;
75	118.93	0;
90	230.55	0];
data.tWde = tWdY(:,[1 2]); % yolk-free embryo dry mass, mg 
units.tWde = {'d','mg'} ; label.tWde = {'time since fertilization', 'yolk-free dry weight'};
bibkey.tWde ={'NinnStev2006'};
comment.tWde ={'Table 1, pp 1878, chorionated'};
temp.tWde = C2K(10); units.temp.tWde = 'K'; label.temp.tWde = 'temperature'; 

data.tWde_E = tWdY(:,[1 3]); % yolk dry mass, mg 
units.tWde_E = {'d','mg'} ; label.tWde_E = {'time since fertilization', 'yolk dry weight'};
bibkey.tWde_E ={'NinnStev2006'};
comment.tWde_E ={'Table 1, pp 1878, chorionated'};
temp.tWde_E = C2K(10); units.temp.tWde_E = 'K'; label.temp.tWde_E = 'temperature'; 

% DaviKenn2014
% Colums of tLW:
%  1  t days post hatch
%  2  L cm, fork length 
%  3 W g, wet weight
%  4 T degC, not constant
tLW = [...
30	NaN 1	12;
61	NaN	2	13.9;
91	NaN	12	14;
122	NaN	37	13.5;
152	NaN	68	13.4;
183	20.14568672	130	13.3;
213	21.5443469	182	12.5;
244	25.04752863	297	12.7;
274	30.37304642	566	13.5;
305	32.15692038	685	13.6;
335	35.37616057	943	14.8;
366	38.472504	1230	15.2;
396	41.5666287	1501	13;
427	43.30037284	1713	13;
457	45.39591421	2002	13.1;
488	47.81766661	2307	13;
518	49.18429576	2570	12.7;
549	51.38700584	2836	12.6;
579	53.4819702	3136	12.5;
610	55.32407387	3556	12.5;
640	57.17912711	4038	12.5
671	59.42611275	4533	12.7;
701	60.72011595	4858	12.9;
732	61.56382501	4970	13.4;
762	62.32757629	5012	13.6;
793	61.98695269	5097	13.4];
tL  = tLW(6:end,[1 2]);
tW = tLW(:,[1 3]) ;

data.tL1  = tL(1:10,:); units.tL1 = {'d', 'cm'}; label.tL1  = {'days post hatch', 'fork length'}; bibkey.tL1 = {'DaviKenn2014'};
comment.tL1 = 'fish reared in water recirculating system, we use the mean temperature';
data.tL2  = tL(11:end,:); units.tL2 = {'d', 'cm'}; label.tL2  = {'days post hatch', 'fork length'}; bibkey.tL2 = {'DaviKenn2014'};
comment.tL1 = 'fish reared in water recirculating system, we use the mean temperature';

data.tW1 = tW(1:10,:) ;   units.tW1 = {'d', 'g'} ; label.tW1 = {'days post hatch', 'wet weight'};  bibkey.tW1 = {'DaviKenn2014'} ;
temp.tW1 = mean(C2K(tLW(:, 4))); units.temp.tW1 = 'K' ;  label.temp.tW1 = 'mean temperature in recirculation system' ; 
comment.tW1 = 'fish reared in water recirculating system, we use the mean temperature' ;

data.tW2 = tW(11:end,:) ;   units.tW2 = {'d', 'g'} ; label.tW2 = {'days post hatch', 'wet weight'};  bibkey.tW2 = {'DaviKenn2014'} ;
temp.tW2 = mean(C2K(tLW(:, 4))); units.temp.tW2 = 'K' ;  label.temp.tW2 = 'mean temperature in recirculation system' ; 
comment.tW2 = 'fish reared in water recirculating system, we use the mean temperature' ;


% t-Ww data from YaniHisa2002 at T = 273 + 8.5
% initial weight 1.54 g
data.tW3 = [... % time (d), wet weight (g)
0.202	1.471
15.609	2.434
31.397	3.448
46.426	4.716
62.217	6.136
77.253	8.215
92.859	10.294
108.271	11.917
123.688	14.148
139.104	16.227
154.725	20.335];
units.tW3 = {'d', 'g'};  label.tW3 = {'time', 'wet weight'};  bibkey.tW3 = {'YaniHisa2002'};
 temp.tW3 = C2K(8.5); units.temp.tW3 = 'K'; label.temp.tW3 = 'temperature';
auxData.W0.tW3 = 1.471;  units.W0.tW3 = 'g';  label.W0.tW3 = 'wet weight at t = 0'; 

%  1 days post fertization, 
%  2 W g, wet weight
%  3 -, maturation stage
%  4 -, actual date given in the article
% tWGomeWeill=[...
% 0	NaN	'Immature'	'OctNov_year1' ;
% 441	508	'Immature'	'Jan_year3' ;
% 517	NaN	'First_signs_of_gametogenesis'	'MarApr_year3' ;
% 745	1755	'Mature'	'Nov_year3'];

% KieAls1998
% Colums of WLO:
%  1 L cm, fork length
%  2 W g, wet weight
%  3 T C, temperature
%  4 Swimming speed (0, 45% of the Ucrit and 75% of the Ucrit
%  5 Oxygen consumption (umol/g/h)
WLO = [...
11.5    15.2	5	0	4.10689479
11.3    18      15	0	5.65205184
% 12.2    20.6    5	45	4.98310919
% 11.6    17.5    15	45	6.58636376
% 11.4    16.8	5   75	6.66874162
% 11.1    15.1    15	75	9.93660271
]; % the data with current are commented
[Y,I]=sort(WLO(:,1)); WLO=WLO(I,:); % sorts by increasing length
%[Y,I]=sort(WLO(:,2)); WLO=WLO(I,:); % sorts by increasing weight
WLO(:,5) = 24 .*1e-3 * WLO(:,5) .* WLO(:,2);  % umol/g/h to mmol/d
data.WJO = [WLO(:, 2) WLO(:, 5)];
units.WJO = {'g', 'mmol/d'}; label.WJO = {'wet weight'; 'oxygen consumption'}; bibkey.WJO = {'KieAls1998'};
 comment.WJO = 'We only use the no current data. The author state said fed to satiation but organisms seem a bit light.';
temp.WJO = C2K(WLO(:,3)); units.temp.WJO = 'K' ;  label.temp.WJO = 'temperature' ; 
forkLength.WJO = WLO(:,1); units.forkLength.WJO = 'cm' ;  label.forkLength.WJO = 'fork length' ; 

% Wieser 1985
% Colums of WLO:
%  1 T C, temperature
%  2 W g, wet weight
%  3 O umol/h,  Oxygen consumption
Wie1985=[...
4	0.07411815		0.29960388
4	0.07570116		0.24584609
4	0.09082459		0.30265766
4	0.09228082		0.25088331
4	0.28830376		1.28400505
4	0.28968824		1.72302612
4	0.40204423		2.66487866
4	0.83844506		2.59816420
4	0.85591715		3.04041960
4	0.99068562		2.50756105
4	1.01043342		2.12238970
4	1.60382871		2.36219106
4	2.66622647		7.11478856
4	2.7348432		8.56534279
4	4.18332479		11.76698089
4	4.43520878		13.71831449
4	4.53190409		18.64541439
4	7.40732736		21.20118267
4	7.4091086		17.92774796
12	0.07149162		0.58938284
12	0.07288104		0.32865209
12	0.09710472		1.00120501
12	0.15387105		1.08926586
12	0.15572133		1.36932418
12	0.18599371		1.80633860
12	0.23481131		2.11247119
12	0.23486762		1.78471567
12	0.64984618		4.79114518
12	0.9625608		4.96740657
12	0.98682057		3.92773573
12	1.15678331		6.47429477
12	1.16328923		9.40429751
12	3.11377238		13.96525595
12	3.14317395		18.82591149
12	4.82204969		37.07129187
12	4.82326097		31.06914016
12	5.94481937		43.66841489
12	7.16758909		32.60241238
12	7.40366066		41.34829138];

Wie1985(:,3) = 24 .*1e-3 * Wie1985(:,3);  % umol/h to mmol/d
[Y,I]=sort(Wie1985(:,2)); Wie1985=Wie1985(I,:); % sorts by increasing weight
data.Wie1985 = [Wie1985(:, 2) Wie1985(:, 3)];
units.Wie1985 = {'g', 'mmol/d'}; label.Wie1985 = {'wet weight'; 'oxygen consumption'}; bibkey.Wie1985 = {'Wie1985'};
comment.Wie1985 = 'data from log scales, --> maybe not very precise';
temp.Wie1985 = C2K(Wie1985(:,1)); units.temp.Wie1985 = 'K' ;  label.temp.Wie1985 = 'temperature' ; 
% temp.tW_ind = 'K'; label.temp.tW_ind = 'temperature';


% %% Grouped plots
% set1 = {'WwJO_1','WwJO_2'}; comment1 = {'O2 uptake SSAF, LSAF'};
% set2 = {'tWw_1','tWw_2'}; comment2 = {'wet weight SSAF, LSAF'};
% set3 = {'tL_1','tL_2'}; comment3 = {'fork length SSAF, LSAF'};
% %set4 = {'LJO5','LJO15'}; comment4 = {'O2 uptake 5 and 15C'};
% metaData.grp.sets = {set1, set2, set3}; metaData.grp.comment = {comment1, comment2, comment3};
% %metaData.grp.sets = {set1, set2, set3, set4}; metaData.grp.comment = {comment1, comment2, comment3, comment4};

%% set weights for all real data
weights = setweights(data, []);

weights.Tah = weights.Tah * 20; 
weights.tW1 = weights.tW1 * 10; 
weights.tL1 = weights.tL1 * 10; 
% growth does something strange after 20 months post hatch, see discussion and paper
weights.tW2(end-7:end) = weights.tW2(end-7:end) * 0; 
weights.tL2(end-7:end) = weights.tL2(end-7:end) * 0;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
auxData.forkLength = forkLength;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Facts
F1 = 'Many subspecies exist, e.g. O. m. irideus  (coastal rainbow trout), O. m. gairdneri (Columbia River redband trout)';
metaData.bibkey.F1 = 'Wiki';
F2 = 'Best culturing temp 15-16 C';
metaData.bibkey.F2 = 'YaniHisa2002';
F3 = 'Able to spawn several times, each time separated by months';
metaData.bibkey.F3 = 'Wiki';
metaData.facts = struct('F1',F1,'F2',F2,'F3',F3);
      
%% Discussion
D1 = 'fish base says that females mature after 3 years and males after 2. However there is no reference to back this up. DaviKlemm2014 observed that rapid egg growth occured after 20 months post hatch, so we assume that this coincided with puberty';
D2 = 'we did a lot of empirical adjustments to the weights which needs further study.';
D3 = 'In mod_3, I removed del_M2 and extra filters; put extra weight on psd.k_J (because the used value was 0.1 * the standard one) and included psd.f_tWL = 1 to prevent that it becomes too high';
D4 = 'previous versions included weight and respiration as function of time for two extreme phenotypic families (large and small size and maturity respectively), from McKenPed2007. These families have different parameters and so should be used to make specific entries for those families. ';
metaData.discussion = struct('D1',D1, 'D2', D2, 'D3', D3, 'D4', D4);

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{https://en.wikipedia.org/wiki/Rainbow_trout}';  
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{http://www.bio.vu.nl/thb/research/bib/Kooy2010.html}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'YaniHisa2002'; type = 'Article'; bib = [ ...  
'author = {T. Yanik, S. A. Hisar and C. B�l�kbas}, ' ...
'year = {2002}, ' ...
'title = {EARLY DEVELOPMENT AND GROWTH OF ARCTIC CHARR (SALVELINUS ALPINUS) AND RAINBOW TROUT (ONCORHYNCHUS MYKISS) AT A LOW WATER TEMPERATURE.}, ' ... 
'journal = {The Israeli Journal of Aquaculture � Bamidgeh}, ' ...
'volume = {54(2)}, '...
'pages = {73}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'fishbase'; type = 'Misc'; bib = ...
'howpublished = {\url{http://www.fishbase.org/summary/239}';  
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% 
bibkey = 'ChenSnow2015'; type = 'Misc'; bib = [...
'author = {Chen, Z., Snow, M., Lawrence, C.S., Church, A.R., Narum, S.R., Devlin R.H., Farrell, A.P}, ' ...
'year = {2015}, ' ...
'title = {Selection for upper thermal tolerance in rainbow trout (Oncorhynchus mykiss Walbaum)}, ' ... 
'journal = {The Journal of Experimental Biology}, ' ...
'volume = {218}, '...
'pages = {803 - 812}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
 %   
bibkey = 'BudyThie2002'; type = 'Misc'; bib = [...
'author = { Budy, P., Thiede, G.P., Haddix, T}, ' ...
'year = {2002}, ' ...
'title = { Rainbow trout growth and survival in Flaming Gorge Reservoir. Project XIV Sport Fisheries Research (USU) Annual Report}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% 
bibkey = 'StraStut1997'; type = 'Misc'; bib = [...
'author = {Straus, D.L., Stuthridge, T.R., Anderson, S.M., Gifford, J.S.}, ' ...
'year = {1997}, ' ...
'title = {Acute toxicity of dehydroabietic acid to rainbow trout: Manipulation of biotransformation.}, ' ... 
'journal = {Australasian Journal of Ecotoxicology}, ' ...
'volume = {3}, '...
'pages = {131 - 139}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% 
bibkey = 'WeatGill1981'; type = 'Misc'; bib = [...
'author = {Weatherly, A.H., Gill, H.S.}, ' ...
'year = {1981}, ' ...
'title = {Recovery growth following periods of restricted rations and starvation in rainbow trout Salmo gairdneri Richardson}, ' ... 
'journal = {J. Fish Biol.}, ' ...
'volume = {18}, '...
'pages = {195 - 208}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% 
bibkey = 'FromRasm1991'; type = 'Misc'; bib = [...
'author = {From, J., Rasmussen, G.}, ' ...
'year = {1991}, ' ...
'title = {Growth of rainbow trout, Oncorhynchus mykiss (Walbaum, 1792) related to egg size and temperature}, ' ... 
'journal = {Dana}, ' ...
'volume = {9}, '...
'pages = {31 - 38}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Vels1987'; type = 'Book'; bib = [ ...  
'author = {F. P. J. Velsen}, ' ...
'year = {1987}, ' ...
'title = {Temperature and Incubation in Pacific Salmon and Rainbow Trout: Compilation of Data on Median Hatching Time, Mortality and Embryonic Staging}, '];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'DaviKenn2014'; type = 'Article'; bib = [ ...  
'author = {J. W. Davidson, P. B. Kenney, M. Manor, C. M. Good, G. M. Weber, A. Aussanasuwannakul, P. J. Turk, C. Welsh, S. T. Summerfelt}, ' ...
'year = {2014}, ' ...
'title = {Growth performance, fillet quality, and reproductive maturity of Rainbow Trout (Oncorhynchus mykiss) cultured to 5 kilograms within freshwater recirculating systems}, ' ... 
'journal = {Journal of Aquaculture Research and Development}, ' ...
'volume = {5}, '...
'number = {4},' ...
'pages = {}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'NinnStev2006'; type = 'Article'; bib = [ ...  
'author = {Ninness, Marcie M. and  Stevens, E. Don and Wright, Patricia A.}, ' ...
'year = {2006}, ' ...
'title = {Removal of the chorion before hatching results in increased movement and accelerated growth in rainbow trout (Oncorhynchus mykiss) embryos}, ' ... 
'journal = {Journal of Experimental Biology}, ' ...
'volume = {209}, '...
'number = {10},' ...
'pages = {1874-1882}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'KieAls1998'; type = 'Article'; bib = [ ...  
'author = {Kieffer, Alsop and  Wood}, ' ...
'year = {1998}, ' ...
'title = {A respirometric analysis of fuel use during aerobic swimming at different temperatures in rainbow trout (Oncorhynchus mykiss)}, ' ... 
'journal = {Journal of Experimental Biology}, ' ...
'volume = {201}, '...
'number = {22},' ...
'pages = {3123-3133}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Wie1985'; type = 'Article'; bib = [ ...  
'author = {Wieser}, ' ...
'year = {1985}, ' ...
'title = {Developmental and metabolic constraints of the scope for activity in young rainbow trout (Salmo Gairdneri)}, ' ... 
'journal = {Journal of Experimental Biology}, ' ...
'volume = {118}, '...
'number = {1},' ...
'pages = {133-142}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
