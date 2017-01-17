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
metaData.date_mod_4 = [2016 04 27];                           
metaData.email_mod_4    = {'bastien.sadoul@hotmail.fr';'starrlight.augustine@akvaplan.niva.no'};                 
metaData.address_mod_4  = {'University of Calgary';'akvaplan-niva'};
% we added respiration data and re-did parameter estimation

metaData.curator     = {'Starrlight Augustine'};
metaData.email_cur   = {'starrlight.augustine@akvaplan.niva.no'}; 
metaData.date_acc    = [2016 02 04]; 

%% set data
% zero-variate data
data.Wd0 = 0.1007 * 0.4111;  units.Wd0 = 'g'; label.Wd0 = 'dry egg weight'; bibkey.Wd0 = {'FromRasm1991'};   
comment.Wd0 = 'wet weight multiplied by percent dry matter of large eggs, see Table 1'; 

data.ah = 19 + 13;  units.ah = 'd'; label.ah = 'age at hatching at 10 degrees'; bibkey.ah = {'FromRasm1991'};   
temp.ah = C2K(10); units.temp.ah = 'K'; label.temp.ah = 'temperature';

data.ah_5 = 33 + 34;  units.ah_5 = 'd'; label.ah_5 = 'age at hatching at 5 degrees'; bibkey.ah_5 = {'FromRasm1991'};   
temp.ah_5 = C2K(5); units.temp.ah_5 = 'K'; label.temp.ah_5 = 'temperature';

data.Wdh = 0.1101 * 0.3523;  units.Wdh = 'g'; label.Wdh = 'weight at hatch'; bibkey.Wdh = {'FromRasm1991'};   
comment.Wdh = 'large eggs, 10 deg C, Table 2, wet weight times percent dry matter';

% data.JO_h = 1.3;   units.JO_h = 'micromol/min/g'; label.JO_h = 'respiration per gram embryo mass at hatch'; bibkey.JO_h = 'NinnStev2006'; % 4-5 month 
%   temp.JO_h = C2K(10);  units.temp.JO_h = 'K'; label.temp.JO_h = 'temperature';
%   comment.JO_h = 'fig 2 pp 1877';

data.ab = 19 + 13 + 22; units.ab = 'd'; label.ab = 'age at birth at 10 degrees'; bibkey.ab = {'FromRasm1991'};   
temp.ab = C2K(10); units.temp.ab = 'K'; label.temp.ab = 'temperature';

data.ab_5 = 33 + 34 + 52; units.ab_5 = 'd'; label.ab_5 = 'age at birth at 5 degrees'; bibkey.ab_5 = {'FromRasm1991'};   
temp.ab_5 = C2K(5); units.temp.ab_5 = 'K'; label.temp.ab_5 = 'temperature';

data.Wdb = 0.1665 * 0.1929;  units.Wdb = 'g'; label.Wdb = 'weight at birth'; bibkey.Wdb = {'FromRasm1991'};   
comment.Wdb = 'large eggs, 10 deg C, Table 2, wet weight times percent dry matter';

data.ap =28 + (20 * 30); units.ap = 'd';    label.ap = 'age at puberty';         bibkey.ap = 'DaviKenn2014';
  temp.ap = C2K(13); units.temp.ap = 'K'; label.temp.ap = 'temperature';
  comment.ap = 'GSI assessment indicated that rapid egg growth started at 20 Mo post hatch, and we assume that it takes about 28 d to hatch according to Vels1987';

data.am = 11*365;  units.am = 'd';    label.am = 'life span';              bibkey.am = {'fishbase'};   
  temp.am = C2K(5); units.temp.am = 'K'; label.temp.am = 'temperature';

data.Lp = 54;      units.Lp = 'cm';   label.Lp = 'forked length at puberty'; bibkey.Lp = 'DaviKenn2014';

data.Wp = 3.5 * 1e3;   units.Wp = 'g';    label.Wp = 'wet weight at puberty';    bibkey.Wp = {'DaviKenn2014'};
comment.Wp = 'weigth at 20 montsh post hatch, when rapid egg growth occurred'; 

data.Li = 120;     units.Li = 'cm';   label.Li = 'ultimate total length';  bibkey.Li = {'fishbase'};

data.Wi = 25400;   units.Wi = 'g';    label.Wi = 'ultimate wet weight';    bibkey.Wi = {'fishbase'};
  
data.Ri = data.Wi * 2.5/ 365; units.Ri = '#/d'; label.Ri = 'maximum reprod rate'; bibkey.Ri = {'Wiki'};   
  temp.Ri = C2K(5); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
  comment.Ri = '2000 till 3000 eggs per kg';
  
  
% uni-variate data
% t-Ww data from YaniHisa2002 at T = 273 + 8.5
% initial weight 1.54 g
data.tW = [... % time (d), wet weight (g)
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
units.tW = {'d', 'g'};  label.tW = {'time', 'wet weight'};  bibkey.tW = {'YaniHisa2002'};
 temp.tW = C2K(8.5); units.temp.tW = 'K'; label.temp.tW = 'temperature';
auxData.W0.tW = 1.471;  units.W0.tW = 'g';  label.W0.tW = 'wet weight at t = 0'; 

 % L_Ww data
 data.LWw = [... % length (cm), wet weight (g)
5.859	2.310; 9.567	9.342; 12.167	14.040; 12.552	43.337
12.889	7.340; 13.178	29.876; 13.756	32.171; 14.334	38.968
16.597	61.652; 24.928	170.345; 25.217	111.830; 27.095	195.273
27.769	251.609; 27.817	213.339; 28.780	224.668; 28.925	179.651
29.888	260.774; 30.610	323.867; 30.658	449.950; 31.332	368.949
32.055	308.216; 32.681	317.268; 33.355	362.347; 35.136	484.057
37.207	504.475; 37.400	738.636; 39.326	612.701; 40.722	716.371
42.119	824.543; 42.648	705.258; 43.274	813.373; 43.467	887.684
43.660	1009.274; 45.345	1173.753; 46.116	1009.459; 48.138	1095.164; 49.920	1376.723];
units.LWw = {'cm', 'g'};  label.LWw = {'length', 'wet weight'};  bibkey.LWw = {'ChenSnow2015', 'BudyThie2002',  'StraStut1997',  'WeatGill1981'};
comment.LWw = 'Compiled from 4 publications';

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
data.tL  = tLW(6:end,[1 2]); units.tL = {'d', 'cm'}; label.tL  = {'days post hatch', 'fork length'}; bibkey.tL = {'DaviKenn2014'};
comment.tL = 'fish reared in water recirculating system, we use the mean temperature';
data.tWw = tLW(:,[1 3]) ;   units.tWw = {'d', 'g'} ; label.tWw = {'days post hatch', 'wet weight'};  bibkey.tWw = {'DaviKenn2014'} ;
temp.tWw = mean(C2K(tLW(:, 4))); units.temp.tWw = 'K' ;  label.temp.tWw = 'mean temperature in recirculation system' ; 
comment.tWw = 'fish reared in water recirculating system, we use the mean temperature' ;

%  1 days post fertization, 
%  2 W g, wet weight
%  3 -, maturation stage
%  4 -, actual date given in the article
% tWGomeWeill=[...
% 0	NaN	'Immature'	'OctNov_year1' ;
% 441	508	'Immature'	'Jan_year3' ;
% 517	NaN	'First_signs_of_gametogenesis'	'MarApr_year3' ;
% 745	1755	'Mature'	'Nov_year3'];


% McKenPed2007 for small size-at-age family (SSAF)
% Colums of wO2:
%  1 W g, wet weight
%  2 mean O2 uptake per day, mmol/d
WwJO_1 = [...
182.022478	35.88060123;
188.53933	36.75507269;
192.359557	36.66290692;
196.179784	37.71930501;
199.325846	38.09071511;
203.146073	38.31084641;
206.516857	37.35688011;
210.337083	38.54067937;
213.932589	30.32143471;
217.752815	32.75766062;
221.797757	38.89625944;
225.393262	41.26175139;
229.662926	43.0810136;
233.483152	44.34460682;
238.202246	44.84228932;
242.022479	44.58948191;
246.516857	45.17000526;
251.348321	45.95614798;
255.617984	46.30142718;
260.345104	46.53687245;
264.919747	46.99361587;
269.951847	48.80583647;
274.831469	49.4383897;
279.711079	48.56930755;
284.743178	48.99043223;
289.927781	49.28986911;
295.11237	50.4058565;
300.601939	50.96802339];

data.WwJO_1=WwJO_1; units.WwJO_1 = {'g', 'mmol/d'}; label.WwJO_1  = {'wet weight', 'O2 uptake'}; bibkey.WwJO_1 = {'McKenPed2007'};
comment.WwJO_1 = 'Water flow of 0.55BL/s, probably at the 182g. At 213.93g food was withdrawn for 2 days.';
temp.WwJO_1 = C2K(14); units.temp.WwJO_1 = 'K' ;  label.temp.WwJO_1 = 'mean temperature' ; 

tWw_1 = [ ...
0   	76.5
21	110.3
42	156.6
63.0	214.4
85	307.8
];
data.tWw_1 = tWw_1; units.tWw_1 = {'d', 'g'}; label.tWw_1  = {'time', 'wet weight'}; bibkey.tWw_1 = {'McKenPed2007'};
comment.tWw_1 = 'SSAF - small size at age family. Diamonds fig. 1';
temp.tWw_1 = C2K(14); units.temp.tWw_1 = 'K' ;  label.temp.tWw_1 = 'mean temperature' ; 

tL_1 = [...
 42 1.45
 63 1.5
 85 1.6
];
tL_1(:,2) = (100 .* tWw_1(3:end,2) ./ tL_1(:,2)).^(1/3);
data.tL_1 = tL_1; units.tL_1 = {'d', 'CM'}; label.tL_1  = {'time', 'fork length'}; bibkey.tL_1 = {'McKenPed2007'};
comment.tL_1 = 'SSAF - small size at age family. Value computed from CF pp 282';
temp.tL_1 = C2K(14); units.temp.WwL_1 = 'K' ;  label.temp.tL_1 = 'mean temperature' ; 


% McKenPed2007 for large size-at-age family (LSAF)
% Colums of wO2:
%  1 W g, wet weight
%  2 mean O2 uptake per day, mmol/d
WwJO_2 = [...
181.966299	34.20024145;
185.016057	35.34066636;
187.608351	37.09269751;
190.505626	39.65565454;
193.402893	40.56620072;
196.147678	41.5873676;
199.044952	42.65367815;
201.94222	44.6045036;
205.144468	46.26686371;
208.194232	46.62372168;
211.0915	46.96094466;
214.293747	47.84367026;
217.343505	48.77135785;
220.545752	49.79044845;
223.900489	49.88679642;
227.255227	49.9633361;
230.609951	50.09862151;
233.812198	50.4756995;
237.471915	51.99381419;
239.911727	43.34750274;
241.894069	44.03527827;
244.638847	48.42395257;
248.298564	51.12195049;
251.805778	52.10137853;
255.770474	51.29533228;
259.430178	50.49747008;
263.394875	50.49157301;
267.054579	51.31441558;
271.171752	51.2125789;
275.136449	51.61768034;
279.253622	52.96085712;
283.523286	53.99596248;
287.792949	55.69142944;
292.367579	55.21555101;
296.332276	54.88756234;
300.906906	54.94902987];
data.WwJO_2 = WwJO_2; units.WwJO_2 = {'g', 'mmol/d'}; label.WwJO_2  = {'wet weight', 'O2 uptake'}; bibkey.WwJO_2 = {'McKenPed2007'};
comment.WwJO_2 = 'LSAF - large size at age family. Water flow of 0.55BL/s, probably at the 182g. At 239.9g food was withdrawn for 2 days.';
temp.WwJO_2 = C2K(14); units.temp.WwJO_2 = 'K' ;  label.temp.WwJO_2 = 'mean temperature' ; 

tWw_2 = [ ...
0    	181.5
21.0	242.0
42	    319.4
63.0	371.9
85   	443.1
];
data.tWw_2 = tWw_2; units.tWw_2 = {'d', 'g'}; label.tWw_2  = {'time', 'wet weight'}; bibkey.tWw_2 = {'McKenPed2007'};
comment.tWw_2 = 'LSAF - large size at age family. Circles fig. 1';
temp.tWw_2 = C2K(14); units.temp.tWw_2 = 'K' ;  label.temp.tWw_2 = 'mean temperature' ; 

tL_2 = [...
 0 1.3
 21 1.4
 42 1.5
];
tL_2(:,2) = (100 .* tWw_2(1:3,2) ./ tL_2(:,2)).^(1/3);
data.tL_2 = tL_2; units.tL_2 = {'d', 'CM'}; label.tL_2  = {'time', 'fork length'}; bibkey.tL_2 = {'McKenPed2007'};
comment.tL_2 = 'SSAF - small size at age family. Value computed from CF pp 282';
temp.tL_2 = C2K(14); units.temp.WwL_2 = 'K' ;  label.temp.tL_2 = 'mean temperature' ; 


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


%------------------------------------------------------------------------------------------------
%     DATA NOT PUBLISHED     (ONLY CONTROLS)
%------------------------------------------------------------------------------------------------

% Our data for control (study gw150b)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names: date   weight.A702_A133	surv.A702_A133	weight.A706_A137	surv.A706_A137	weight.A713_A144	surv.A713_A144

tW_gw150=[...
    0.0        0.1301370      1.0000000        0.1280822      1.0000000        0.1264407      1.0000000
   14.0        0.2308772      0.9760274        0.2238596      0.9760274        0.2456349      0.8542373
   28.0        0.4258993      0.9520548        0.4257143      0.9589041        0.4395062      0.8237288
   42.0        0.7931408      0.9486301        0.7700000      0.9589041        0.7590909      0.8203390
   56.0        1.3334545      0.9417808        1.2082143      0.9589041        1.2599174      0.8203390
   70.0        2.1609489      0.9383562        1.9246429      0.9589041        2.0731405      0.8203390
   84.0        2.9242424      0.9349315        2.6357843      0.9554795        2.7353293      0.8169492
   98.0        4.1257576      0.9349315        3.8529412      0.9554795        4.0053892      0.8169492
  112.0        6.2828283      0.9349315        5.5490196      0.9554795        5.7710843      0.8120572
 126.0        8.0085859      0.9349315        7.7745098      0.9554795        7.8078313      0.8120572
 140.0       10.4131980      0.9302096       10.2745098      0.9554795        9.7969880      0.8120572
 154.0       14.6192893      0.9302096       13.6764706      0.9554795       13.3734940      0.8120572
 175.0       22.0558376      0.9302096       21.8137255      0.9554795       21.3253012      0.8120572
 196.0       31.8622449      0.9254878       31.8137255      0.9554795       31.4457831      0.8120572
 217.0       42.5765306      0.9254878       42.0833333      0.9554795       41.9578313      0.8120572
 245.0       59.9744898      0.9254878       58.0392157      0.9554795       59.5481928      0.8120572
 273.0       76.7692308      0.9207659       76.5686275      0.9554795       76.7575758      0.8120572
 357.0       114.1857143      0.9113221      120.3918033      0.9414283      116.1598639      0.8120572
%  357.5      119.2000000      0.9113221      124.8000000      0.9414283     118.0000000      0.8120572                 %  because of cull effect
%  385.0      163.6111111      0.9113221      163.5000000      0.9414283      167.5675676      0.7958161
%  412.0      187.9166667      0.9113221      189.1250000      0.9414283      191.7567568      0.7958161
];

% data.tW_gw150(:,1)=tW_gw150(:,1)+64;         % to put in dpf
data.tW_gw150(:,1)=tW_gw150(:,1)+42;         % to put in dph
data.tW_gw150(:,2) = (1/3)*(tW_gw150(:,[2])+ tW_gw150(:,[4])+ tW_gw150(:,[6]));
units.tW_gw150 = {'d', 'g'};  label.tW_gw150 = {'age since hatch', 'wet weight'};  bibkey.tW_gw150 = {'gw150_Control'};
auxData.t0.tW_gw150  = 'dpf';

% Our data for control (study gw124b)
%  days from first feeding (don't know the dpf...)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names:   date weight.C57_C10_B17	surv.C57_C10_B17	weight.C63_C17_B17	surv.C63_C17_B17	weight.C69_C22_B17 surv.C69_C22_B17

tW_gw124ini=[...
    0          0.1181818        1.0000000          0.1026596        1.0000000          0.1088398        1.0000000
   13          0.1785276        0.9431818          0.1675978        0.9680851          0.1748538        0.9613260
   34          0.4401361        0.9142499          0.4166667        0.9194104          0.4314103        0.9388388
   59          1.2429630        0.9018111          1.1756579        0.9194104          1.2558621        0.9328206
   76          2.1733333        0.9018111          2.0960526        0.9194104          2.1944828        0.9328206
   97          3.8037037        0.9018111          3.6315789        0.9194104          3.9806897        0.9328206
  118          6.6153846        0.8951310          5.9246575        0.9073129          6.6906475        0.9135209
  139         10.6201550        0.8882454         10.0344828        0.9010985         11.2773723        0.9003767
  160         17.1317829        0.8882454         16.1724138        0.9010985         18.1751825        0.9003767
  181         26.7460317        0.8882454         25.1428571        0.8886695         27.3880597        0.9003767
  202         38.7301587        0.8882454         34.3795620        0.8696266         38.8721805        0.8936575
  223         53.3730159        0.8882454         47.6642336        0.8696266         53.6090226        0.8936575
  244         68.1200000        0.8811958         59.8888889        0.8569313         67.5939850        0.8936575
  286        112.2400000        0.8811958         96.4661654        0.8442360        109.6240602        0.8936575
%   327                 NaN        0.8741462                 NaN        0.8315408                 NaN        0.8936575
];

data.tW_gw124ini(:,1) = tW_gw124ini(:,1)+42;     % days post hatch
data.tW_gw124ini(:,2) = (1/3)*(tW_gw124ini(:,[2]) + tW_gw124ini(:,[4])+ tW_gw124ini(:,[6]));

units.tW_gw124ini = {'d', 'g'};  label.tW_gw124ini = {'age since hatch', 'wet weight'};  bibkey.tW_gw124ini = {'gw124ini_Control'};

auxData.t0.tW_gw124ini  = 'dpb';

% Our data for control (study gw124b)
%  1  t days from first feeding (don't know the dpf)
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding (surv taken using the mean of the 
%  column names:         date     weight      surv 

tW_gw124fin=[...
    375.0  235.53333 0.8664482
    431.0  358.33333 0.8664482
    488.0  522.00000 0.8664482
    552.0  724.05405 0.8548955
    552.5  688.49315 0.8548955
    608.0  816.91781 0.8548955
    664.0  901.87500 0.8431846
    720.0 1049.45652 0.8197628
    720.5 1097.15909 0.8197628
    781.0 1332.95455 0.8197628
    837.0 1688.97727 0.8197628
    894.0 1853.86364 0.8197628
    894.5 1853.86364 0.8197628
    949.0 2145.22727 0.8197628
   1005.0 2118.02326 0.8011318
%    1013.5 2118.02326 0.8011318
];

data.tW_gw124fin = tW_gw124fin(:,[1 2]);
data.tW_gw124fin(:,1) = tW_gw124fin(:,1)+42;     % days post hatch
units.tW_gw124fin = {'d', 'g'};  label.tW_gw124fin = {'age since hatch', 'wet weight'};  bibkey.tW_gw124fin = {'gw124fin_Control'};
 
auxData.t0.tW_gw124fin  = 'dpb';


% % Our data, weight and length individually
% % 1 is day post fertilization
% % 2 is length (mm)
% % 3 is weight (g)
% 
% tLW_ind=[...
% 85	25.7	0.1354
% 85	26.1	0.1554
% 85	24.6	0.1448
% 85	28	0.2262
% 85	26	0.1425
% 85	29.1	0.1972
% 85	27.5	0.1817
% 85	26.7	0.1489
% 85	32.4	0.2906
% 85	30.8	0.2605
% 85	27.8	0.1706
% 85	26.5	0.1615
% 85	28.8	0.195
% 85	28.1	0.1783
% 85	25.7	0.1556
% 85	29.1	0.2655
% 85	30.3	0.2419
% 85	30.8	0.2612
% 112	44	0.91
% 112	44	0.904
% 112	38	0.6
% 112	40	0.619
% 112	34	0.506
% 112	46	0.9
% 112	44	0.81
% 112	41	0.715
% 112	36	0.505
% 112	34	0.443
% 112	34	0.462
% 112	44	0.878
% 112	47	1.004
% 112	44	0.9
% 112	32	0.323
% 112	44	0.655
% 112	41	0.751
% 112	39	0.601
% 112	42	0.824
% 112	41	0.693
% 112	45	0.926
% 112	40	0.643
% 112	45	0.814
% 112	45	0.862
% 112	47	0.911
% 112	42	0.74
% 112	42	0.746
% 112	44	0.805
% 112	41	0.691
% 112	36	0.477
% 112	34	0.443
% 112	42	0.782
% 112	45	0.756
% 112	39	0.518
% 112	42	0.638
% 112	41	0.58
% 112	45	0.755
% 112	44	0.767
% 112	35	0.331
% 112	25	0.193
% 112	41	0.606
% 112	43	0.833
% 112	44	0.761
% 112	39	0.58
% 112	47	0.859
% 112	47	0.985
% 112	43	0.667
% 112	42	0.646
% 112	41	0.711
% 112	35	0.398
% 112	39	0.624
% 112	41	0.657
% 112	44	0.834
% 112	38	0.466
% 112	44	0.837
% 112	44	0.822
% 112	44	0.812
% 112	38	0.497
% 112	39	0.574
% 112	38	0.502
% 112	45	0.853
% 112	44	0.799
% 112	47	0.967
% 112	42	0.642
% 112	37	0.442
% 112	40	0.633
% 112	40	0.531
% 112	40	0.542
% 112	35	0.379
% 112	26	0.186
% 112	41	0.673
% 112	42	0.653
% 112	42	0.695
% 112	42	0.64
% 112	47	0.949
% 112	47	0.918
% 112	44	0.837
% 112	39	0.634
% 112	47	0.961
% 112	40	0.65
% 112	38	0.543
% 112	35	0.63
% 112	48	0.985
% 112	44	0.795
% 112	42	0.695
% 112	39	0.601
% 112	42	0.642
% 112	43	0.67
% 112	34	0.526
% 112	44	0.819
% 112	43	0.743
% 112	42	0.626
% 112	44	0.847
% 112	43	0.81
% 112	45	0.819
% 112	43	0.769
% 112	40	0.617
% 112	41	0.624
% 112	40	0.483
% 112	42	0.687
% 112	41	0.586
% 112	45	0.782
% 112	45	0.82
% 112	40	0.415
% 112	42	0.725
% 112	40	0.636
% 112	43	0.771
% 112	35	0.398
% 112	36	0.508
% 112	39	0.531
% 112	38	0.513
% 112	43	0.78
% 112	43	0.694
% 112	34	0.46
% 112	39	0.467
% 112	40	0.532
% 112	43	0.756
% 112	45	0.791
% 112	44	0.684
% 112	41	0.617
% 140	65	2.9
% 140	63	2.03
% 140	71	3.38
% 140	63	2.13
% 140	61	2.27
% 140	74	4.11
% 140	70	2.96
% 140	67	2.77
% 140	58	1.81
% 140	66	2.92
% 140	58	1.75
% 140	65	2.44
% 140	66	2.59
% 140	59	1.96
% 140	59	1.89
% 140	65	2.33
% 140	64	2.35
% 140	55	1.97
% 140	60	1.98
% 140	65	2.43
% 140	61	1.87
% 140	63	2.07
% 140	67	2.89
% 140	58	1.68
% 140	62	2.31
% 140	60	1.79
% 140	57	1.54
% 140	70	3.06
% 140	60	1.9
% 140	61	2.03
% 140	58	1.82
% 140	65	2.17
% 140	65	2.49
% 140	60	1.96
% 140	57	1.68
% 140	64	2.27
% 140	56	1.57
% 140	63	2.2
% 140	69	3.2
% 140	65	2.63
% 140	64	2.6
% 140	63	2.67
% 140	68	2.93
% 140	71	3.11
% 140	65	2.41
% 140	60	1.98
% 140	64	2.79
% 140	62	2.01
% 140	60	1.85
% 140	57	1.67
% 140	59	1.93
% 140	62	2.42
% 140	51	1.33
% 140	55	1.49
% 140	59	1.87
% 140	43	0.81
% 365	212	94.59
% 365	195	76.24
% 365	210	103.3
% 365	166	42.78
% 365	210	96.26
% 365	166	45.35
% 365	220	118.91
% 365	190	73.36
% 365	210	95.51
% 365	195	80.78
% 365	155	72.97
% 365	180	65.12
% 365	204	88.91
% 365	196	84.93
% 365	234	123.97
% 365	196	80.71
% 365	197	77.31
% 365	214	108.75
% 365	157	37.92
% 365	205	94.09
% 365	200	78.34
% 365	200	92.81
% 365	195	76.01
% 365	210	93.72
% 365	210	103.21
% 365	180	55.31
% 365	149	30.79
% 365	178	60.83
% 365	212	104.59
% 365	195	81.08
% 365	205	94.95
% 365	194	86.46
% 365	195	82.55
% 365	185	64.72
% 365	210	107.06
% 365	215	109.44
% 720	430	1109
% 720	453	1241
% 720	390	772
% 720	445	1225.7
% 720	410	956.5
% 720	429	1032.7
% 720	410	965
% 720	404	770
% 720	424	1126.2
% 720	455	1155.4
% 720	375	634.7
% 720	425	934.1
% 720	474	1478
% 720	446	1272
% 720	460	1682
% 720	417	1152
% 720	375	675
% 720	420	952
% 720	447	1180
% 720	452	1425
% 720	409	1061
% 720	444	1430
% 720	400	945
% 720	400	847
% ];
% 
% data.tL_ind = tLW_ind(:,[1 2]);
% data.tL_ind(:,2) = data.tL_ind(:,2)/10;   % transform in cm 
% units.tL_ind = {'d', 'cm'};  label.tL_ind = {'age since fertilization', 'total length'};  bibkey.tL_ind = {'individual data'};
%  temp.tL_ind = C2K(8.5); units.temp.tL_ind = 'K'; label.temp.tL_ind = 'temperature';
% 
% data.tW_ind = tLW_ind(:,[1 3]);
% units.tW_ind = {'d', 'g'};  label.tW_ind = {'age since fertilization', 'weight'};  bibkey.tW_ind = {'individual data'};
%  temp.tW_ind = C2K(8.5); units.temp.tW_ind = 'K'; label.temp.tW_ind = 'temperature';


% %% Grouped plots
set1 = {'WwJO_1','WwJO_2'}; comment1 = {'O2 uptake SSAF, LSAF'};
set2 = {'tWw_1','tWw_2'}; comment2 = {'wet weight SSAF, LSAF'};
set3 = {'tL_1','tL_2'}; comment3 = {'fork length SSAF, LSAF'};
%set4 = {'LJO5','LJO15'}; comment4 = {'O2 uptake 5 and 15C'};
metaData.grp.sets = {set1, set2, set3}; metaData.grp.comment = {comment1, comment2, comment3};
%metaData.grp.sets = {set1, set2, set3, set4}; metaData.grp.comment = {comment1, comment2, comment3, comment4};

%% set weights for all real data
weights = setweights(data, []);

% growth does something strange after 20 m post hatch, see discussion and paper
% weights.LWw = weights.LWw * 0.01; 
% weights.tW = weights.tW * 10; 
% weights.tWw = weights.tWw * 200; 
% weights.tL = weights.tL * 30; 
weights.tWw(end-7:end) = weights.tWw(end-7:end) * 0; 
weights.tL(end-7:end) = weights.tL(end-7:end) * 0;
% weights.Tah = weights.Tah * 60; % this is empirical, it just helped
% weights.tWde_E = weights.tWde_E * 200; % this is empirical, it just helped
% weights.tWde = weights.tWde * 200; % this is empirical, it just helped
% weights.Wi= weights.Wi * 80; % this is empirical, it just helped
% weights.Wd0= weights.Wd0 * 800; % this is empirical, it just helped
% % weights.WJO = weights.WJO * 50;  % this is empirical, it just helped
% % weights.LJO5 = weights.LJO5 * 50;  % this is empirical, it just helped
% % weights.LJO15 = weights.LJO15 * 50;  % this is empirical, it just helped
weights.tW_gw150= weights.tW_gw150 * 0; % Try without using non published data. this is empirical, it just helped
weights.tW_gw124ini= weights.tW_gw124ini * 0; % Try without using non published data. this is empirical, it just helped
weights.tW_gw124fin= weights.tW_gw124fin * 0; % Try without using non published data. this is empirical, it just helped


%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
weights.psd.kap = weights.psd.kap * 0;
weights.psd.k_J = weights.psd.k_J * 50;

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
metaData.discussion = struct('D1',D1, 'D2', D2, 'D3', D3);

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
'author = {T. Yanik, S. A. Hisar and C. Bölükbas}, ' ...
'year = {2002}, ' ...
'title = {EARLY DEVELOPMENT AND GROWTH OF ARCTIC CHARR (SALVELINUS ALPINUS) AND RAINBOW TROUT (ONCORHYNCHUS MYKISS) AT A LOW WATER TEMPERATURE.}, ' ... 
'journal = {The Israeli Journal of Aquaculture – Bamidgeh}, ' ...
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
