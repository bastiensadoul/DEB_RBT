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
metaData.data_1     = {'L-Ww'; 't-Ww'; 't-L'; 'tWde'; 'tWde_E'; 'T-ah'}; % tags for different types of uni-variate data

metaData.COMPLETE = 2.4; % using criteria of LikaKear2011

metaData.author   = {'Bas Kooijman'};        
metaData.date_subm = [2014 09 26];                           
metaData.email    = {'bas.kooijman@vu.nl'};                 
metaData.address  = {'VU University Amsterdam'}; 

metaData.author_mod_1   = {'Starrlight Augustine'};        
metaData.date_mod_1 = [2016 01 26];                           
metaData.email_mod_1    = {'starrlight.augustine@akvaplan.niva.no'};                 
metaData.address_mod_1  = {'akvaplan-niva'};

metaData.author_mod_2   = {'Bastien Sadoul';'Starrlight Augustine'; };        
metaData.date_mod_2 = [2016 02 01];                           
metaData.email_mod_2    = {'bastien.sadoul@hotmail.fr';'starrlight.augustine@akvaplan.niva.no'};                 
metaData.address_mod_2  = {'University of Calgary';'akvaplan-niva'};

metaData.author_mod_3   = {'Bas Kooijman'};        
metaData.date_mod_3 = [2016 04 07];                           
metaData.email_mod_3    = {'bas.kooijman@vu.nl'};                 
metaData.address_mod_3  = {'VU University Amsterdam'};

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


% % McKenPed2007 for small size-at-age family (SSAF)
% % Colums of wO2:
% %  1 W g, wet weight
% %  2 mean O2 uptake per day, mmol/d
% wO2_1 = [...
% 182.022478	35.88060123;
% 188.53933	36.75507269;
% 192.359557	36.66290692;
% 196.179784	37.71930501;
% 199.325846	38.09071511;
% 203.146073	38.31084641;
% 206.516857	37.35688011;
% 210.337083	38.54067937;
% 213.932589	30.32143471;
% 217.752815	32.75766062;
% 221.797757	38.89625944;
% 225.393262	41.26175139;
% 229.662926	43.0810136;
% 233.483152	44.34460682;
% 238.202246	44.84228932;
% 242.022479	44.58948191;
% 246.516857	45.17000526;
% 251.348321	45.95614798;
% 255.617984	46.30142718;
% 260.345104	46.53687245;
% 264.919747	46.99361587;
% 269.951847	48.80583647;
% 274.831469	49.4383897;
% 279.711079	48.56930755;
% 284.743178	48.99043223;
% 289.927781	49.28986911;
% 295.11237	50.4058565;
% 300.601939	50.96802339];
% 
% data.wO2_1=wO2_1; units.wO2_1 = {'g', 'mmol/d'}; label.wO2_1  = {'wet weight', 'O2 uptake'}; bibkey.wO2_1 = {'McKenPed2007 SSAF'};
% comment.wO2_1 = 'Water flow of 0.55BL/s, probably at the 182g. At 213.93g food was withdrawn for 2 days.';
% temp.wO2_1 = C2K(14); units.temp.wO2_1 = 'K' ;  label.temp.wO2_1 = 'mean temperature' ; 
% 
% 
% % McKenPed2007 for large size-at-age family (LSAF)
% % Colums of wO2:
% %  1 W g, wet weight
% %  2 mean O2 uptake per day, mmol/d
% wO2_2 = [...
% 181.966299	34.20024145;
% 185.016057	35.34066636;
% 187.608351	37.09269751;
% 190.505626	39.65565454;
% 193.402893	40.56620072;
% 196.147678	41.5873676;
% 199.044952	42.65367815;
% 201.94222	44.6045036;
% 205.144468	46.26686371;
% 208.194232	46.62372168;
% 211.0915	46.96094466;
% 214.293747	47.84367026;
% 217.343505	48.77135785;
% 220.545752	49.79044845;
% 223.900489	49.88679642;
% 227.255227	49.9633361;
% 230.609951	50.09862151;
% 233.812198	50.4756995;
% 237.471915	51.99381419;
% 239.911727	43.34750274;
% 241.894069	44.03527827;
% 244.638847	48.42395257;
% 248.298564	51.12195049;
% 251.805778	52.10137853;
% 255.770474	51.29533228;
% 259.430178	50.49747008;
% 263.394875	50.49157301;
% 267.054579	51.31441558;
% 271.171752	51.2125789;
% 275.136449	51.61768034;
% 279.253622	52.96085712;
% 283.523286	53.99596248;
% 287.792949	55.69142944;
% 292.367579	55.21555101;
% 296.332276	54.88756234;
% 300.906906	54.94902987];
% data.wO2_2=wO2_2; units.wO2_2 = {'g', 'mmol/d'}; label.wO2_2  = {'wet weight', 'O2 uptake'}; bibkey.wO2_2 = {'McKenPed2007 LSAF'};
% comment.wO2_2 = 'Water flow of 0.55BL/s, probably at the 182g. At 239.9g food was withdrawn for 2 days.';
% temp.wO2_2 = C2K(14); units.temp.wO2_2 = 'K' ;  label.temp.wO2_2 = 'mean temperature' ; 
% 


% %% Grouped plots
% set1 = {'tWde_E','tWde'}; comment1 = {'mg, dry weight of yolk, embryo'};
% metaData.grp.sets = {set1}; metaData.grp.comment = {comment1};

%% set weights for all real data
weights = setweights(data, []);

% growth does something strange after 20 m post hatch, see discussion and paper
weights.LWw = weights.LWw * 0.01; 
% weights.tW = weights.tW * 10; 
weights.tWw = weights.tWw * 200; 
weights.tL = weights.tL * 30; 
weights.tWw(end-7:end) = weights.tWw(end-7:end) * 0; 
weights.tL(end-7:end) = weights.tL(end-7:end) * 0; % 
weights.Tah = weights.Tah * 60; % this is empirical, it just helped
weights.tWde_E = weights.tWde_E * 200; % this is empirical, it just helped
weights.tWde = weights.tWde * 200; % this is empirical, it just helped
weights.Wi= weights.Wi * 80; % this is empirical, it just helped
weights.Wd0= weights.Wd0 * 800; % this is empirical, it just helped

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
data.psd.f_tWL = 1;  units.psd.f_tWL = '-';  label.psd.f_tWL = 'scaled functional response for puberty data';
weights.psd.f_tWL = weights.psd.kap * 20;
weights.psd.kap = weights.psd.kap * 0;
weights.psd.k_J = weights.psd.k_J * 50;

%% pack auxData and txtData for output
auxData.temp = temp;
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
