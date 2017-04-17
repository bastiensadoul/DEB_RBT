function [data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss_SSAF

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

metaData.author   = {'Bastien Sadoul';'Starrlight Augustine' };        
metaData.date_subm = [2017 04 15];                           
metaData.email    = {'bastien.sadoul@hotmail.fr';'starrlight.augustine@akvaplan.niva.no'};                 
metaData.address  = {'University of Calgary';'Akvaplan-niva'};

% metaData.curator     = {'Bas Kooijman'};
% metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
% metaData.date_acc    = [2017 04 15]; 

%% set data
% zero-variate data
data.Wd0 =;  units.Wd0 = 'g'; label.Wd0 = ''; bibkey.Wd0 = '';   
comment.Wd0 = ''; 

data.ah = ;  units.ah = 'd'; label.ah = 'age at hatch'; bibkey.ah = '';   
temp.ah = C2K(); units.temp.ah = 'K'; label.temp.ah = 'temperature';

data.Wdh = ;  units.Wdh = 'g'; label.Wdh = 'weight at hatch'; bibkey.Wdh = '';   
comment.Wdh = '';

data.ab = ; units.ab = 'd'; label.ab = 'age at birth'; bibkey.ab = '';   
temp.ab = C2K(); units.temp.ab = 'K'; label.temp.ab = 'temperature';

data.Wdb = ;  units.Wdb = 'g'; label.Wdb = 'weight at birth'; bibkey.Wdb = '';   
comment.Wdb = 'large eggs, 10 deg C, Table 2, wet weight times percent dry matter';

data.ap = ; units.ap = 'd';    label.ap = 'age at puberty';         bibkey.ap = '';
  temp.ap = C2K(13); units.temp.ap = 'K'; label.temp.ap = 'temperature';
  comment.ap = '';

data.am = 11*365;  units.am = 'd';    label.am = 'life span';              bibkey.am = 'fishbase';   
  temp.am = C2K(5); units.temp.am = 'K'; label.temp.am = 'temperature';

data.Lp = ;      units.Lp = 'cm';   label.Lp = 'length at puberty'; bibkey.Lp = '';

data.Wp = ;   units.Wp = '';    label.Wp = '';    bibkey.Wp = {''};

data.Li = 120;     units.Li = 'cm';   label.Li = 'ultimate total length';  bibkey.Li = 'fishbase';
data.Wi = 25400;   units.Wi = 'g';    label.Wi = 'ultimate wet weight';    bibkey.Wi = 'fishbase';
data.Ri = data.Wi * 2.5/ 365; units.Ri = '#/d'; label.Ri = 'maximum reprod rate'; bibkey.Ri = 'Wiki';   
  temp.Ri = C2K(5); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
  comment.Ri = '2000 till 3000 eggs per kg';
  


% McKenPed2007 for small size-at-age family (SSAF)
% Colums of wO2:
%  1 W g, wet weight
%  2 mean O2 uptake per day, mmol/d
WwJO = [...
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

data.WwJO=WwJO; units.WwJO = {'g', 'mmol/d'}; label.WwJO  = {'wet weight', 'O2 uptake'}; bibkey.WwJO = 'McKenPed2007';
comment.WwJO = 'Water flow of 0.55BL/s, probably at the 182g. At 213.93g food was withdrawn for 2 days.';
temp.WwJO = C2K(14); units.temp.WwJO = 'K' ;  label.temp.WwJO = 'mean temperature' ; 

tWw = [ ...
0   	76.5
21	110.3
42	156.6
63.0	214.4
85	307.8
];
data.tWw = tWw; units.tWw = {'d', 'g'}; label.tWw  = {'time', 'wet weight'}; bibkey.tWw = {'McKenPed2007'};
comment.tWw = 'SSAF - small size at age family. Diamonds fig. 1';
temp.tWw = C2K(14); units.temp.tWw = 'K' ;  label.temp.tWw = 'mean temperature' ; 

tL = [...
 42 1.45
 63 1.5
 85 1.6
];
tL(:,2) = (100 .* tWw(3:end,2) ./ tL(:,2)).^(1/3);
data.tL = tL; units.tL = {'d', 'CM'}; label.tL  = {'time', 'fork length'}; bibkey.tL = {'McKenPed2007'};
comment.tL = 'SSAF - small size at age family. Value computed from CF pp 282';
temp.tL = C2K(14); units.temp.WwL = 'K' ;  label.temp.tL = 'mean temperature' ; 


% % %% Grouped plots
% set1 = {'WwJO_1','WwJO_2'}; comment1 = {'O2 uptake SSAF, LSAF'};
% metaData.grp.sets = {set1, set2, set3}; metaData.grp.comment = {comment1, comment2, comment3};

%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp    = temp;
txtData.units   = units;
txtData.label   = label;
txtData.bibkey  = bibkey;
txtData.comment = comment;

%% Facts
F1 = 'Many subspecies exist, e.g. O. m. irideus  (coastal rainbow trout), O. m. gairdneri (Columbia River redband trout)';
metaData.bibkey.F1 = 'Wiki';
F2 = 'Able to spawn several times, each time separated by months';
metaData.bibkey.F2 = 'Wiki';
metaData.facts = struct('F1',F1,'F2',F2);
      
%% Discussion
D1 = 'data from SSAF - extreme phenotype';
metaData.discussion = struct('D1',D1);

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
