function [data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss

%% set metadata
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Salmoniformes'; 
metaData.family     = 'Salmonidae';
metaData.species    = 'Oncorhynchus_mykiss'; 
metaData.species_en = 'rainbow trout'; 
metaData.T_typical  = C2K(15.5); % K, body temp
metaData.data_0     = {'ah'; 'ab'; 'ap'; 'am'; 'Lb'; 'Lp'; 'Li'; 'Wwi'; 'Ri'};  % tags for different types of zero-variate data
metaData.data_1     = {'L-Ww'}; % tags for different types of uni-variate data

metaData.COMPLETE = 2.4; % using criteria of LikaKear2011

metaData.author   = {'Bas Kooijman'};        
metaData.date_acc = [2014 09 26];                           
metaData.email    = {'bas.kooijman@vu.nl'};                 
metaData.address  = {'VU University Amsterdam'}; 

%% set data
% zero-variate data
data.ah_8_5 = 33;      units.ah_8_5 = 'd';    label.ah_8_5 = 'age at hatch';           bibkey.ah_8_5 = 'YaniHisa2002'; % 4-5 month 
  temp.ah_8_5 = C2K(8.5);  units.temp.ah_8_5 = 'K'; label.temp.ah_8_5 = 'temperature';
  comment.ah_8_5 = '30-36 d';
data.ab_8_5 = data.ah_8_5 + 20; units.ab_8_5 = 'd'; label.ab_8_5 = 'age at birth';        bibkey.ab_8_5 = 'YaniHisa2002'; 
  temp.ab_8_5 = C2K(7); units.temp.ab_8_5 = 'K'; label.temp.ab_8_5 = 'temperature';
  comment.ab_8_5 = 'ah + 19-21 d';
data.ap = 2.5*365; units.ap = 'd';    label.ap = 'age at puberty';         bibkey.ap = 'fishbase';
  temp.ap = C2K(5); units.temp.ap = 'K'; label.temp.ap = 'temperature';
data.am = 11*365;  units.am = 'd';    label.am = 'life span';              bibkey.am = {'fishbase'};   
  temp.am = C2K(5); units.temp.am = 'K'; label.temp.am = 'temperature';

%data.L0 = 0.45;    units.L0 = 'cm';   label.L0 = 'egg diameter';           bibkey.Lb = 'Wiki'; 
data.Lb = 2;       units.Lb = 'cm';   label.Lb = 'total length at birth';  bibkey.Lb = 'Kooy2014'; 
data.Lp = 15;      units.Lp = 'cm';   label.Lp = 'total length at puberty';bibkey.Lp = 'fishbase';
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
 
% T-ah data from From1991
% given as the 50% value
data.Tah = [... % Temperature (K), age at hatch (d) --> fertilization to hatching 
5	61
10 31.3
];
units.Tah = {'K', 'd'};  label.Tah = {'Temperature', 'age at hatch'};  bibkey.Tah = {'From1991'};
comment.Tah = 'age is given as the 50% value';
data.Tah(:,1) = C2K(data.Tah(:,1)); % convert C to K

% T-ab data from From1991
% given as the 50% value
data.Tab = [... % Temperature (K), age at birth (d) --> fertilization to yolk absorption 
5	113
10 53.3
];
units.Tab = {'K', 'd'};  label.Tab = {'Temperature', 'age at birth'};  bibkey.Tab = {'From1991'}; 
comment.Tab = 'age is given as the 50% value';
data.Tab(:,1) = C2K(data.Tab(:,1)); % convert C to K

% Davidson2014
% Colums of tLW:
%  1  t days post hatch
%  2  L cm, fork length 
%  3 W g, wet weight
%  4 T degC, not constant

tLW = [...
30	NaN 1	12
61	NaN	2	13.9
91	NaN	12	14
122	NaN	37	13.5
152	NaN	68	13.4
183	20.14568672	130	13.3
213	21.5443469	182	12.5
244	25.04752863	297	12.7
274	30.37304642	566	13.5
305	32.15692038	685	13.6
335	35.37616057	943	14.8
366	38.472504	1230	15.2
396	41.5666287	1501	13
427	43.30037284	1713	13
457	45.39591421	2002	13.1
488	47.81766661	2307	13
518	49.18429576	2570	12.7
549	51.38700584	2836	12.6
579	53.4819702	3136	12.5
610	55.32407387	3556	12.5
640	57.17912711	4038	12.5
671	59.42611275	4533	12.7
701	60.72011595	4858	12.9
732	61.56382501	4970	13.4
762	62.32757629	5012	13.6
793	61.98695269	5097	13.4];

data.tL_Davidson2014 = tLW(:,[1 2]) ; data.tW_Davidson2014 = tLW(:,[1 3]) ; auxData.tT_Davidson2014 = tLW(:,[1 4]); auxData.tT_Davidson2014(:,2) = C2K(auxData.tT_Davidson2014(:,2));
units.tL_Davidson2014 = {'d', 'cm'}; units.tW_Davidson2014 = {'d', 'g'} ;
label.tL_Davidson2014 = {'time', 'length'}; label.tW_Davidson2014 = {'time', 'wet weight'} ;
bibkey.tL_Davidson2014 = {'Davidson2014'}; bibkey.tW_Davidson2014 = {'Davidson2014'} ;
comment.tL_Davidson2014 = 'fish reared in water recirculating system'; comment.tW_Davidson2014 = 'fish reared in water recirculating system' ;

%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Grouped plots
set1 = {'Tah','Tab'};
comment1 = {'Age at hatch, at birth'};
metaData.grp.sets = {set1};
metaData.grp.comment = {comment1};

%% Facts
F1 = 'Many subspecies exist, e.g. O. m. irideus  (coastal rainbow trout), O. m. gairdneri (Columbia River redband trout)';
metaData.bibkey.F1 = 'Wiki';
F2 = 'Best culturing temp 15-16 C';
metaData.bibkey.F2 = 'YaniHisa2002';
F3 = 'Able to spawn several times, each time separated by months';
metaData.bibkey.F3 = 'Wiki';
metaData.facts = struct('F1',F1,'F2',F2,'F3',F3);
                                 
%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/wiki/Oncorhynchus_mykiss}';  
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
bibkey = 'From1991'; type = 'Article'; bib = [ ...  
'author = {J. From, G. Rasmussen}, ' ...
'year = {1991}, ' ...
'title = {Growth of rainbow trout, Oncorhynchus mykiss (Walbaum, 1792) related to egg size and temperature}, ' ... 
'journal = {Dana}, ' ...
'volume = {9}, '...
'pages = {31-38}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Davidson2014'; type = 'Article'; bib = [ ...  
'author = {J. W. Davidson, P. B. Kenney, M. Manor, C. M. Good, G. M. Weber, A. Aussanasuwannakul, P. J. Turk, C. Welsh, S. T. Summerfelt}, ' ...
'year = {2014}, ' ...
'title = {Growth performance, fillet quality, and reproductive maturity of Rainbow Trout (Oncorhynchus mykiss) cultured to 5 kilograms within freshwater recirculating systems}, ' ... 
'journal = {Journal of Aquaculture Research and Development}, ' ...
'volume = {5(4)}, '...
'pages = {}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'fishbase'; type = 'Misc'; bib = ...
'howpublished = {\url{http://www.fishbase.org/summary/239}';  
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2014'; type = 'Misc'; bib = ...
'note = {taken from from Salmo trutta}';  
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
  
% we need the reference for From1991