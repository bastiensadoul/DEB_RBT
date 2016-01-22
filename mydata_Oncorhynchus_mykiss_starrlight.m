function [data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss

%% set metadata
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Salmoniformes'; 
metaData.family     = 'Salmonidae';
metaData.species    = 'Oncorhynchus_mykiss';    % previously called Salmo gairdneri (see Billard 1989)
metaData.species_en = 'rainbow trout'; 
metaData.T_typical  = C2K(15.5); % K, body temp
metaData.data_0     = {'ah'; 'ab'; 'ap'; 'am'; 'Lb'; 'Lp'; 'Li'; 'Wwi'; 'Ri'};  % tags for different types of zero-variate data
metaData.data_1     = {'L-Ww'; 't-Ww'; 't-Wde_E'; 't-Wde'; }; % tags for different types of uni-variate data

metaData.COMPLETE = 2.4; % using criteria of LikaKear2011

metaData.author   = {'Bastien Sadoul, Starrlight Augustine'};        
metaData.date_acc = [2015 09 01];                           
metaData.email    = {'bastien.sadoul@hotmail.fr'};                 
metaData.address  = {'University of Calgary'}; 

%% set data
% zero-variate data
data.ah = 3;      units.ah = 'd';    label.ah = 'age at hatch';           bibkey.ah = 'NinnStev2006'; % 4-5 month 
  temp.ah = C2K(10);  units.temp.ah = 'K'; label.temp.ah = 'temperature';
  comment.ah = 'authors report hatching between 30–34 d.p.f';
data.ab = 50; units.ab = 'd'; label.ab = 'age at birth';        bibkey.ab = 'NinnStev2006'; 
  temp.ab = C2K(10); units.temp.ab = 'K'; label.temp.ab = 'temperature';
  comment.ab = 'First feeding';
data.ap = 2.5*365; units.ap = 'd';    label.ap = 'age at puberty';         bibkey.ap = 'fishbase';
  temp.ap = C2K(5); units.temp.ap = 'K'; label.temp.ap = 'temperature';
data.am = 11*365;  units.am = 'd';    label.am = 'life span';              bibkey.am = 'fishbase';   
  temp.am = C2K(5); units.temp.am = 'K'; label.temp.am = 'temperature';

% data.JO_h = 1.3;   units.JO_h = 'micromol/min/g'; label.JO_h = 'respiration per gram embryo mass at hatch'; bibkey.JO_h = 'NinnStev2006'; % 4-5 month 
%   temp.JO_h = C2K(10);  units.temp.JO_h = 'K'; label.temp.JO_h = 'temperature';
%   comment.JO_h = 'fig 2 pp 1877';
  
%data.L0 = 0.45;    units.L0 = 'cm';   label.L0 = 'egg diameter';           bibkey.Lb = 'Wiki'; 
data.Lb = 2;       units.Lb = 'cm';   label.Lb = 'total length at birth';  bibkey.Lb = 'Kooy2014'; 
data.Lp = 15;      units.Lp = 'cm';   label.Lp = 'total length at puberty';bibkey.Lp = 'fishbase';
data.Li = 120;     units.Li = 'cm';   label.Li = 'ultimate total length';  bibkey.Li = 'fishbase';

data.Wi = 25400;   units.Wi = 'g';    label.Wi = 'ultimate wet weight';    bibkey.Wi = 'fishbase';
  
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
 auxData.initWeight.tW = 1.471; units.initWeight.tW = 'g'; label.initWeight.tW = 'initial weight'; 
 
% T-ah data from From1991
% given as the 50% value
data.Tah = [... % Temperature (C), age at hatch (d) --> fertilization to hatching 
5	61
10 31.3
];
units.Tah = {'K', 'd'};  label.Tah = {'Temperature', 'age at hatch'};  bibkey.Tah = {'From1991'};
comment.Tah = 'age is given as the 50% value';
data.Tah(:,1) = C2K(data.Tah(:,1)); % convert C to K

% T-ah data from From1991
% given as the 50% value
data.Tab = [... % Temperature (C), age at birth (d) --> fertilization to yolk absorption 
5	113
10 53.3
];
units.Tab = {'K', 'd'};  label.Tab = {'Temperature', 'age at birth'};  bibkey.Tab = {'From1991'}; 
comment.Tab = 'age at which 50%hatched';
data.Tab(:,1) = C2K(data.Tab(:,1)); % convert C to K

% T-ab data from Velsen 1987
% given as the 50% value

Tah_Velsen = [...
2	115
2.5	106
2.8	93
3	111
3.2	101
4.5	72.9
4.8	75
5	72
5	68
5	64
5.8	63
6.1	61
6.2	61
6.5	57.5
7	56
7	60
7.2	45
7.5	43
7.5	43
7.7	48
7.7	44
7.7	46.5
7.8	48
7.8	49
7.8	44
7.9	42
7.9	43
7.9	48
7.9	46
7.9	46
8	41
8.7	40.3
9.2	41
9.3	35
9.5	36
9.5	36
9.5	36
9.5	36
10	38
10	34
10	33
10	35.5
10.3	29.6
10.3	28
10.4	32.1
10.7	27
10.7	29.2
10.8	29.5
11.3	30
11.5	30.3
11.5	27
11.5	28
11.7	24
12	25
12.2	26
12.2	23
12.4	24
12.5	27
12.8	24.5
12.9	18
13	28
14.5	21
15	26
15	22
15.5	18
17.5	18];

Tah_Velsen(:,1)=C2K(Tah_Velsen(:,1));
data.Tah_Velsen=Tah_Velsen;
units.Tah_Velsen={'K','d'}; label.Tah_Velsen={'Temperature', 'age at hatch'};
bibkey.Tah_Velsen={'Velsen1987'};
comment.Tah_Velsen = 'age is given as the 50% value';

% Davidson2014
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

data.tL_Davidson2014 = tLW(6:end,[1 2]) ;  % I think we would run into problems if we compare predictions with NaN values
data.tW_Davidson2014 = tLW(:,[1 3]) ; 
temp.tT_Davidson2014 = tLW(:,[1 4]);  % for the auxData we need to put temp.dataLabel or it was other info we could use another name than temp
temp.tT_Davidson2014(:,2) = C2K(temp.tT_Davidson2014(:,2));
units.tL_Davidson2014 = {'d', 'cm'}; units.tW_Davidson2014 = {'d', 'g'} ; 
units.temp.tW_Davidson2014 = {'d', 'K'} ;  % we need units and label for temp info ...
label.tL_Davidson2014 = {'days post hatch', 'fork length'}; label.tW_Davidson2014 = {'days post hatch', 'wet weight'} ; label.temp.tW_Davidson2014 = {'days post hatch', 'K'} ;
bibkey.tL_Davidson2014 = {'Davidson2014'}; bibkey.tW_Davidson2014 = {'Davidson2014'} ;
comment.tL_Davidson2014 = 'fish reared in water recirculating system'; comment.tW_Davidson2014 = 'fish reared in water recirculating system' ;

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
comment.tWdE ={'Table 1, pp 1878, chorionated'};
temp.tWde = C2K(10); units.temp.tWde = 'K'; label.temp.tWde = 'temperature'; 

data.tWde_E = tWdY(:,[1 3]); % yolk dry mass, mg 
units.tWde_E = {'d','mg'} ; label.tWde_E = {'time since fertilization', 'yolk dry weight'};
bibkey.tWde_E ={'NinnStev2006'};
comment.tWde_E ={'Table 1, pp 1878, chorionated'};
temp.tWde_E = C2K(10); units.temp.tWde_E = 'K'; label.temp.tWde_E = 'temperature'; 

% Our data for control (study gw150)
%  1  t days from first feeding
%  2 W g, wet weight

tW_gw150meancontrol=[...
0  0.1282200
14  0.2334573
28  0.4303732
42  0.7740772
56  1.2671954
70  2.0529108
84  2.7651187
98  3.9946960
112  5.8676441
126  7.8636423
140 10.1615652
154 13.8897513
175 21.7316214
196 31.7072512
217 42.2058984
245 59.1872994
273 76.6984780
];

data.tW_gw150meancontrol = tW_gw150meancontrol;
units.tW_gw150meancontrol = {'d', 'g'};  label.tW_gw150meancontrol = {'age since birth', 'wet weight'};  bibkey.tW_gw150meancontrol = {'gw150meancontrol'};
 temp.tW_gw150meancontrol = C2K(8.5); units.temp.tW_gw150meancontrol = 'K'; label.temp.tW_gw150meancontrol = 'temperature';

 
 % Our data for control (study gw124bvar)
%  1  t days from first feeding
%  2 W g, wet weight

tW_gw124bvarmeancontrol=[...
0   0.1098937
13   0.1736597
34   0.4294043
59   1.2248276
76   2.1546229
97   3.8053241
118   6.4102299
139  10.6440034
160  17.1597931
181  26.4256495
202  37.3273004
223  51.5487573
244  65.2009580
286 106.1100752
];

data.tW_gw124bvarmeancontrol = tW_gw124bvarmeancontrol;
units.tW_gw124bvarmeancontrol = {'d', 'g'};  label.tW_gw124bvarmeancontrol = {'age since birth', 'wet weight'};  bibkey.tW_gw124bvarmeancontrol = {'gw124bvarmeancontrol'};
 temp.tW_gw124bvarmeancontrol = C2K(8.5); units.temp.tW_gw124bvarmeancontrol = 'K'; label.temp.tW_gw124bvarmeancontrol = 'temperature';

% Our data, weight and length individually
% 1 is day post fertilization
% 2 is length (mm)
% 3 is weight (g)

tLW_ind=[...
85	25.7	0.1354
85	26.1	0.1554
85	24.6	0.1448
85	28	0.2262
85	26	0.1425
85	29.1	0.1972
85	27.5	0.1817
85	26.7	0.1489
85	32.4	0.2906
85	30.8	0.2605
85	27.8	0.1706
85	26.5	0.1615
85	28.8	0.195
85	28.1	0.1783
85	25.7	0.1556
85	29.1	0.2655
85	30.3	0.2419
85	30.8	0.2612
112	44	0.91
112	44	0.904
112	38	0.6
112	40	0.619
112	34	0.506
112	46	0.9
112	44	0.81
112	41	0.715
112	36	0.505
112	34	0.443
112	34	0.462
112	44	0.878
112	47	1.004
112	44	0.9
112	32	0.323
112	44	0.655
112	41	0.751
112	39	0.601
112	42	0.824
112	41	0.693
112	45	0.926
112	40	0.643
112	45	0.814
112	45	0.862
112	47	0.911
112	42	0.74
112	42	0.746
112	44	0.805
112	41	0.691
112	36	0.477
112	34	0.443
112	42	0.782
112	45	0.756
112	39	0.518
112	42	0.638
112	41	0.58
112	45	0.755
112	44	0.767
112	35	0.331
112	25	0.193
112	41	0.606
112	43	0.833
112	44	0.761
112	39	0.58
112	47	0.859
112	47	0.985
112	43	0.667
112	42	0.646
112	41	0.711
112	35	0.398
112	39	0.624
112	41	0.657
112	44	0.834
112	38	0.466
112	44	0.837
112	44	0.822
112	44	0.812
112	38	0.497
112	39	0.574
112	38	0.502
112	45	0.853
112	44	0.799
112	47	0.967
112	42	0.642
112	37	0.442
112	40	0.633
112	40	0.531
112	40	0.542
112	35	0.379
112	26	0.186
112	41	0.673
112	42	0.653
112	42	0.695
112	42	0.64
112	47	0.949
112	47	0.918
112	44	0.837
112	39	0.634
112	47	0.961
112	40	0.65
112	38	0.543
112	35	0.63
112	48	0.985
112	44	0.795
112	42	0.695
112	39	0.601
112	42	0.642
112	43	0.67
112	34	0.526
112	44	0.819
112	43	0.743
112	42	0.626
112	44	0.847
112	43	0.81
112	45	0.819
112	43	0.769
112	40	0.617
112	41	0.624
112	40	0.483
112	42	0.687
112	41	0.586
112	45	0.782
112	45	0.82
112	40	0.415
112	42	0.725
112	40	0.636
112	43	0.771
112	35	0.398
112	36	0.508
112	39	0.531
112	38	0.513
112	43	0.78
112	43	0.694
112	34	0.46
112	39	0.467
112	40	0.532
112	43	0.756
112	45	0.791
112	44	0.684
112	41	0.617
140	65	2.9
140	63	2.03
140	71	3.38
140	63	2.13
140	61	2.27
140	74	4.11
140	70	2.96
140	67	2.77
140	58	1.81
140	66	2.92
140	58	1.75
140	65	2.44
140	66	2.59
140	59	1.96
140	59	1.89
140	65	2.33
140	64	2.35
140	55	1.97
140	60	1.98
140	65	2.43
140	61	1.87
140	63	2.07
140	67	2.89
140	58	1.68
140	62	2.31
140	60	1.79
140	57	1.54
140	70	3.06
140	60	1.9
140	61	2.03
140	58	1.82
140	65	2.17
140	65	2.49
140	60	1.96
140	57	1.68
140	64	2.27
140	56	1.57
140	63	2.2
140	69	3.2
140	65	2.63
140	64	2.6
140	63	2.67
140	68	2.93
140	71	3.11
140	65	2.41
140	60	1.98
140	64	2.79
140	62	2.01
140	60	1.85
140	57	1.67
140	59	1.93
140	62	2.42
140	51	1.33
140	55	1.49
140	59	1.87
140	43	0.81
365	212	94.59
365	195	76.24
365	210	103.3
365	166	42.78
365	210	96.26
365	166	45.35
365	220	118.91
365	190	73.36
365	210	95.51
365	195	80.78
365	155	72.97
365	180	65.12
365	204	88.91
365	196	84.93
365	234	123.97
365	196	80.71
365	197	77.31
365	214	108.75
365	157	37.92
365	205	94.09
365	200	78.34
365	200	92.81
365	195	76.01
365	210	93.72
365	210	103.21
365	180	55.31
365	149	30.79
365	178	60.83
365	212	104.59
365	195	81.08
365	205	94.95
365	194	86.46
365	195	82.55
365	185	64.72
365	210	107.06
365	215	109.44
720	430	1109
720	453	1241
720	390	772
720	445	1225.7
720	410	956.5
720	429	1032.7
720	410	965
720	404	770
720	424	1126.2
720	455	1155.4
720	375	634.7
720	425	934.1
720	474	1478
720	446	1272
720	460	1682
720	417	1152
720	375	675
720	420	952
720	447	1180
720	452	1425
720	409	1061
720	444	1430
720	400	945
720	400	847
];

data.tL_ind = tLW_ind(:,[1 2]);
data.tL_ind(:,2) = data.tL_ind(:,2)/10;   % transform in cm 
units.tL_ind = {'d', 'cm'};  label.tL_ind = {'age since fertilization', 'total length'};  bibkey.tL_ind = {'individual data'};
 temp.tL_ind = C2K(8.5); units.temp.tL_ind = 'K'; label.temp.tL_ind = 'temperature';

data.tW_ind = tLW_ind(:,[1 3]);
units.tW_ind = {'d', 'g'};  label.tW_ind = {'age since fertilization', 'weight'};  bibkey.tW_ind = {'individual data'};
 temp.tW_ind = C2K(8.5); units.temp.tW_ind = 'K'; label.temp.tW_ind = 'temperature';

%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);


weights.Tah_Velsen = weights.Tah_Velsen * 5;
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
bibkey = 'Velsen1987'; type = 'Book'; bib = [ ...  
'author = {F. P. J. Velsen}, ' ...
'year = {1987}, ' ...
'title = {Temperature and Incubation in Pacific Salmon and Rainbow Trout: Compilation of Data on Median Hatching Time, Mortality and Embryonic Staging}, '];
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
