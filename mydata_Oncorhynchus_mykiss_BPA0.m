function [data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss_BPA0

%% set metadata
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Salmoniformes'; 
metaData.family     = 'Salmonidae';
metaData.species    = 'Oncorhynchus_mykiss';    % previously called Salmo gairdneri (see Billard 1989)
metaData.species_en = 'rainbow trout'; 
metaData.T_typical  = C2K(15.5); % K, body temp
metaData.data_0     = {'ah'; 'ab'; 'ap'; 'am'; 'Lb'; 'Lp'; 'Li'; 'Wwi'; 'Ri'};  % tags for different types of zero-variate data
metaData.data_1     = {'L-Ww'}; % tags for different types of uni-variate data

metaData.COMPLETE = 2.4; % using criteria of LikaKear2011

metaData.author   = {'Bastien Sadoul, Starrlight Augustine'};        
metaData.date_subm = [2016 01 27];                           
metaData.email    = {'bastien.sadoul@hotmail.fr'};                 
metaData.address  = {'University of Calgary'}; 

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
 auxData.initWeight.tW = 1.471; units.initWeight.tW = 'g'; label.initWeight.tW = 'initial weight'; 
 
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

data.tL_Davidson2014 = tLW(5:end,[1 2]) ;  % I think we would run into problems if we compare predictions with NaN values
data.tW_Davidson2014 = tLW(:,[1 3]) ; 
temp.tT_Davidson2014 = tLW(:,[1 4]);  % for the auxData we need to put temp.dataLabel or it was other info we could use another name than temp
temp.tT_Davidson2014(:,2) = C2K(temp.tT_Davidson2014(:,2));
units.tL_Davidson2014 = {'d', 'cm'}; units.tW_Davidson2014 = {'d', 'g'} ; 
units.temp.tW_Davidson2014 = {'d', 'K'} ;  % we need units and label for temp info ...
label.tL_Davidson2014 = {'days post hatch', 'fork length'}; label.tW_Davidson2014 = {'days post hatch', 'wet weight'} ; label.temp.tW_Davidson2014 = {'days post hatch', 'K'} ;
bibkey.tL_Davidson2014 = {'Davidson2014'}; bibkey.tW_Davidson2014 = {'Davidson2014'} ;
comment.tL_Davidson2014 = 'fish reared in water recirculating system'; comment.tW_Davidson2014 = 'fish reared in water recirculating system' ;

% t-Wd-Wdyolk from Ninness 2006
% age dpf, dry W mg of embryo, dry W mg of yolk
% tWd_Ninness=[...
% 24	1.01	32.74;
% 30	3.08	31.8;
% 45	12.6	15.55;
% 50	22.91	8.24;
% 60	46.55	2.53;
% 75	118.93	0;
% 90	230.55	0];
% 
% tWd_Ninness(:,[2 3]) = tWd_Ninness(:,[2 3])/1000; % convert mg to g
% 
% data.tWd_Ninness=tWd_Ninness(:, [1 2]);
% units.tWd_Ninness={'d', 'g'} ; label.tWd_Ninness={'time since fecundation', 'dry weight'};
% bibkey.tWd_Ninness={'Ninness2006'};



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
 357.0      109.1714286      0.9113221      115.9836066      0.9414283      114.3197279      0.8120572
 357.5      119.2000000      0.9113221      124.8000000      0.9414283      118.0000000      0.8120572
 385.0      163.6111111      0.9113221      163.5000000      0.9414283      167.5675676      0.7958161
 412.0      187.9166667      0.9113221      189.1250000      0.9414283      191.7567568      0.7958161
];

tW_gw150(:,1)=tW_gw150(:,1)+64;         % to put in dpf

data.tW_gw150A = tW_gw150(:,[1 2]);
data.tW_gw150B = tW_gw150(:,[1 4]);
data.tW_gw150C = tW_gw150(:,[1 6]);

units.tW_gw150A = {'d', 'g'};  label.tW_gw150A = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150A = {'gw150A'};
units.tW_gw150B = {'d', 'g'};  label.tW_gw150B = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150B = {'gw150B'};
units.tW_gw150C = {'d', 'g'};  label.tW_gw150C = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150C = {'gw150C'};

temp.tW_gw150A = C2K(8.5); units.temp.tW_gw150A = 'K'; label.temp.tW_gw150A = 'temperature';
temp.tW_gw150B = C2K(8.5); units.temp.tW_gw150B = 'K'; label.temp.tW_gw150B = 'temperature';
temp.tW_gw150C = C2K(8.5); units.temp.tW_gw150C = 'K'; label.temp.tW_gw150C = 'temperature';

 
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
  327                 NaN        0.8741462                 NaN        0.8315408                 NaN        0.8936575
];

data.tW_gw124iniA = tW_gw124ini(1:14,[1 2]);
data.tW_gw124iniB = tW_gw124ini(1:14,[1 4]);
data.tW_gw124iniC = tW_gw124ini(1:14,[1 6]);

units.tW_gw124iniA = {'d', 'g'};  label.tW_gw124iniA = {'age since first feeding', 'wet weight'};  bibkey.tW_gw124iniA = {'gw124iniA'};
units.tW_gw124iniB = {'d', 'g'};  label.tW_gw124iniB = {'age since first feeding', 'wet weight'};  bibkey.tW_gw124iniB = {'gw124iniA'};
units.tW_gw124iniC = {'d', 'g'};  label.tW_gw124iniC = {'age since first feeding', 'wet weight'};  bibkey.tW_gw124iniC = {'gw124iniA'};

temp.tW_gw124iniA = C2K(8.5); units.temp.tW_gw124iniA = 'K'; label.temp.tW_gw124iniA = 'temperature';
temp.tW_gw124iniB = C2K(8.5); units.temp.tW_gw124iniB = 'K'; label.temp.tW_gw124iniB = 'temperature';
temp.tW_gw124iniC = C2K(8.5); units.temp.tW_gw124iniC = 'K'; label.temp.tW_gw124iniC = 'temperature';

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
   1013.5 2118.02326 0.8011318
];

data.tW_gw124fin = tW_gw124fin;
units.tW_gw124fin = {'d', 'g'};  label.tW_gw124fin = {'age since first feeding', 'wet weight'};  bibkey.tW_gw124fin = {'gw124fin'};
 temp.tW_gw124fin = C2K(8.5); units.temp.tW_gw124fin = 'K'; label.temp.tW_gw124fin = 'temperature';

 
 
 
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
