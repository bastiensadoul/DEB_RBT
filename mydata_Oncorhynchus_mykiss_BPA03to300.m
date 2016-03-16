function [data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss_BPA03to300

%% set metadata
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Salmoniformes'; 
metaData.family     = 'Salmonidae';
metaData.species    = 'Oncorhynchus_mykiss_BPA3';    % previously called Salmo gairdneri (see Billard 1989)
metaData.species_en = 'rainbow trout'; 
metaData.T_typical  = C2K(15.5); % K, body temp
metaData.data_0     = {''};  % tags for different types of zero-variate data
metaData.data_1     = {'a-Ww';'L-Ww'}; % tags for different types of uni-variate data

metaData.COMPLETE   = 1.5;        


metaData.author   = {'Bastien Sadoul, Starrlight Augustine'};        
metaData.date_subm = [2016 01 27];                           
metaData.email    = {'bastiensadoul@hotmail.fr'};                 
metaData.address  = {'University of Calgary'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set data

%---------------------------------------------Study gw150

% Our data for BPA3 (study gw150b)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names:  date weight.A703_A134 surv.A703_A134 weight.A708_A139 surv.A708_A139 weight.A712_A143 surv.A712_A143

tW_gw150_BPA3=[...
     0.0        0.1304348      1.0000000        0.1309764      1.0000000        0.1300676      1.0000000
   14.0        0.2377104      0.9933110        0.2422680      0.9797980        0.2311644      0.9864865
   28.0        0.4255102      0.9832776        0.4283737      0.9730640        0.4136364      0.9662162
   42.0        0.7924915      0.9799331        0.7768166      0.9730640        0.7088028      0.9594595
   56.0        1.3000000      0.9799331        1.3072664      0.9730640        1.2239437      0.9594595
   70.0        2.0955631      0.9799331        2.0312500      0.9696970        1.9770318      0.9560811
   84.0        2.8876147      0.9799331        2.7957547      0.9663300        2.6417476      0.9425676
   98.0        4.2460465      0.9709429        3.9622642      0.9663300        3.7752427      0.9425676
 112.0        5.9534884      0.9709429        5.9056604      0.9663300        5.3009709      0.9425676
 126.0        7.9051402      0.9709429        7.9264151      0.9663300        7.4805825      0.9425676
 140.0       10.3823810      0.9618687       10.7445498      0.9617718        9.8409756      0.9379920
 154.0       14.0000000      0.9618687       14.7630332      0.9617718       13.1463415      0.9379920
 175.0       21.9285714      0.9618687       22.5476190      0.9572136       20.7317073      0.9379920
 196.0       31.6190476      0.9618687       32.1770335      0.9526555       30.3902439      0.9379920
 217.0       41.8095238      0.9618687       42.3923445      0.9526555       39.9268293      0.9379920
 245.0       58.7380952      0.9618687       60.0478469      0.9526555       55.9313725      0.9334164
 273.0       73.7380952      0.9618687       77.6315789      0.9526555       73.4313725      0.9334164
 357.0      112.3682540      0.9481277      118.6894737      0.9526555      113.0021189      0.9196897
%  357.5      116.8000000      0.9481277      121.8000000      0.9526555    116.6326531      0.9196897                             % because of cull effect
%  385.0      154.8684211      0.9481277      153.9473684      0.9526555      162.0833333      0.9196897
%  412.0      176.9736842      0.9481277      182.5000000      0.9526555      185.1388889      0.9196897
];

tW_gw150_BPA3(:,1)=tW_gw150_BPA3(:,1)+64;         % to put in dpf

data.tW_gw150A_BPA3 = tW_gw150_BPA3(:,[1 2]);
data.tW_gw150B_BPA3 = tW_gw150_BPA3(:,[1 4]);
data.tW_gw150C_BPA3 = tW_gw150_BPA3(:,[1 6]);

units.tW_gw150A_BPA3 = {'d', 'g'};  label.tW_gw150A_BPA3 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150A_BPA3 = {'gw150A-BPA3'};
units.tW_gw150B_BPA3 = {'d', 'g'};  label.tW_gw150B_BPA3 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150B_BPA3 = {'gw150B-BPA3'};
units.tW_gw150C_BPA3 = {'d', 'g'};  label.tW_gw150C_BPA3 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150C_BPA3 = {'gw150C-BPA3'};

temp.tW_gw150A_BPA3 = C2K(8.5); units.temp.tW_gw150A_BPA3 = 'K'; label.temp.tW_gw150A_BPA3 = 'temperature';
temp.tW_gw150B_BPA3 = C2K(8.5); units.temp.tW_gw150B_BPA3 = 'K'; label.temp.tW_gw150B_BPA3 = 'temperature';
temp.tW_gw150C_BPA3 = C2K(8.5); units.temp.tW_gw150C_BPA3 = 'K'; label.temp.tW_gw150C_BPA3 = 'temperature';


% Our data for BPA30 (study gw150b)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names: date weight.A705_A136 surv.A705_A136 weight.A709_A140 surv.A709_A140 weight.A711_A142 surv.A711_A142

tW_gw150_BPA30=[...
    0.0        0.1163194      1.0000000        0.1155235      1.0000000        0.1118881      1.0000000
   14.0        0.1992674      0.9479167        0.1991935      0.8953069        0.1972868      0.9020979
   28.0        0.3555133      0.9131944        0.3693277      0.8592058        0.3653386      0.8776224
   42.0        0.6513308      0.9131944        0.6777311      0.8592058        0.6596000      0.8741259
   56.0        1.0704981      0.9062500        1.1445378      0.8592058        1.0943775      0.8706294
   70.0        1.6735632      0.9062500        1.8620253      0.8555957        1.7599190      0.8636364
   84.0        2.2680851      0.9062500        2.4761006      0.8519856        2.3616279      0.8601399
   98.0        3.3755319      0.9062500        3.5635220      0.8519856        3.4465116      0.8601399
 112.0        5.1170213      0.9062500        5.0691824      0.8519856        4.8837209      0.8601399
 126.0        6.5037234      0.9062500        7.1402516      0.8519856        6.7732558      0.8601399
 140.0        8.4465241      0.9014295        9.6331210      0.8412688        9.1639535      0.8601399
 154.0       11.3903743      0.9014295       12.9166667      0.8359104       11.9590643      0.8551390
 175.0       18.3957219      0.9014295       20.0641026      0.8359104       18.3625731      0.8551390
 196.0       26.7379679      0.9014295       28.9423077      0.8359104       25.9064327      0.8551390
 217.0       36.1497326      0.9014295       38.5576923      0.8359104       35.2631579      0.8551390
 245.0       51.0695187      0.9014295       53.7500000      0.8359104       49.9707602      0.8551390
 273.0       69.5989305      0.9014295       70.3205128      0.8359104       66.5497076      0.8551390
 357.0      108.9285714      0.8966090      116.4978261      0.8359104      105.1212418      0.8551390
%  357.5      110.0000000      0.8966090      124.3000000      0.8359104     110.7000000      0.8551390                    %because of cull effect
%  385.0      148.9473684      0.8966090      160.9210526      0.8359104      147.4285714      0.8551390
%  412.0      170.5263158      0.8966090      184.7368421      0.8359104      167.5000000      0.8551390
];

tW_gw150_BPA30(:,1)=tW_gw150_BPA30(:,1)+64;         % to put in dpf

data.tW_gw150A_BPA30 = tW_gw150_BPA30(:,[1 2]);
data.tW_gw150B_BPA30 = tW_gw150_BPA30(:,[1 4]);
data.tW_gw150C_BPA30 = tW_gw150_BPA30(:,[1 6]);

units.tW_gw150A_BPA30 = {'d', 'g'};  label.tW_gw150A_BPA30 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150A_BPA30 = {'gw150A-BPA30'};
units.tW_gw150B_BPA30 = {'d', 'g'};  label.tW_gw150B_BPA30 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150B_BPA30 = {'gw150B-BPA30'};
units.tW_gw150C_BPA30 = {'d', 'g'};  label.tW_gw150C_BPA30 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150C_BPA30 = {'gw150C-BPA30'};

temp.tW_gw150A_BPA30 = C2K(8.5); units.temp.tW_gw150A_BPA30 = 'K'; label.temp.tW_gw150A_BPA30 = 'temperature';
temp.tW_gw150B_BPA30 = C2K(8.5); units.temp.tW_gw150B_BPA30 = 'K'; label.temp.tW_gw150B_BPA30 = 'temperature';
temp.tW_gw150C_BPA30 = C2K(8.5); units.temp.tW_gw150C_BPA30 = 'K'; label.temp.tW_gw150C_BPA30 = 'temperature';



%---------------------------------------------Study gw124

% Our data for BPA100 (study gw124)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names:  date	weight.C55_C8_B15	surv.C55_C8_B15	weight.C61_C14_B15	surv.C61_C14_B15	weight.C67_C20_B15	surv.C67_C20_B15

tW_gw124_BPA100=[...
0   0.1118056	1	0.1045455	1	0.09874214	1
13	0.1415929	0.8055556	0.1472603	0.9675325	0.1527027	0.9496855
34	0.36375	0.6558505	0.3403101	0.9277709	0.3379845	0.8919344
59	1	0.6394543	0.9576271	0.9205788	0.93706897	0.8711917
76	1.7794118	0.6394543	1.6940678	0.9205788	1.675	0.8711917
97	3.2735294	0.6394543	3.0466102	0.9205788	3.05172414	0.8711917
118	5.5555556	0.6300505	5.3982301	0.9127773	5.3539823	0.8711917
139	9.0322581	0.6200497	8.4513274	0.9127773	8.125	0.8634821
160	14.7580645	0.6200497	14.3243243	0.896622	13.48214286	0.8634821
181	22.0338983	0.6200497	21.588785	0.8885443	20.63636364	0.8634821
202	32.7118644	0.6200497	32.3364486	0.8885443	30.77272727	0.8634821
223	46.2711864	0.6200497	47.9906542	0.8885443	43.62385321	0.8556322
244	61.1864407	0.6200497	60.9433962	0.8802401	55.97222222	0.8477824
286	100.3389831	0.6200497	105.5714286	0.871936	95.9047619	0.8242329
];

tW_gw124_BPA100(:,1)=tW_gw124_BPA100(:,1)+64;         % to put in dpf             !!!!!!!!!!!!!!!!!!!!!!!!!!!!    not sure....

data.tW_gw124A_BPA100 = tW_gw124_BPA100(:,[1 2]);
data.tW_gw124B_BPA100 = tW_gw124_BPA100(:,[1 4]);
data.tW_gw124C_BPA100 = tW_gw124_BPA100(:,[1 6]);

units.tW_gw124A_BPA100 = {'d', 'g'};  label.tW_gw124A_BPA100 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw124A_BPA100 = {'gw124A-BPA100'};
units.tW_gw124B_BPA100 = {'d', 'g'};  label.tW_gw124B_BPA100 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw124B_BPA100 = {'gw124B-BPA100'};
units.tW_gw124C_BPA100 = {'d', 'g'};  label.tW_gw124C_BPA100 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw124C_BPA100 = {'gw124C-BPA100'};

temp.tW_gw124A_BPA100 = C2K(8.5); units.temp.tW_gw124A_BPA100 = 'K'; label.temp.tW_gw124A_BPA100 = 'temperature';
temp.tW_gw124B_BPA100 = C2K(8.5); units.temp.tW_gw124B_BPA100 = 'K'; label.temp.tW_gw124B_BPA100 = 'temperature';
temp.tW_gw124C_BPA100 = C2K(8.5); units.temp.tW_gw124C_BPA100 = 'K'; label.temp.tW_gw124C_BPA100 = 'temperature';

% Our data for BPA100 end (study gw124)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names:  date	weight surv

tW_gw124_BPA100end=[...
375	238.1208	0.7589778
431	355.6376	0.7589778
488	510.6376	0.7589778
552	761.3087	0.7589778
552.5	769.3919	0.7589778
608	839.0541	0.7589778
664	1019.7297	0.7589778
720	1148.3	0.7589778
720.5	1171.7347	0.7589778
781	1360	0.7589778
837	1567.3469	0.7589778
894	1782.7551	0.7589778
894.5	1855.814	0.7589778
949	2014.3023	0.7589778
1005	2006.6279	0.7589778
1013.5	2006.6279	0.7589778
];


tW_gw124_BPA100end(:,1)=tW_gw124_BPA100end(:,1)+64;         % to put in dpf             !!!!!!!!!!!!!!!!!!!!!!!!!!!!    not sure....

data.tW_gw124_BPA100end = tW_gw124_BPA100end(:,[1 2]);

units.tW_gw124_BPA100end = {'d', 'g'};  label.tW_gw124_BPA100end = {'age since fertilization', 'wet weight'};  bibkey.tW_gw124_BPA100end = {'gw124-BPA100end'};

temp.tW_gw124_BPA100end = C2K(8.5); units.temp.tW_gw124_BPA100end = 'K'; label.temp.tW_gw124_BPA100end = 'temperature';



% Our data for BPA03 (study gw124)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names:  date	weight surv

tW_gw124_BPA03=[...
    0	0.1150735	1
13	0.172428	0.9411765
34	0.4146226	0.8753328
59	1.25	0.8753328
76	2.1742574	0.8753328
97	3.7572139	0.8709995
118	6.4736842	0.8493329
139	9.7326203	0.8359223
160	14.8918919	0.826982
181	22.4858757	0.826982
202	31.9318182	0.8223098
223	43.4571429	0.8176376
244	60.1436782	0.8129654
286	98.2369942	0.8082931
298.5	112.5333333	0.8082931
327	138.8	0.8082931
375	218.1208054	0.8029045
431	343.1208054	0.8029045
488	515.2013423	0.8029045
552	754.7315436	0.8029045
552.5	728.4459459	0.8029045
608	831.7567568	0.8029045
664	992.6351351	0.8029045
720	1181.8	0.8029045
781	1394.4	0.8029045
837	1587.083333	0.7707883
893	1843.369565	0.7386722
893.5	1960.131579	0.7386722
949	2264.210526	0.7386722
1005	2303.289474	0.7386722
1013.5	2303.289474	0.7386722
1019.5	2318.890454	0.7386722
];

tW_gw124_BPA03(:,1)=tW_gw124_BPA03(:,1)+64;         % to put in dpf             !!!!!!!!!!!!!!!!!!!!!!!!!!!!    not sure....

data.tW_gw124_BPA03 = tW_gw124_BPA03(:,[1 2]);

units.tW_gw124_BPA03 = {'d', 'g'};  label.tW_gw124_BPA03 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw124_BPA03 = {'gw124-BPA0.3'};

temp.tW_gw124_BPA03 = C2K(8.5); units.temp.tW_gw124_BPA03 = 'K'; label.temp.tW_gw124_BPA03 = 'temperature';



% Our data for BPA3 (study gw124)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names:  date	weight surv

tW_gw124_BPA3=[...
0	0.1114983	1
13	0.1673004	0.9512195
34	0.3970954	0.9186683
59	1.1731602	0.9186683
76	2.1038961	0.9186683
97	3.7884956	0.8987837
118	6.5366972	0.8908298
139	10.2790698	0.8785707
160	16.3785047	0.8744843
181	24.7087379	0.8744843
202	34.6097561	0.8702393
223	43.4068627	0.8659942
244	57.4137931	0.8617491
286	96.4179104	0.853259
298.5	105.9333333	0.853259
327	133.6577181	0.8475706
375	211.1486486	0.8418822
431	335.7094595	0.8418822
488	515.3040541	0.8418822
552	768.6148649	0.8418822
552.5	751.0273973	0.8418822
608	834.5205479	0.8418822
664	999.3150685	0.8418822
720	1192.755102	0.8418822
781	1355.204082	0.8418822
837	1508.673469	0.8418822
893	1690.106383	0.8075197
893.5	1827.794118	0.8075197
949	2160.735294	0.8075197
1005	2203.529412	0.8075197
1013.5	2203.529412	0.8075197
1019.5	2218.67711	0.8075197
];

tW_gw124_BPA3(:,1)=tW_gw124_BPA3(:,1)+64;         % to put in dpf             !!!!!!!!!!!!!!!!!!!!!!!!!!!!    not sure....

data.tW_gw124_BPA3 = tW_gw124_BPA3(:,[1 2]);

units.tW_gw124_BPA3 = {'d', 'g'};  label.tW_gw124_BPA3 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw124_BPA3 = {'gw124-BPA3'};

temp.tW_gw124_BPA3 = C2K(8.5); units.temp.tW_gw124_BPA3 = 'K'; label.temp.tW_gw124_BPA3 = 'temperature';


% Our data for BPA30 (study gw124)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names:  date	weight surv

tW_gw124_BPA30=[...
0	0.09906832	1
13	0.1511785	0.9534161
34	0.3774545	0.9245248
59	1.06654	0.9178009
76	1.94981	0.9178009
97	3.451923	0.9073317
118	5.885827	0.9073317
139	8.924303	0.8966152
160	14.2	0.893043
181	21.98347	0.893043
202	31.59751	0.8893528
223	40.46025	0.8819723
244	53.17597	0.8598307
286	87.48918	0.8524502
298.5	101.1	0.8450696
327	128.9667	0.8338021
375	209.094	0.8282434
431	330.6376	0.8282434
488	478.3893	0.8282434
552	714.9329	0.8282434
552.5	722.7703	0.8282434
608	815.3425	0.8170509
664	1006.301	0.8170509
720	1206.042	0.8058584
781	1390.312	0.8058584
837	1583.125	0.8058584
893	1751.489	0.7890697
893.5	1915	0.7890697
949	2169.73	0.7890697
1005	2177.432	0.7890697
1013.5	2177.432	0.7890697
1019.5	2192.46	0.7890697
];

tW_gw124_BPA30(:,1)=tW_gw124_BPA30(:,1)+64;         % to put in dpf             !!!!!!!!!!!!!!!!!!!!!!!!!!!!    not sure....

data.tW_gw124_BPA30 = tW_gw124_BPA30(:,[1 2]);

units.tW_gw124_BPA30 = {'d', 'g'};  label.tW_gw124_BPA30 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw124_BPA30 = {'gw124-BPA30'};

temp.tW_gw124_BPA30 = C2K(8.5); units.temp.tW_gw124_BPA30 = 'K'; label.temp.tW_gw124_BPA30 = 'temperature';

% Our data for BPA300 (study gw124)
%  days from first feeding (64dpf)
%  on 3 different tanks
%  weight : wet weight in g
%  surv : cumulated survival in the tank from first feeding
%  column names:  date	weight surv

tW_gw124_BPA300=[...
0	0.09329609	1
13	0.1023077	0.7877095
34	0.2214953	0.7210572
59	0.6319588	0.7210572
76	1.129167	0.7136236
97	1.969892	0.6913228
118	3.235294	0.6764557
139	5.180723	0.6605391
160	8.951613	0.6525808
181	13.16327	0.6525808
202	19.18367	0.6525808
223	28.26531	0.6525808
244	40.72917	0.6392628
286	73.4375	0.6392628
298.5	87	0.6392628
327	110.4167	0.6392628
375	180.8333	0.6392628
431	301.6667	0.6392628
488	480.4167	0.6392628
552	706.7021	0.6259448
552.5	706.7021	0.6259448
608	844.6739	0.6126269
664	1042.717	0.6126269
720	1304.286	0.6126269
781	1598.393	0.6126269
837	1881.964	0.6126269
893	2138.75	0.6126269
893.5	2138.75	0.6126269
949	2473.393	0.6126269
1005	2517.5	0.6126269
1013.5	2517.5	0.6126269
1019.5	2534.053	0.6126269
];

tW_gw124_BPA300(:,1)=tW_gw124_BPA300(:,1)+64;         % to put in dpf             !!!!!!!!!!!!!!!!!!!!!!!!!!!!    not sure....

data.tW_gw124_BPA300 = tW_gw124_BPA300(:,[1 2]);

units.tW_gw124_BPA300 = {'d', 'g'};  label.tW_gw124_BPA300 = {'age since fertilization', 'wet weight'};  bibkey.tW_gw124_BPA300 = {'gw124-BPA300'};

temp.tW_gw124_BPA300 = C2K(8.5); units.temp.tW_gw124_BPA300 = 'K'; label.temp.tW_gw124_BPA300 = 'temperature';


%% set weights for all real data
weights = setweights(data, []);

% %% set pseudodata and respective weights
% [data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
%txtData.comment = comment;


%% Facts
                                 
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
%
bibkey = 'Kooy2014'; type = 'Misc'; bib = ...
'note = {taken from from Salmo trutta}';  
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
