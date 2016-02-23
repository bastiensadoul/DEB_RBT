function [data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss_BPA3and30

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
metaData.email    = {'bastien.sadoul@hotmail.fr'};                 
metaData.address  = {'University of Calgary'}; 

%% set data

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
