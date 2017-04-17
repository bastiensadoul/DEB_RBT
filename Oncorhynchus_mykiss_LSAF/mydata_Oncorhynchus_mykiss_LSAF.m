function [data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss_LSAF

%% set metadata
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Salmoniformes'; 
metaData.family     = 'Salmonidae';
metaData.species    = 'Oncorhynchus_mykiss_LSAF'; 
metaData.species_en = 'Rainbow trout, large size at age family'; 
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
  


% McKenPed2007 for large size-at-age family (LSAF)
% Colums of wO2:
%  1 W g, wet weight
%  2 mean O2 uptake per day, mmol/d
WwJO = [...
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
data.WwJO = WwJO; units.WwJO = {'g', 'mmol/d'}; label.WwJO  = {'wet weight', 'O2 uptake'}; bibkey.WwJO = {'McKenPed2007'};
comment.WwJO = 'LSAF - large size at age family. Water flow of 0.55BL/s, probably at the 182g. At 239.9g food was withdrawn for 2 days.';
temp.WwJO = C2K(14); units.temp.WwJO = 'K' ;  label.temp.WwJO = 'mean temperature' ; 

tWw = [ ...
0    	181.5
21.0	242.0
42	    319.4
63.0	371.9
85   	443.1
];
data.tWw = tWw; units.tWw = {'d', 'g'}; label.tWw  = {'time', 'wet weight'}; bibkey.tWw = {'McKenPed2007'};
comment.tWw = 'LSAF - large size at age family. Circles fig. 1';
temp.tWw = C2K(14); units.temp.tWw = 'K' ;  label.temp.tWw = 'mean temperature' ; 

tL = [...
 0 1.3
 21 1.4
 42 1.5
];
tL(:,2) = (100 .* tWw(1:3,2) ./ tL(:,2)).^(1/3);
data.tL = tL; units.tL = {'d', 'cm'}; label.tL  = {'time', 'fork length'}; bibkey.tL = 'McKenPed2007';
comment.tL = 'LSAF - large size at age family. Value computed from CF pp 282';
temp.tL = C2K(14); units.temp.tL = 'K' ;  label.temp.tL = 'mean temperature' ; 


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


