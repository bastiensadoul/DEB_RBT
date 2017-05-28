function [data, auxData, metaData, txtData, weights] = mydata_reconstruction

%% set metadata
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Actinopterygii'; 
metaData.order      = 'Salmoniformes'; 
metaData.family     = 'Salmonidae';
metaData.species    = 'Oncorhynchus_mykiss'; 
metaData.species_en = 'Rainbow trout'; 

metaData.author   = {'Bastien Sadoul'};    


%% set data

% BPA0 gw150, time in dpf, weight in g and survival in 3 tanks
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
tW_gw150(:,1)=tW_gw150(:,1)+64;         % to put in dpf
% 
% temp.tW_gw150mean=8.5;      % °C
% data.tW_gw150mean (:,1) = tW_gw150(:,1);
% data.tW_gw150mean (:,2) = 1/3*(tW_gw150(:,1)+tW_gw150(:,2)+tW_gw150(:,4)+tW_gw150(:,6));
% 
% units.tW_gw150mean = {'d', 'g'};  label.tW_gw150mean = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150mean = {'gw150A'};
% 
% auxData.t0.tW_gw150mean  = 'dpf';

data.tW_gw150A = tW_gw150(:,[1 2]);
data.tW_gw150B = tW_gw150(:,[1 4]);
data.tW_gw150C = tW_gw150(:,[1 6]);

temp.tW_gw150 = C2K(8.5);

units.tW_gw150A = {'d', 'g'};  label.tW_gw150A = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150A = {'gw150A'};
units.tW_gw150B = {'d', 'g'};  label.tW_gw150B = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150B = {'gw150B'};
units.tW_gw150C = {'d', 'g'};  label.tW_gw150C = {'age since fertilization', 'wet weight'};  bibkey.tW_gw150C = {'gw150C'};

auxData.t0.tW_gw150A  = 'dpf';
auxData.t0.tW_gw150B  = 'dpf';
auxData.t0.tW_gw150C  = 'dpf';
  
%% set weights for all real data
weights = setweights(data, []);
weights.tW_gw150A(1:10) = weights.tW_gw150A(1:10) * 100;
weights.tW_gw150B(1:10) = weights.tW_gw150B(1:10) * 100;
weights.tW_gw150C(1:10) = weights.tW_gw150C(1:10) * 100;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;

%% References

bibkey.tW_gw150A = {'gw150A-BPA0'};
bibkey.tW_gw150B = {'gw150B-BPA0'};
bibkey.tW_gw150C = {'gw150C-BPA0'};

