% calculate the residues (difference between model and predictions)

clear all; clc; close all;

%% about the parameters
% these are new parameters assuming k_J = 0.002 ; but the control
% maintenance is now MUCH MUCH higher

% anyways I got f values for all of the tanks and in this script you will
% see the prediced vs observed  values which are loglogted ....


%% controls :
load('results_Oncorhynchus_mykiss_BPA0.mat');
[data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss_BPA0;
[prdData, info] = predict_Oncorhynchus_mykiss_BPA0(par, data, auxData);
% just to get the names shorter :
d = data;
pD = prdData;



%% exposed
load('results_Oncorhynchus_mykiss_BPA03to300.mat');
[data, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss_BPA03to300;
[prdData, info] = predict_Oncorhynchus_mykiss_BPA03to300 (par, data, auxData);
% just to get the names shorter :
dC = data;
pDc = prdData;



%% makes loglogs of residues 
% tell me if you think of better way I can do this ...
% one example has it in log-log, the others are not in log scale


figure(1)
subplot(131)
loglog(d.tW_gw150A(:,2),pD.tW_gw150A,'.b'); hold on
% loglog([-2 6],[-2 6],'k')
% xlim([-2 6]); ylim([-2 6]); 
title('tW_gw150A') ; xlabel('observed'); ylabel('predicted')


subplot(132)
loglog(d.tW_gw150B(:,2),pD.tW_gw150B,'.b'); hold on
loglog([1e-2 140],[1e-2 140],'k')
% xlim([1e-2 140]); ylim([0 140]); 
title('tW_gw150B') ; xlabel('observed'); ylabel('predicted')
subplot(133)
loglog(d.tW_gw150C(:,2),pD.tW_gw150C,'.b'); hold on
% loglog([1e-2 140],[1e-2 140],'k')
% xlim([0 140]); ylim([0 140]); 
title('tW_gw150C') ; xlabel('observed'); ylabel('predicted')


figure(2)
subplot(131)
loglog(d.tW_gw124iniA(:,2),pD.tW_gw124iniA,'.b'); hold on
% loglog([1e-2 140],[1e-2 140],'k')
% xlim([0 60]); ylim([0 60]); 
title('tW_gw124iniA') ; xlabel('observed'); ylabel('predicted')
% the tox data :
loglog(dC.tW_gw124A_BPA100(:,2),pDc.tW_gw124A_BPA100,'.r')

subplot(132)
loglog(d.tW_gw124iniB(:,2),pD.tW_gw124iniB,'.b'); hold on
% loglog([1e-2 140],[1e-2 140],'k')
% xlim([0 60]); ylim([0 60]); 
title('tW_gw124iniB ') ; xlabel('observed'); ylabel('predicted')
loglog(dC.tW_gw124B_BPA100(:,2),pDc.tW_gw124B_BPA100,'.r')


subplot(133)
loglog(d.tW_gw124iniC(:,2),pD.tW_gw124iniC,'.b'); hold on
% loglog([1e-2 140],[1e-2 140],'k')
% xlim([0 60]); ylim([0 60]); 
title('tW_gw124iniC') ; xlabel('observed'); ylabel('predicted')
loglog(dC.tW_gw124C_BPA100(:,2),pDc.tW_gw124C_BPA100,'.r')

figure(3)
loglog(d.tW_gw124fin(:,2),pD.tW_gw124fin,'.b'); hold on
% loglog([1e-4 2500],[1e-4 2500],'k')
% xlim([1e-4 2500]); ylim([1e-4 2500]); 
title('tW_gw124fin') ; xlabel('observed'); ylabel('predicted')

loglog(dC.tW_gw124_BPA100end(:,2),pDc.tW_gw124_BPA100end,'or'); hold on


%% exposed

names=fieldnames(pDc);
for i=1:length(names) 
    name=names{i};
    field=getfield(pDc,name);
    dCfield=getfield(dC,name);
    field(:,2)=dCfield(:,2);
    field(:,3)=dCfield(:,1);
    pDc = setfield(pDc,name,field);
end
    
names=fieldnames(pD);
for i=1:length(names) 
    name=names{i};
    field=getfield(pD,name);
    dfield=getfield(d,name);
    field(:,2)=dfield(:,2);
    field(:,3)=dfield(:,1);
    pD = setfield(pD,name,field);
end
    

save('prdDataControls.mat','pD');
save('prdDataExposed.mat','pDc');


