clear all; clc;

[d0, auxData, metaData, txtData0, weights] = mydata_Oncorhynchus_mykiss_BPA0;
[d, auxData, metaData, txtData, weights] = mydata_Oncorhynchus_mykiss_BPA3and30;


figure(1)
hold on;
plot(d0.tW_gw150A(:,1), d0.tW_gw150A(:,2), 'bo','markerfacecolor','b','linestyle','none')
plot(d.tW_gw150A_BPA3(:,1), d.tW_gw150A_BPA3(:,2), 'md','markerfacecolor','none','linestyle','none')
plot(d.tW_gw150A_BPA30(:,1), d.tW_gw150A_BPA30(:,2), 'rs','markerfacecolor','none','linestyle','none')