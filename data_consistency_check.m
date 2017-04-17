clear all; clc;

[d, auxData, txtData, weights] = mydata_BPA;

figure(1)
hold on;
plot(d.tW_gw150A(:,1), d.tW_gw150A(:,2), 'bo','markerfacecolor','b','linestyle','none')
plot(d.tW_gw150A_BPA3(:,1), d.tW_gw150A_BPA3(:,2), 'md','markerfacecolor','none','linestyle','none')
plot(d.tW_gw150A_BPA30(:,1), d.tW_gw150A_BPA30(:,2), 'rs','markerfacecolor','none','linestyle','none')