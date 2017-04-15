clc; clear all; close all;

radius=6371220;
N=[32 64 128];
dx=radius./(4*N);
e1=[1.1422*10^-6 7.1216*10^-8 4.4469*10^-9];
e2=[1.3885*10^-6 8.6513*10^-8 5.4018*10^-9];
e3=[2.4469*10^-6 1.5229*10^-7 9.5186*10^-9];
ldx=log(dx);
le1=log(e1);
le2=log(e2);
le3=log(e3);



%% plot 1
figure(11); 
%% curve 1
hl1=loglog(dx,e1,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');
set(hl1,'MarkerFaceColor','r');

hold on;
[a11,b1]=polyfit(ldx,le1,1);
ltt=linspace(min(ldx),max(ldx),1000);
lsq1=a11(2)+a11(1)*ltt;

loglog(dx,exp(a11(2).*dx.^a11(1)),'x')