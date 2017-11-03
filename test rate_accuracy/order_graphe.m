clc; clear all; close all;

radius=6371220;
N=[40 50 60 80 100 150];
%N=[32 64 128];
% N=[16 32 64 128];
dx=2*pi*radius./(4*N);

%% Nair-Machenhauer
e1=[1.989*10^-3 7.638*10^-4 3.023*10^-4 5.2979*10^-5 1.5036*10^-5 1.9244*10^-6];
e2=[7.255*10^-3 3.161*10^-3 1.313*10^-3 2.391*10^-4 6.4568*10^-5 9.2082*10^-6];
e3=[4.039*10^-2 1.918*10^-2 7.556*10^-3 1.561*10^-3 4.329*10^-4 7.6848*10^-5];
%% Nair-Jablonowski
% e1=[2.2199*10^-3 9.9676*10^-4 4.4566*10^-4 1.3189*10^-4 5.3380*10^-5 1.0401*10^-5];
% e2=[8.1592*10^-3 3.9476*10^-3 1.9545*10^-3 5.8682*10^-4 2.3789*10^-4 4.7173*10^-5];
% e3=[4.9298*10^-2 2.8948*10^-2 1.6474*10^-2 5.6752*10^-3 2.3177*10^-3 4.9192*10^-4];
%% Williamson 2 (alpha=0)
% e1=[1.1422*10^-6 7.1216*10^-8 4.4469*10^-9];
% e2=[1.3885*10^-6 8.6513*10^-8 5.4018*10^-9];
% e3=[2.4469*10^-6 1.5229*10^-7 9.5186*10^-9];
%% Williamson 2 (alpha=pi/4)
% e1=[7.5712*10^-7 4.7213*10^-8 2.9487*10^-9];
% e2=[1.0446*10^-6 6.5124*10^-8 4.0672*10^-9];
% e3=[2.7809*10^-6 1.7387*10^-7 1.0858*10^-8];

% e1=[4.4750*10^-2 5.7696*10^-3 3.4225*10^-4 1.6969*10^-5];

ldx=log10(dx);
le1=log10(e1);
le2=log10(e2);
le3=log10(e3);



%% plot 1
figure(11); 
%% curve 1
hl1=plot(ldx,le1,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');
%set(hl1,'MarkerFaceColor','r');

hold on;
[a11,b1]=polyfit(ldx,le1,1);
tt=linspace(min(ldx),max(ldx),1000);
lsq1=a11(2)+a11(1)*tt;


%% curve 2
hl2=plot(ldx,le2,'+k');
set(hl2,'LineWidth',2.0);
set(hl2,'MarkerSize',10);
set(hl2,'MarkerEdgeColor','k');
set(hl2,'MarkerFaceColor','green');

hold on;
[a12,b1]=polyfit(ldx,le2,1);
tt=linspace(min(ldx),max(ldx),1000);
lsq2=a12(2)+a12(1)*tt;


%% curve 3
hl3=plot(ldx,le3,'xk');
set(hl3,'LineWidth',2.0);
set(hl3,'MarkerSize',10);
set(hl3,'MarkerEdgeColor','k');
set(hl3,'MarkerFaceColor','magenta');

hold on;
[a13,b1]=polyfit(ldx,le3,1);
tt=linspace(min(ldx),max(ldx),1000);
lsq3=a13(2)+a13(1)*tt;


%% legend
legend([hl1,hl2,hl3],{['norm 1   : slope = ' num2str(a11(1))],['norm 2   : slope = ' num2str(a12(1))],['norm inf : slope = ' num2str(a13(1))]},'Location','NorthWest')
% legend(['norm inf : slope = ' num2str(a11(1))],'Location','NorthWest')
% xlabel
ha=gca;
xa=get(gca,'Xlabel');
set(ha,'Xgrid','on');
set(xa,'String','Log10(\Delta)');
set (xa,'FontName','Calibri');
set(xa,'FontSize',12);

% ylabel
ha=gca;
ya=get(gca,'Ylabel');
set(ha,'Ygrid','on');
set(ya,'String','Log10(error_{max})');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

% title
ha=gca;
xt=get(gca,'Title');
title('Convergence rate');
set (xt,'FontSize',12);


% %% courbes...
hlf1=plot(tt,lsq1,'--k');
set(hlf1,'LineWidth',2.0);

hlf1=plot(tt,lsq2,'-k');
set(hlf1,'LineWidth',2.0);

hlf1=plot(tt,lsq3,'-.k');
set(hlf1,'LineWidth',2.0);

%% texte
% ht=text('Position',[4.1,-8.65,0],'String','N=128');
% set(ht,'FontSize',12);
% ht=text('Position',[4.1,-8.8,0],'String','\Deltat=2.5min');
% set(ht,'FontSize',12);
% 
% ht=text('Position',[4.38,-7.65,0],'String','N=64');
% set(ht,'FontSize',12);
% ht=text('Position',[4.38,-7.8,0],'String','\Deltat=5min');
% set(ht,'FontSize',12);
% 
% ht=text('Position',[4.6,-6.6,0],'String','N=32');
% set(ht,'FontSize',12);
% ht=text('Position',[4.6,-6.75,0],'String','\Deltat=10min');
% set(ht,'FontSize',12);
