clc; clear all; close all;

radius=6371220;
N=[16 32 64 128 256];
dx=2*pi*radius./(4*N);
e=[1.8636*10^-4 1.3545*10^-5 9.7564*10^-7 6.5592*10^-8 4.2563*10^-9];
%e=[5.0425*10^-3 3.9159*10^-4 3.1916*10^-5 2.2341*10^-6 1.4700*10^-7];

ldx=log10(dx);
leh=log10(e);



%% plot 1
figure(11); 
%% curve 1
hl1=plot(ldx,leh,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');

hold on;
[a11,b1]=polyfit(ldx,leh,1);
tt=linspace(min(ldx),max(ldx),1000);
lsq1=a11(2)+a11(1)*tt;

% %% courbes...
hlf1=plot(tt,lsq1,'-k');
set(hlf1,'LineWidth',2.0);

%% legend
legend(hl1,{['norm on height   : slope = ' num2str(a11(1))]},'Location','SouthEast')

ha=gca;
xa=get(gca,'Xlabel');
set(ha,'Xgrid','on');
set(xa,'String','Log_{10}(dx)');
set (xa,'FontName','Calibri');
set(xa,'FontSize',12);

% ylabel
ha=gca;
ya=get(gca,'Ylabel');
set(ha,'Ygrid','on');
set(ya,'String','Log_{10}(error)');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

% title
ha=gca;
xt=get(gca,'Title');
title('Convergence rate');
set (xt,'FontSize',12);


