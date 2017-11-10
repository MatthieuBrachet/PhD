clc; clear all; close all;

radius=6371220;
N=[16 32 64 128 256];
dx=2*pi*radius./(4*N);
%e=[1.0556*10^-3 7.9898*10^-5 5.8219*10^-6 3.9327*10^-7 2.5537*10^-8];
%e=[1.8636*10^-4 1.3545*10^-5 9.7564*10^-7 6.5592*10^-8 4.2563*10^-9];
%e=[5.0425*10^-3 3.9160*10^-4 3.1917*10^-5 2.2341*10^-6 1.4700*10^-7];
%e=[1.8475*10^-5 2.0130*10^-6 2.1919*10^-7 2.4346*10^-8 2.9159*10^-9];

%e=[1.7775*10^-3 1.4246*10^-4 1.0058*10^-5 6.6737*10^-7 4.2949*10^-8];
%e=[3.1368*10^-4 2.3910*10^-5 1.6716*10^-6 1.1088*10^-7 7.1437*10^-9];
%e=[1.1437*10^-2 9.7399*10^-4 6.6478*10^-5 4.2363*10^-6 2.6521*10^-7];
%e=[3.1349*10^-5 3.4020*10^-6 3.9216*10^-7 4.7945*10^-8 3.8473*10^-9]

%e=[3.1553*10^-4 2.6911*10^-5 1.9635*10^-6 1.3267*10^-7 8.6173*10^-9];
%e=[5.0919*10^-5 4.4192*10^-6 3.2504*10^-7 2.2033*10^-8 1.4340*10^-9];
%e=[1.3447*10^-3 1.4283*10^-4 1.1211*10^-5 7.7742*10^-7 5.1040*10^-8];
e=[4.3631*10^-6 2.2195*10^-7 3.3456*10^-8 4.5292*10^-9 5.8348*10^-10];

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


