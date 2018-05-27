clc; clear all; close all;

radius=6371220;
N=[16 32 64 128 256 512];
dx=2*pi*radius./(4*N);
%e1=[1.8475*10^-5 2.0130*10^-6 2.1919*10^-7 2.4346*10^-8 2.9159*10^-9 3.5643*10^-10];
e1=[5.0425*10^-3 3.9160*10^-4 3.1917*10^-5 2.2341*10^-6 1.4700*10^-7 9.4165*10^-9];
e2=[1.8636*10^-4 1.3545*10^-5 9.7564*10^-7 6.5592*10^-8 4.2563*10^-9 2.7115*10^-10];
%e1=[3.1349*10^-5 3.4020*10^-6 3.9216*10^-7 4.7945*10^-8 3.8473*10^-9 7.2161*10^-10];
% e1=[1.1437*10^-2 9.7399*10^-4 6.6478*10^-5 4.2363*10^-6 2.6521*10^-7 1.6554*10^-8];
% e2=[3.1368*10^-4 2.3910*10^-5 1.6716*10^-6 1.1088*10^-7 7.1437*10^-9 4.5361*10^-10];

ldx=log10(dx);
le1=log10(e1);
le2=log10(e2);

%% plot 1
figure(11); 
%% curve 1
hl1=plot(ldx,le1,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');

hold on;
[a1,b1]=polyfit(ldx,le1,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a1(2)+a1(1)*tt;
% %% courbes...
hlf1=plot(tt,lsq1,'k--');
set(hlf1,'LineWidth',2.0);


%% curve 2
hl1=plot(ldx,le2,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');

hold on;
[a2,b1]=polyfit(ldx,le2,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a2(2)+a2(1)*tt;
% %% courbes...
hlf2=plot(tt,lsq1,'-k');
set(hlf2,'LineWidth',2.0);

ha=gca;
xa=get(gca,'Xlabel');
set(ha,'Xgrid','on');
set(xa,'String','Log_{10}(\Delta)');
set (xa,'FontName','Calibri');
set(xa,'FontSize',12);

% ylabel
ha=gca;
ya=get(gca,'Ylabel');
set(ha,'Ygrid','on');
set(ya,'String','Log_{10}(erreur)');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

%% legend
legend([hlf1,hlf2],{['pente = ' num2str(a1(1))],['pente = ' num2str(a2(1))]},'Location','SouthEast')
