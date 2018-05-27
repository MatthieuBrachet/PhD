clc; clear all; close all;

radius=6371220;
N=[40 50 60 80 100 150];
dx=2*pi*radius./(4*N);
% e1=[3.7638e-2 2.1323e-2 1.3546e-2 6.6905e-3 3.9119e-3 1.4922e-3];
% e2=[2.0633e-2 1.2496e-2 8.4518e-3 4.6505e-3 2.9448e-3 1.3019e-3];
% e3=[1.4639e-2 1.0056e-2 7.3167e-3 4.6060e-3 3.1809e-3 1.6341e-3];

e1=[4.3043e-2 2.4403e-2 1.5367e-2 7.5508e-3 4.3709e-3 1.6538e-3];
e2=[2.4784e-2 1.4917e-2 9.9131e-3 5.3960e-3 3.3958e-3 1.4917e-3];
e3=[2.0921e-2 1.3748e-2 1.0476e-2 6.2646e-3 4.4360e-3 2.2885e-3]; 

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
%set(hl1,'MarkerFaceColor','m');
%set(hl1,'MarkerEdgeColor','k');

hold on;
[a1,b1]=polyfit(ldx,le1,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a1(2)+a1(1)*tt;
%a1=[3.8523 -35.1389]; lsq1=a1(2)+a1(1)*tt;
% %% courbes...
%hlf1=plot(tt,lsq1,'-r');
hlf1=plot(tt,lsq1,':k');
set(hlf1,'LineWidth',2.0);


%% curve 2
hl2=plot(ldx,le2,'ok');
set(hl2,'LineWidth',2.0);
set(hl2,'MarkerSize',10);
%set(hl2,'MarkerFaceColor','b');
%set(hl2,'MarkerEdgeColor','k');

hold on;
[a2,b1]=polyfit(ldx,le2,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a2(2)+a2(1)*tt;
%a2=[3.8805 -35.1747]; lsq1=a2(2)+a2(1)*tt;
% %% courbes...
%hlf2=plot(tt,lsq1,'-b');
hlf2=plot(tt,lsq1,'--k');
set(hlf2,'LineWidth',2.0);

%% curve 3
hl3=plot(ldx,le3,'ok');
set(hl3,'LineWidth',2.0);
set(hl3,'MarkerSize',10);
%set(hl3,'MarkerFaceColor','r');
%set(hl3,'MarkerEdgeColor','k');

hold on;
[a3,b1]=polyfit(ldx,le3,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a3(2)+a3(1)*tt;
%a3=[3.9806 -35.3017]; lsq1=a3(2)+a3(1)*tt;
% %% courbes...
%hlf3=plot(tt,lsq1,'-g');
hlf3=plot(tt,lsq1,'-k');
set(hlf3,'LineWidth',2.0);

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
legend([hlf1,hlf2,hlf3],{['norme 1 - pente = ' num2str(a1(1))],['norme 2 - pente = ' num2str(a2(1))],['norme \infty - pente = ' num2str(a3(1))]},'Location','NorthWest')

%% texte
ht=text('Position',[4.5,-8.3,0],'String','N=256');
set(ht,'FontSize',12);
ht=text('Position',[4.8,-9,0],'String','N=128');
set(ht,'FontSize',12);
ht=text('Position',[5.1,-7.9,0],'String','N=64');
set(ht,'FontSize',12);
ht=text('Position',[5.4,-6.5,0],'String','N=32');
set(ht,'FontSize',12);
ht=text('Position',[5.7,-5.5,0],'String','N=16');
set(ht,'FontSize',12);
ht=text('Position',[6,-4.2,0],'String','N=8');
set(ht,'FontSize',12);