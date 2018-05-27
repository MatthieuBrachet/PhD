clc; clear all; close all;

radius=6371220;
N=[32 64 128];
dx=2*pi*radius./(4*N);

e1=[1.1422e-6 7.1216e-8 4.4469e-9];
e2=[1.0446e-6 6.5124e-8 4.0672e-9];
e3=[2.4469e-6 1.5229e-7 9.5186e-9];

e1=[7.5712e-7 4.7213e-8 2.9487e-9];
e2=[1.0446e-6 6.5124e-8 4.0672e-9];
e3=[2.7809e-6 1.7387e-7 1.0858e-8];


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
ht=text('Position',[4.95, -8.4,0],'String','N=128');
set(ht,'FontSize',12);
ht=text('Position',[5.18,-7.6,0],'String','N=64');
set(ht,'FontSize',12);
ht=text('Position',[5.45,-6.5,0],'String','N=32');
set(ht,'FontSize',12);