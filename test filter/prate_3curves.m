clc; clear all; close all;

radius=6371220;
N=[8 16 32 64 128 256];
dx=2*pi*radius./(4*N);
e1=[6.0168e-03   1.5281e-03   3.8365e-04   9.6037e-05   2.4020e-05   6.0062e-06];
e2=[6.0690e-03   1.5225e-03   3.8031e-04   9.4973e-05   2.3726e-05   5.9290e-06];
e3=[2.4651e-6 8.4471e-9 6.1885e-10 3.9985e-11 2.5708e-12 1.6163e-13];

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

hold on;
[a1,b1]=polyfit(ldx,le1,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a1(2)+a1(1)*tt;
% %% courbes...
%hlf1=plot(tt,lsq1,'-r');
hlf1=plot(tt,lsq1,':k');
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
%hlf2=plot(tt,lsq1,'-b');
hlf2=plot(tt,lsq1,'--k');
set(hlf2,'LineWidth',2.0);

%% curve 3
hl1=plot(ldx,le3,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');

hold on;
[a3,b1]=polyfit(ldx,le3,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a3(2)+a3(1)*tt;
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
legend([hlf1,hlf2,hlf3],{['norm 1 - slope = ' num2str(a1(1))],['norm 2 - slope = ' num2str(a2(1))],['norm \infty - slope = ' num2str(a3(1))]},'Location','SouthEast')
