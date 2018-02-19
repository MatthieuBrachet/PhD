clc; clear all; close all;

N=[16 32 64 128];
dx=1./N;
e1=[5.0991e-1 1.4492e-1 1.1084e-2 3.8979e-4];
e2=[5.1131e-1 1.4528e-1 1.1094e-2 3.8982e-4];
e3=[5.1918e-1 1.4859e-1 1.1304e-2 9.9221e-4];

e1=[5.9120e-1 1.8406e-1 1.6794e-2 6.1900e-4];
e2=[7.3391e-1 4.5397e-1 1.5584e-1 2.8292e-2];
e3=[9.3716e-1 9.3614e-1 9.2169e-1 8.6938e-1];


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
set(hl1,'MarkerFaceColor','m');
set(hl1,'MarkerEdgeColor','k');

hold on;
[a1,b1]=polyfit(ldx,le1,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a1(2)+a1(1)*tt;
% %% courbes...
%hlf1=plot(tt,lsq1,'-r');
hlf1=plot(tt,lsq1,'-k');
set(hlf1,'LineWidth',2.0);


%% curve 2
hl2=plot(ldx,le2,'ok');
set(hl2,'LineWidth',2.0);
set(hl2,'MarkerSize',10);
set(hl2,'MarkerFaceColor','b');
set(hl2,'MarkerEdgeColor','k');

hold on;
[a2,b1]=polyfit(ldx,le2,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a2(2)+a2(1)*tt;
% %% courbes...
%hlf2=plot(tt,lsq1,'-b');
hlf2=plot(tt,lsq1,'-k');
set(hlf2,'LineWidth',2.0);

%% curve 3
hl3=plot(ldx,le3,'ok');
set(hl3,'LineWidth',2.0);
set(hl3,'MarkerSize',10);
set(hl3,'MarkerFaceColor','r');
set(hl3,'MarkerEdgeColor','k');

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
set(ya,'String','Log_{10}(error)');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

%% legend
legend([hl1,hl2,hl3],{['sans filtre           - rate : ' num2str(a1(1))],['Filtre d''ordre 10 - rate : ' num2str(a2(1))],['Filtre d''ordre 8   - rate : ' num2str(a3(1))]},'Location','SouthEast')
