clc; clear all; close all;

radius=6371220;
N=[16 32 64 128 256 512];
dx=2*pi*radius./(4*N);
e1=[1.0377e-4 6.3236e-6 3.9444e-7 2.4726e-8 1.5500e-9 9.7139e-11];
e2=[1.2588e-4 7.4682e-6 4.6278e-7 2.8931e-8 1.8111e-9 1.1342e-10];
e3=[3.6670e-4 2.2222e-5 1.3713e-6  8.5415e-8 5.3339e-9 3.3330e-10];

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
hlf3=plot(tt,lsq1,'-k');
set(hlf3,'LineWidth',2.0);

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

%% legend
legend([hlf1,hlf2,hlf3],{['norme 1 - pente = ' num2str(a1(1))],['norme 2 - pente = ' num2str(a2(1))],['norme \infty - pente = ' num2str(a3(1))]},'Location','SouthEast')
