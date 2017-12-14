clc; clear all; close all;

radius=6371220;
N=[8 16 32 64 128 256];
dx=2*pi*radius./(4*N);
% e1=[1.9936e-2 1.1325e-3 6.8262e-5 4.2292e-6 2.6413e-7 1.6559e-8];
% e2=[1.8717e-2 1.0216e-3 6.1471e-5 3.8166e-6 2.3881e-7 1.4975e-8];
% e3=[2.5623e-2 1.4497e-3 8.9084e-5 5.5817e-6 3.5265e-7 2.2194e-8];

% e1=[1.9936e-2 1.1325e-3 6.8262e-5 4.2292e-6 2.6413e-7 1.6559e-8];
% e2=[1.8717e-2 1.0216e-3 6.1471e-5 3.8166e-6 2.3881e-7 1.4975e-8];
% e3=[2.5623e-2 1.4497e-3 8.9084e-5 5.5817e-6 3.5265e-7 2.2194e-8];

e1=[3.4770e-4 1.6686e-5 1.0089e-6 6.3688e-8 4.0283e-9 2.5396e-10];
e2=[1.4286e-4 7.0948e-6 4.1258e-7 2.5581e-8 1.6143e-9 1.0283e-10];
e3=[1.1470e-4 1.0256e-5 6.2778e-7 3.7403e-8 2.2653e-9 1.9459e-10];


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
hlf1=plot(tt,lsq1,'-b');
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
hlf2=plot(tt,lsq1,'-g');
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
hlf3=plot(tt,lsq1,'-r');
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
legend([hlf1,hlf2,hlf3],{['norm 1 - slope = ' num2str(a1(1))],['norm 2 - slope = ' num2str(a2(1))],['norm \infty - slope = ' num2str(a3(1))]},'Location','SouthEast')
