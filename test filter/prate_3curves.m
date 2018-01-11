clc; clear all; close all;

radius=6371220;
N=[8 16 32 64 128 256];
dx=2*pi*radius./(4*N);
% e1=[2.0963e-2 1.5863e-2 4.9861e-4 1.2394e-6 1.7658e-9 9.5777e-12];
% e2=[2.6194e-2 2.3535e-2 7.8589e-4 1.9919e-6 3.6229e-9 1.1097e-10];
% e3=[4.4752e-2 6.6690e-2 2.9727e-3 8.5586e-6 7.6882e-8 4.5050e-9];

% e1=[8.5678e-5 2.4618e-7 1.0687e-9 2.2554e-11 6.6050e-13 2.0170e-14];
% e2=[1.3701e-4 3.8854e-7 4.7661e-9 1.8082e-10 7.4291e-12 3.1657e-13];
% e3=[4.5859e-4 1.5631e-6 8.6449e-8 4.4339e-9 2.5436e-10 1.5149e-11];

% e1=[1.3632e-7 5.7406e-10 1.1462e-11 3.2764e-13 1.0064e-14 3.3383e-16];
% e2=[3.1191e-7 1.1452e-9 3.8740e-11 1.7535e-12 7.8953e-14 3.5222e-15];
% e3=[2.4651e-6 8.4471e-9 6.1885e-10 3.9985e-11 2.5708e-12 1.6163e-13];

% e1=[1.1315e-3 2.8629e-4 7.1814e-5 1.7973e-5 4.4950e-6 1.1239e-6];
% e2=[1.1413e-3 2.8525e-4 7.1189e-5 1.7773e-5 4.4398e-6 1.1095e-6];
% e3=[1.6770e-3 4.3815e-4 1.1085e-4 2.7766e-5 6.9449e-6 1.7366e-6];

e1=[2.7388e-5 1.7510e-6 1.1017e-7 6.9021e-9 4.3178e-10 2.6997e-11];
e2=[3.2726e-5 2.1184e-6 1.3355e-7 8.3682e-9 5.2347e-10 3.2729e-11];
e3=[7.0509e-5 4.7655e-6 3.0421e-7 1.9135e-8 1.2152e-9 7.8438e-11];

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
set(ya,'String','Log_{10}(error)');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

%% legend
legend([hlf1,hlf2,hlf3],{['norm 1 - slope = ' num2str(a1(1))],['norm 2 - slope = ' num2str(a2(1))],['norm \infty - slope = ' num2str(a3(1))]},'Location','SouthEast')
