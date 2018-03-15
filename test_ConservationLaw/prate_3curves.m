clc; clear all; close all;

radius=1;
N=[16 32 64 128];
dx=2*pi*radius./(4*N);
% e1=[3.1446e-5 2.2000e-6 1.4092e-7 8.7856e-9];
% e2=[1.5759e-5 1.1752e-6 7.9823e-8 5.0291e-9];
% e3=[1.4251e-5 1.0776e-6 7.7308e-8 4.5510e-9];

% e1=[6.9356e-6 4.4315e-7 2.8338e-8 1.7687e-9];
% e2=[4.3733e-6 2.9918e-7 1.9272e-8 1.2084e-9];
% e3=[6.2809e-6 4.3611e-7 2.7803e-8 1.7354e-9];

e1=[1.3615 0.2034 2.8087e-2 4.0114e-3];
e2=[0.2499 5.0362e-2 1.0272e-2 2.1867e-3];
e3=[.0570 1.6202e-2 4.8284e-3 1.4308e-3];

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

ht=text('Position',[-1.95,-7.5,0],'String','N=128');
set(ht,'FontSize',12);
ht=text('Position',[-1.65,-6.5,0],'String','N=64');
set(ht,'FontSize',12);
ht=text('Position',[-1.35,-5.3,0],'String','N=32');
set(ht,'FontSize',12);
ht=text('Position',[-1.05,-4.2,0],'String','N=16');
set(ht,'FontSize',12);

