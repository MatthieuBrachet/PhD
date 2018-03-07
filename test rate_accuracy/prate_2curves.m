clc; clear all; close all;

radius=6371220;
N=[32 64 128];
dx=2*pi*radius./(4*N);
e1=[1.2725*10^-1 7.0431*10^-3 2.7248*10^-4];
e2=[1.5352*10^-3 1.5060*10^-4 6.6822*10^-6];

%  e1=[5.1609*10^-4 2.2834*10^-5 9.9035*10^-7];
%  e2=[1.5944*10^-3 1.5444*10^-4 6.8467*10^-6];


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

% ylabel
ha=gca;
ya=get(gca,'Ylabel');
set(ha,'Ygrid','on');
set(ya,'String','Log_{10}(error)');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

%% legend
legend([hlf1,hlf2],{['error on height - rate : ' num2str(a1(1))],['error on velocity - rate : ' num2str(a2(1))]},'Location','SouthEast')

ht=text('Position',[5.42,-3.3,0],'String','N=32');
set(ht,'FontSize',12);
ht=text('Position',[5.17,-4.3,0],'String','N=64');
set(ht,'FontSize',12);
ht=text('Position',[4.85,-4.6,0],'String','N=128');
set(ht,'FontSize',12);
