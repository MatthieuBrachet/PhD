clc; clear all; close all;

N=[50 100 500 1000];
dx=1./N;
e1=[1.0222e-3 6.5636e-5 1.0551e-7 6.5999e-9];
e2=[1.3268e-3 8.4902e-5 1.3640e-7 8.5324e-9];

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
hlf1=plot(tt,lsq1,'-k');
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


% ylabel
ha=gca;
ya=get(gca,'Ylabel');
set(ha,'Ygrid','on');
set(ya,'String','Log_{10}(error)');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

%% legend
legend([hlf1,hlf2],{['norm 1 - slope = ' num2str(a1(1))],['norm 2 - slope = ' num2str(a2(1))]},'Location','SouthEast')
