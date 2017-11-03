clc; clear all; close all;

radius=6.37122d+06;
N=log10(.5*pi*radius./[16 32 64 128]);
ev=log10([3.5840 1.1995*10^-1 2.07*10^-2 2.1465*10^-3]);



%% curve
hl1=plot(N,ev,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerFaceColor','r');
set(hl1,'MarkerEdgeColor','k'); hold on


hold on;
[a11,b1]=polyfit(N,ev,1);
tt=linspace(min(N),max(N),1000);
lsq1=a11(2)+a11(1)*tt;

legend(['norm inf : slope = ' num2str(a11(1))],'Location','SouthEast')


%% legend and axis

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
set(ya,'String','Log_{10}(error_{max})');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

% title
ha=gca;
xt=get(gca,'Title');
title('Convergence rate');
set (xt,'FontSize',12);


%% courbes...
hlf1=plot(tt,lsq1,'-k');
set(hlf1,'LineWidth',2.0);

print('-dpng', ['rate_v_LSWE1.png'])
