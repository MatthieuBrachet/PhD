clc; clear all; close all;
radius=6371220;
N=log10(.5*pi*radius./[16 32 64 128]);
eh1=log10([.77612 .055566 .0030761 2.1494*10^-4]);
eh2=log10([.90855 .091356 .0065756 3.0062*10^-4]);
ehi=log10([2.6219 .38964 .039668 2.1993*10^-3]);


%% curve
hl1=plot(N,eh1,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerFaceColor','r');
set(hl1,'MarkerEdgeColor','k'); hold on

hl1=plot(N,eh2,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerFaceColor','b');
set(hl1,'MarkerEdgeColor','k');

hl1=plot(N,ehi,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerFaceColor','m');
set(hl1,'MarkerEdgeColor','k');

hold on;
[a11,b1]=polyfit(N,eh1,1);
tt=linspace(min(N),max(N),1000);
lsq1=a11(2)+a11(1)*tt;

hold on;
[a12,b1]=polyfit(N,eh2,1);
tt=linspace(min(N),max(N),1000);
lsq2=a12(2)+a12(1)*tt;


hold on;
[a13,b1]=polyfit(N,ehi,1);
tt=linspace(min(N),max(N),1000);
lsq3=a13(2)+a13(1)*tt;

legend(['norm 1   : slope = ' num2str(a11(1))],['norm 2   : slope = ' num2str(a12(1))],['norm inf : slope = ' num2str(a13(1))],'Location','SouthEast')


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

hlf1=plot(tt,lsq2,'-k');
set(hlf1,'LineWidth',2.0);

hlf1=plot(tt,lsq3,'-k');
set(hlf1,'LineWidth',2.0);

print('-dpng', ['rate_h_LSWE1.png'])