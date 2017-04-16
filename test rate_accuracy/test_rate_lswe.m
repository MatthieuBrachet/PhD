clc; clear all; close all;

radius=6371220;
N=[32 64 128];
dx=radius./(4*N);

%% LSWE test stationnaire
% eh=[3.8953*10^-4 2.3646*10^-5 6.7734*10^-7];
% ev=[2.0264*10^-3 2.0962*10^-4 1.0382*10^-5];
%% LSWE test exp(-sigma*t)
eh=[2.4655*10^-3 1.2002*10^-3 5.9203*10^-4];
ev=[1.3575*10^-3 1.4654*10^-4 7.2371*10^-6];

ldx=log10(dx);
leh=log10(eh);
lev=log10(ev);



%% plot 1
figure(11); 
%% curve 1
hl1=plot(ldx,leh,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');

hold on;
[a11,b1]=polyfit(ldx,leh,1);
tt=linspace(min(ldx),max(ldx),1000);
lsq1=a11(2)+a11(1)*tt;


%% curve 2
hl1=plot(ldx,lev,'xk');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');

hold on;
[a12,b1]=polyfit(ldx,lev,1);
tt=linspace(min(ldx),max(ldx),1000);
lsq2=a12(2)+a12(1)*tt;



%% legend
legend(['norm on height   : slope = ' num2str(a11(1))],['norm on velocity : slope = ' num2str(a12(1))],'Location','SouthEast')

ha=gca;
xa=get(gca,'Xlabel');
set(ha,'Xgrid','on');
set(xa,'String','Log10(\Delta)');
set (xa,'FontName','Calibri');
set(xa,'FontSize',12);

% ylabel
ha=gca;
ya=get(gca,'Ylabel');
set(ha,'Ygrid','on');
set(ya,'String','Log10(error_{max})');
set (ya,'FontName','Calibri');
set(ya,'FontSize',12);

% title
ha=gca;
xt=get(gca,'Title');
title('Convergence rate');
set (xt,'FontSize',12);


% %% courbes...
hlf1=plot(tt,lsq1,'--k');
set(hlf1,'LineWidth',2.0);

hlf1=plot(tt,lsq2,'-k');
set(hlf1,'LineWidth',2.0);
