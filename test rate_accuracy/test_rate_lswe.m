clc; clear all; close all;

radius=6371220;
N=[32 64 128];
dx=2*pi./(4*N);
test = 1
if test == 1
    %% LSWE test stationnaire h=10000
    eh=[5.1609*10^-4 2.2834*10^-5 9.9035*10^-7];
    ev=[1.5944*10^-3 1.5444*10^-4 6.8467*10^-6];
elseif test == 2
    %% LSWE test exp(-sigma*t) h=10
    eh=[1.5202*10^-3 1.4494*10^-4 1.3400*10^-5];
    ev=[1.4808*10^-3 1.4663*10^-4 6.4788*10^-6];
elseif test == 3
    %% LSWE test exp(-sigma*t) h=10000
    eh=[1.2725*10^-1 7.0431*10^-3 2.7248*10^-4];
    ev=[1.5352*10^-3 1.5060*10^-4 6.6822*10^-6];
end

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
