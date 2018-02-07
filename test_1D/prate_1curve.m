clc; clear all; close all;

N=[10 100 500 1000];
dx=1./N;
e=[4.3740e-1 8.0641e-4 1.2998e-6 8.1303e-8]

ldx=log10(dx);
leh=log10(e);



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

% %% courbes...
hlf1=plot(tt,lsq1,'-r');
set(hlf1,'LineWidth',2.0);

%% legend
legend(hlf1,{['slope = ' num2str(a11(1))]},'Location','SouthEast')

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