clc; clear all; close all;

radius=1;
N=[16 32 64 128];
dx=2*pi*radius./(4*N);
%e=[3.038e-7 8.3644e-9 9.6391e-11 8.4850e-12];
e=[5.5343e-8 5.8053e-10 6.3007e-12 7.0442e-14];


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
hlf1=plot(tt,lsq1,'-k');
set(hlf1,'LineWidth',2.0);

%% legend
legend(hlf1,{['pente = ' num2str(a11(1))]},'Location','SouthEast')

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


