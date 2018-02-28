clc; clear all; close all;

radius=6371220;
N=[16 32 64 128 256 512];
dx=2*pi*radius./(4*N);
e=[1.8636*10^-4 1.3545*10^-5 9.7564*10^-7 6.5592*10^-8 4.2563*10^-9 2.7115*10^-10];
%e=[3.1368*10^-4 2.3910*10^-5 1.6716*10^-6 1.1088*10^-7 7.1437*10^-9 4.5361*10^-10];


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


