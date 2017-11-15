clc; clear all; close all;

radius=6371220;
N=[16 32 64 128 256 512];
dx=2*pi*radius./(4*N);
e=[1.6810*10^-6 1.1083*10^-7 7.0272*10^-9 4.4093*10^-10 2.7597*10^-11 1.7256*10^-12];
e=[7.2201*10^-6 4.4825*10^-7 2.7976*10^-8 1.7475*10^-9 1.0921*10^-10 6.8256*10^-12];
e=[4.5330*10^-5 2.9233*10^-6 1.8413*10^-7 1.1526*10^-8 7.2074*10^-10 4.5050*10^-11];
e=[9.1275*10^-7 5.4669*10^-8 3.3774*10^-9 2.1023*10^-10 1.3121*10^-11 8.1956*10^-13];
e=[9.2327*10^-7 5.4798*10^-8 3.3815*10^-9 2.1037*10^-10 1.3125*10^-11 8.1976*10^-13]


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
hlf1=plot(tt,lsq1,'-m');
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