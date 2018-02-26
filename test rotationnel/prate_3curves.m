clc; clear all; close all;

radius=6371220;
N=[8 16 32 64 128 256];
% N=[40 50 60 80 100 150];
dx=2*pi*radius./(4*N);
% e1=[1.9936e-2 1.1325e-3 6.8262e-5 4.2292e-6 2.6413e-7 1.6559e-8];
% e2=[1.8717e-2 1.0216e-3 6.1471e-5 3.8166e-6 2.3881e-7 1.4975e-8];
% e3=[2.5623e-2 1.4497e-3 8.9084e-5 5.5817e-6 3.5265e-7 2.2194e-8];
% 
% e1=[1.6305e-2 9.0026e-4 5.4527e-5 3.3854e-6 2.1153e-7 1.3242e-8];
% e2=[1.5494e-2 8.1833e-4 4.9215e-5 3.0551e-6 1.9108e-7 1.1976e-8];
% e3=[1.7473e-2 9.5002e-4 7.3232e-5 5.1367e-6 3.3731e-7 2.1624e-8];
% 
% e1=[3.4770e-4 1.6686e-5 1.0089e-6 6.3688e-8 4.0283e-9 2.5396e-10];
% e2=[1.4286e-4 7.0948e-6 4.1258e-7 2.5581e-8 1.6143e-9 1.0283e-10];
% e3=[1.1470e-4 1.0256e-5 6.2778e-7 3.7403e-8 2.2653e-9 1.9459e-10];
% 
% 
% %% projetée
% e1=[1.0377e-4 6.3236e-6 3.9444e-7 2.4726e-8 1.5500e-9 9.7139e-11];
% e2=[1.2588e-4 7.4682e-6 4.6278e-7 2.8931e-8 1.8111e-9 1.1342e-10];
% e3=[3.6670e-4 2.2222e-5 1.3713e-6 8.5415e-8 5.3339e-9 3.3330e-10];
% 
% %% zonale
% e1=[2.9158e-4 1.7719e-5 1.1025e-6 6.9056e-8 4.3244e-9 2.7061e-10];
% e2=[3.3039e-4 1.9906e-5 1.2416e-6 7.7821e-8 4.8755e-9 3.0522e-10];
% e3=[6.7103e-4 4.0648e-5 2.5207e-6 1.6433e-7 1.0822e-8 6.9474e-10];
% 
% e1=[4.3694e-3 2.0968e-4 .2678e-5 8.0033e-7 5.0621e-8 3.1914e-9];
% e2=[2.8177e-10 1.3994e-11 8.1375e-13 5.0456e-14 3.1840e-15 2.0282e-16];
% e3=[3.5509e-17 3.1748e-18 1.9435e-19 1.1579e-20 7.0129e-22 6.0241e-23];
% 
% e1=[4.4832e-3 2.3244e-4 1.4639e-5 9.3537e-7 5.9638e-8 3.7867e-9];
% e2=[2.9091e-10 1.5727e-11 9.3911e-13 5.8935e-14 3.7494e-15 2.4082e-16];
% e3=[3.3243e-17 2.4699e-18 1.8005e-19 1.1002e-20 6.8031e-22 7.5741e-23];
N=[8 16 32 64 128 256];
dx=2*pi*radius./(4*N);
e1=[6.7655e-1 5.5017e+1 8.0953e+0 4.3791e-1 2.5533e-2 1.5590e-3]./(4*pi*radius*radius);
e2=[5.1432e-8 3.3097e-6 4.7068e-7 2.4642e-8 1.4237e-9 8.6829e-11]./sqrt(4*pi*radius*radius);
e3=[6.6711e-15 4.9935e-13 4.5428e-14 2.6199e-15 1.3325e-16 9.4033e-18];

% N=[16 32 64 128 256];
% dx=2*pi*radius./(4*N);
% e1=[5.5017e+1 8.0953e+0 4.3791e-1 2.5533e-2 1.5590e-3]./(4*pi*radius*radius);
% e2=[3.3097e-6 4.7068e-7 2.4642e-8 1.4237e-9 8.6829e-11]./sqrt(4*pi*radius*radius);
% e3=[4.9935e-13 4.5428e-14 2.6199e-15 1.3325e-16 9.4033e-18];

ldx=log10(dx);
le1=log10(e1);
le2=log10(e2);
le3=log10(e3);

%% plot 1
figure(11); 
%% curve 1
hl1=plot(ldx,le1,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerFaceColor','m');
%set(hl1,'MarkerEdgeColor','k');

hold on;
[a1,b1]=polyfit(ldx,le1,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a1(2)+a1(1)*tt;
a1=[3.8523 -35.1389]; lsq1=a1(2)+a1(1)*tt;
% %% courbes...
%hlf1=plot(tt,lsq1,'-r');
hlf1=plot(tt,lsq1,'k');
set(hlf1,'LineWidth',2.0);


%% curve 2
hl2=plot(ldx,le2,'ok');
set(hl2,'LineWidth',2.0);
set(hl2,'MarkerSize',10);
set(hl2,'MarkerFaceColor','b');
%set(hl2,'MarkerEdgeColor','k');

hold on;
[a2,b1]=polyfit(ldx,le2,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a2(2)+a2(1)*tt;
a2=[3.8805 -35.1747]; lsq1=a2(2)+a2(1)*tt;
% %% courbes...
%hlf2=plot(tt,lsq1,'-b');
hlf2=plot(tt,lsq1,'k');
set(hlf2,'LineWidth',2.0);

%% curve 3
hl3=plot(ldx,le3,'ok');
set(hl3,'LineWidth',2.0);
set(hl3,'MarkerSize',10);
set(hl3,'MarkerFaceColor','r');
%set(hl3,'MarkerEdgeColor','k');

hold on;
[a3,b1]=polyfit(ldx,le3,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a3(2)+a3(1)*tt;
a3=[3.9806 -35.3017]; lsq1=a3(2)+a3(1)*tt;
% %% courbes...
%hlf3=plot(tt,lsq1,'-g');
hlf3=plot(tt,lsq1,'k');
set(hlf3,'LineWidth',2.0);

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

%% legend
legend([hl1,hl2,hl3],{['norm 1 - slope = ' num2str(a1(1))],['norm 2 - slope = ' num2str(a2(1))],['norm \infty - slope = ' num2str(a3(1))]},'Location','NorthWest')

%% texte
ht=text('Position',[4.5,-16,0],'String','N=256');
set(ht,'FontSize',12);
ht=text('Position',[4.8,-15,0],'String','N=128');
set(ht,'FontSize',12);
ht=text('Position',[5.1,-14,0],'String','N=64');
set(ht,'FontSize',12);
ht=text('Position',[5.4,-12.5,0],'String','N=32');
set(ht,'FontSize',12);
ht=text('Position',[5.7,-11.6,0],'String','N=16');
set(ht,'FontSize',12);
ht=text('Position',[6,-13.5,0],'String','N=8');
set(ht,'FontSize',12);