clc; clear all; close all; format shorte
global n nn radius
global scheme
global opt_ftr

scheme='compact4';
opt_ftr='redonnet2';

NN=[31 63 127 255 511];
for l=1:length(NN)
    n=NN(l);
    mod101
    disp('mod101 : ok')

    [mfunfI]   = fun(x_fI,y_fI,z_fI);
    [mfunfII]  = fun(x_fII,y_fII,z_fII);
    [mfunfIII] = fun(x_fIII,y_fIII,z_fIII);
    [mfunfIV]  = fun(x_fIV,y_fIV,z_fIV);
    [mfunfV]   = fun(x_fV,y_fV,z_fV);
    [mfunfVI]  = fun(x_fVI,y_fVI,z_fVI);

    [funfI,funfII,funfIII,funfIV,funfV,funfVI]=ftr_mixte101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI, n, nn);

    errI   = abs(funfI-mfunfI);
    errII  = abs(funfII-mfunfII);
    errIII = abs(funfIII-mfunfIII);
    errIV  = abs(funfIV-mfunfIV);
    errV   = abs(funfV-mfunfV);
    errVI  = abs(funfVI-mfunfVI);


    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'1');
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'1');
    e1(l)=nrmg./MMM;

    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'2');
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'2');
    e2(l)=nrmg./MMM;

    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'infty');
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'infty');
    e3(l)=nrmg./MMM;
end






N=NN+1;
dx=2*pi*radius./(4*N);

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
set(hl1,'MarkerEdgeColor','k');

hold on;
[a1,b1]=polyfit(ldx,le1,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a1(2)+a1(1)*tt;
% %% courbes...
%hlf1=plot(tt,lsq1,'-r');
hlf1=plot(tt,lsq1,':k');
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
%hlf2=plot(tt,lsq1,'-b');
hlf2=plot(tt,lsq1,'--k');
set(hlf2,'LineWidth',2.0);

%% curve 3
hl1=plot(ldx,le3,'ok');
set(hl1,'LineWidth',2.0);
set(hl1,'MarkerSize',10);
set(hl1,'MarkerEdgeColor','k');

hold on;
[a3,b1]=polyfit(ldx,le3,1);
tt=linspace(min(ldx),max(ldx),10);
lsq1=a3(2)+a3(1)*tt;
% %% courbes...
%hlf3=plot(tt,lsq1,'-g');
hlf3=plot(tt,lsq1,'-k');
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
legend([hlf1,hlf2,hlf3],{['norm 1 - slope = ' num2str(a1(1))],['norm 2 - slope = ' num2str(a2(1))],['norm \infty - slope = ' num2str(a3(1))]},'Location','SouthEast')

