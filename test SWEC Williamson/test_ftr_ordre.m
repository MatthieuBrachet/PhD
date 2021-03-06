%% ************************************************************************
% verification de l'ordre des filtreages
% *************************************************************************
clc; clear all; close all;
format long

global n dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr opt_ftr1 test scheme
global gp h0 u0 radius
global alpha


test=1;
opt_ftr='redonnet10';
opt_ftr1='redonnet10';
scheme='compact4';

NN=[10;20;40;80;160];
E=[];
H=[];
for i=1:length(NN)
    clc; n=NN(i)
    mod74
    %% *** test data **********************************************************

    if test == 0
        alpha=0;
        u0=2*pi*radius/(12*24*3600);
        h0=2.94*10^4/gp;
    elseif test == 1
        alpha=0;
        u0=20;
        h0=5960;
    end

    %% *** initial data *******************************************************
    t=0;
    [ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
    [ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
    [ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
    [ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
    [ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
    [ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);
    
    
    [funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
            ftr_mixte74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI,n,nn);
        
    ee_fI=abs(funftI-ht_fI)./abs(ht_fI);
    ee_fII=abs(funftII-ht_fII)./abs(ht_fII);
    ee_fIII=abs(funftIII-ht_fIII)./abs(ht_fIII);
    ee_fIV=abs(funftIV-ht_fIV)./abs(ht_fIV);
    ee_fV=abs(funftV-ht_fV)./abs(ht_fV);
    ee_fVI=abs(funftVI-ht_fVI)./abs(ht_fVI);
    
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
            nrm74(ee_fI, ee_fII, ee_fIII, ee_fIV, ee_fV, ee_fVI ,n,nn,'infty');
        
    E=[E nrmg];
    H=[H dxi];
        
end

figure(1)
loglog(H,E,'o-',H,H.^10,'x-')
legend('numerical','theorical')
grid on

P=polyfit(log(H),log(E),1);
disp(P(1))

figure(2)
plot_cs11(n,nn,ee_fI,ee_fII,ee_fIII,ee_fIV,ee_fV,ee_fVI);
title('error')

figure(3)
plot_cs11(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI);
title('initial map')

figure(4)
plot_cs11(n,nn,funftI,funftII,funftIII,funftIV,funftV,funftVI);
title('filtered map')

fig_placier;
