% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
clc; close all; clear all;

global n nn;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global scheme opt_ftr opt_detec opt_ftr1


scheme='compact4';
opt_ftr='redonnet10';
opt_detec='redonnet10';
opt_ftr1='redonnet6';

n=63;
mod101; 

[funfI]=fun4(x_fI,y_fI,z_fI);
[funfII]=fun4(x_fII,y_fII,z_fII);
[funfIII]=fun4(x_fIII,y_fIII,z_fIII);
[funfIV]=fun4(x_fIV,y_fIV,z_fIV);
[funfV]=fun4(x_fV,y_fV,z_fV);
[funfVI]=fun4(x_fVI,y_fVI,z_fVI);

[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_mixte101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

errg_fI=funftI-funfI;
errg_fII=funftII-funfII;
errg_fIII=funftIII-funfIII;
errg_fIV=funftIV-funfIV;
errg_fV=funftV-funfV;
errg_fVI=funftVI-funfVI;


figure(1)
surf(errg_fI)
title('panel I')

figure(2)
surf(errg_fII)
title('panel II')

figure(3)
surf(errg_fIII)
title('panel III')

figure(4)
surf(errg_fIV)
title('panel IV')

figure(5)
surf(errg_fV)
title('panel V')

figure(6)
surf(errg_fVI)
title('panel VI')

figure(7)
plot_cs11(n,nn,errg_fI,errg_fII,errg_fIII,errg_fIV,errg_fV,errg_fVI)
colorbar

figure(8)
plot_cs11(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI)
colorbar

figure(9)
plot_cs11(n,nn,funftI,funftII,funftIII,funftIV,funftV,funftVI)
colorbar
fig_placier