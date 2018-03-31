clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global dxi corner
global opt_ftr scheme filtre
global ang1 ang2

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';
corner = 1;

test=1;
n=31;
ang1=2*pi*rand(1);
ang2=2*pi*rand(1);

mod101;

[fun_I    ,int] = fun_quad2(x_fI  ,y_fI  ,z_fI   ,test);
[fun_II   ,int] = fun_quad2(x_fII ,y_fII ,z_fII  ,test);
[fun_III  ,int] = fun_quad2(x_fIII,y_fIII,z_fIII ,test);
[fun_IV   ,int] = fun_quad2(x_fIV ,y_fIV ,z_fIV  ,test);
[fun_V    ,int] = fun_quad2(x_fV  ,y_fV  ,z_fV   ,test);
[fun_VI   ,int] = fun_quad2(x_fVI ,y_fVI ,z_fVI  ,test);
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,'infty');
nrmg

str='int';
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg1]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
abs((nrmg1-int))/int;

figure(2)
plot_cs11(n,nn,fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI)
colorbar

ang1=0;
ang2=0;
[fun_I    ,int] = fun_quad2(x_fI  ,y_fI  ,z_fI   ,test);
[fun_II   ,int] = fun_quad2(x_fII ,y_fII ,z_fII  ,test);
[fun_III  ,int] = fun_quad2(x_fIII,y_fIII,z_fIII ,test);
[fun_IV   ,int] = fun_quad2(x_fIV ,y_fIV ,z_fIV  ,test);
[fun_V    ,int] = fun_quad2(x_fV  ,y_fV  ,z_fV   ,test);
[fun_VI   ,int] = fun_quad2(x_fVI ,y_fVI ,z_fVI  ,test);

figure(3)
plot_cs11(n,nn,fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI)
colorbar
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,'infty');
nrmg