clc; clear all; close all;
format long

global n nn radius
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr scheme filtre

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';

n=1023; % for snapshot and better spherical integration (B. Portenelle works), n must be odd !
mod101
disp('mod74 : ok')

test=2;
[fun_I    ,int] = fun_quad(x_fI  ,y_fI  ,z_fI   ,test);
[fun_II   ,int] = fun_quad(x_fII ,y_fII ,z_fII  ,test);
[fun_III  ,int] = fun_quad(x_fIII,y_fIII,z_fIII ,test);
[fun_IV   ,int] = fun_quad(x_fIV ,y_fIV ,z_fIV  ,test);
[fun_V    ,int] = fun_quad(x_fV  ,y_fV  ,z_fV   ,test);
[fun_VI   ,int] = fun_quad(x_fVI ,y_fVI ,z_fVI  ,test);

figure(1)
plot_cs11(n,nn,(fun_I),(fun_II),(fun_III),(fun_IV),(fun_V),(fun_VI))
colorbar

figure(2)
plot_cs102(n,nn,(fun_I),(fun_II),(fun_III),(fun_IV),(fun_V),(fun_VI))
colorbar

figure(3)
plot_cs104(n,nn,(fun_I),(fun_II),(fun_III),(fun_IV),(fun_V),(fun_VI))
colorbar

figure(4)
surf(fun_I)

str='int'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg1]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
(nrmg1-int)/int

str='test'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg2]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
(nrmg2-int)/int

str='simpson'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg3]=nrm101(fun_I,fun_II,fun_III,fun_IV,fun_V,fun_VI,n,nn,str);
(nrmg3-int)/int

fig_placier




