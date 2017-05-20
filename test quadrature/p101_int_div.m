clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr scheme filtre test

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';
test=1;

n=31; % for snapshot and better spherical integration (B. Portenelle works), n must be odd !
mod101
disp('mod74 : ok')

t=0;
[vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI);
[vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII);
[vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII);
[vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV);
[vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV);
[vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI);

[div_I,div_II,div_III,div_IV,div_V,div_VI]=div102(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,n,nn);

%[div_I,div_II,div_III,div_IV,div_V,div_VI]=ftr_mixte101(div_I,div_II,div_III,div_IV,div_V,div_VI,n,nn);

figure(1)
plot_cs11(n,nn,div_I,div_II,div_III,div_IV,div_V,div_VI)

str='int'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg1]=nrm101(div_I,div_II,div_III,div_IV,div_V,div_VI,n,nn,str);
nrmg1

str='cons_int'
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg2]=nrm101(div_I,div_II,div_III,div_IV,div_V,div_VI,n,nn,str);
nrmg2

