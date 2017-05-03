clc; clear all; close all;timest=clock;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr scheme filtre test

filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';
test=1;

n=15; % for snapshot and better spherical integration (B. Portenelle works), n must be odd !
mod101
disp('mod74 : ok')

t=0;
[vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI);
[vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII);
[vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII);
[vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV);
[vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV);
[vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI);

tic
[div1_fI,div1_fII,div1_fIII,div1_fIV,div1_fV,div1_fVI]=div101(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,n,nn);
toc

tic
[div1_fI,div1_fII,div1_fIII,div1_fIV,div1_fV,div1_fVI]=div102(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,n,nn);
toc

tic
[div2_fI,div2_fII,div2_fIII,div2_fIV,div2_fV,div2_fVI]=div103(vt_fI,vt_fII,vt_fIII,vt_fIV,vt_fV,vt_fVI,n,nn);
toc

err_fI=abs(div1_fI-div2_fI);
err_fII=abs(div1_fII-div2_fII);
err_fIII=abs(div1_fIII-div2_fIII);
err_fIV=abs(div1_fIV-div2_fIV);
err_fV=abs(div1_fV-div2_fV);
err_fVI=abs(div1_fVI-div2_fVI);

figure(1)
plot_cs11(n,nn,err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI)

