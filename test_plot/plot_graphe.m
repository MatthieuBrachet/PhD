clc; clear all; close all;timest=clock;
format long

global n
global opt_ftr scheme
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;

opt_ftr='redonnet10';
scheme='compact4';

n=31;
mod103

[ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI);
[ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII);
[ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV);
[ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV);
[ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI);
vmin=min(min([ht_fI ht_fII ht_fIII ht_fIV ht_fV ht_fVI]));
vmax=max(max([ht_fI ht_fII ht_fIII ht_fIV ht_fV ht_fVI]));
v=linspace(vmin,vmax,10);

figure(1)
plot_cs100(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI)
colorbar

figure(2)
plot_cs101(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI,v)

figure(3)
plot_cs102(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI)
colorbar

figure(4)
plot_cs103(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI,v)
colorbar

figure(5)
plot_cs11(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI)
colorbar

figure(6)
plot_cs104(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI)
colorbar

figure(7)
plot_cs105(n,nn,ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI)
colorbar


fig_placier