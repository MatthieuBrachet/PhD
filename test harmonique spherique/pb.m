clc; clear all; close all

global n nn;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI;
global opt_ftr;

%% *** Benchmarks data ****************************************************
n=31;
nn=n+2;
l=12;
m=4;   % be carreful : [-l <= m <= l]
%% *** Opitons ************************************************************
opt_ftr=10;

%% ************************************************************************
mod_1b

 
[h_fI]=fun(x_fI,y_fI,z_fI,m,l);
[h_fII]=fun(x_fII,y_fII,z_fII,m,l);
[h_fIII]=fun(x_fIII,y_fIII,z_fIII,m,l);
[h_fIV]=fun(x_fIV,y_fIV,z_fIV,m,l);
[h_fV]=fun(x_fV,y_fV,z_fV,m,l);
[h_fVI]=fun(x_fVI,y_fVI,z_fVI,m,l);

figure(1);
plot_cs5(n,nn,h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI);
colorbar;
title(['harmonique l=' num2str(l) '; m=' num2str(m)])