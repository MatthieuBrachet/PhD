%% ************************************************************************
% test divergence on sphere (radius=1)
%% ************************************************************************
clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr scheme


opt_ftr='redonnet10';
scheme='compact4';

n=31;
mod101

[mfunfI,dfunI]=fun6(x_fI,y_fI,z_fI);
[mfunfII,dfunII]=fun6(x_fII,y_fII,z_fII);
[mfunfIII,dfunIII]=fun6(x_fIII,y_fIII,z_fIII);
[mfunfIV,dfunIV]=fun6(x_fIV,y_fIV,z_fIV);
[mfunfV,dfunV]=fun6(x_fV,y_fV,z_fV);
[mfunfVI,dfunVI]=fun6(x_fVI,y_fVI,z_fVI);

[div_fI,div_fII,div_fIII,div_fIV,div_fV,div_fVI]=...
    div101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);

% calcul erreur finale
errg_fI=max(max(abs(div_fI-dfunI)));
errg_fII=max(max(abs(div_fII-dfunII)));
errg_fIII=max(max(abs(div_fIII-dfunIII)));
errg_fIV=max(max(abs(div_fIV-dfunIV)));
errg_fV=max(max(abs(div_fV-dfunV)));
errg_fVI=max(max(abs(div_fVI-dfunVI)));

err_div=max([errg_fI errg_fII errg_fIII errg_fIV errg_fV errg_fVI])./max(max([dfunI dfunII dfunIII dfunIV dfunV dfunVI]))

figure(1)
plot_cs11(n,nn,dfunI,dfunII,dfunIII,dfunIV,dfunV,dfunVI)

figure(2)
plot_cs102(n,nn,dfunI,dfunII,dfunIII,dfunIV,dfunV,dfunVI)

fig_placier

