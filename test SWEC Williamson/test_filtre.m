%% ************************************************************************
% Resolution de SWEC sur la Cubed-Sphere.
% *** options :
% test = 0 : test 2 of Williamson & al.,
%        1 : test 5 of Williamson & al..
% scheme : numerical spatial scheme used. 
% video : 'yes' ou 'no', do a video or not,
% nper  :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet) or implicit (visbal) filtering.
% alfa_ftr : parameter for implicit fliter (type visbal only, 
%            if alpha_ftr=0, the filter is equivalent to Redonnet filter,
%            if alpha_ftr=0.5, the filter is inexistant).
% delta_ftr : is between 0 and 1. If the filter of u is note Fu, then, the
%             filtering action is :
%                        (1-delta_ftr)*u + delta_ftr*Fu
%
%% ************************************************************************
clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr alfa_ftr test scheme
global gp h0 u0 radius
global alpha

test=1;
video = 'no';
nper=1;
sauvegarde = 1;
opt_ftr='redonnet10';
alfa_ftr=0;
delta_ftr=1;
scheme='compact4';
snapshot='yes';

n=31; % for snapshot, n must be in the form 2^m-1 !
mod74

%% *** test data **********************************************************

if test == 0
    alpha=pi/7;
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

%% Filtrage
[htf2_fI, htf2_fII, htf2_fIII, htf2_fIV, htf2_fV, htf2_fVI]=ftr74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);

figure(2)
plot_cs15(n,nn,ht_fI-htf2_fI,ht_fII-htf2_fII,ht_fIII-htf2_fIII,ht_fIV-htf2_fIV,ht_fV-htf2_fV,ht_fVI-htf2_fVI);
title('old method')

