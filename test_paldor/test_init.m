%% ************************************************************************
% Solve LSWEC on the Cubed-Sphere.
% *** options :
% scheme : numerical spatial scheme used. 
% video  : 'yes' ou 'no', do a video or not,
% nper   :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet10, redonnet8, ...).
%% ************************************************************************
clc; clear all; close all;timest=clock;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr scheme nrm
global waveFlag

comment='.';
video = 'no';
sauvegarde = 1;
opt_ftr='redonnet10';
scheme='compact4';
snapshot='yes';
nrm='int';
waveFlag='EIG';

n=63;
ndaymax=.1;
cfl=0.9;
mod101
disp('mod101 : ok')

%% ************************************************************************
time=0;
[ht_fI,vt_fI] = sol_exacte(x_fI,y_fI,z_fI,time);
[ht_fII,vt_fII] = sol_exacte(x_fII,y_fII,z_fII,time);
[ht_fIII,vt_fIII] = sol_exacte(x_fIII,y_fIII,z_fIII,time);
[ht_fIV,vt_fIV] = sol_exacte(x_fIV,y_fIV,z_fIV,time);
[ht_fV,vt_fV] = sol_exacte(x_fV,y_fV,z_fV,time);
[ht_fVI,vt_fVI] = sol_exacte(x_fVI,y_fVI,z_fVI,time);


figure(1)
plot_cs100(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
colorbar



[up_fI,up_fII,up_fIII,up_fIV,up_fV,up_fVI] = plot_u(vt_fI,vt_fII,vt_fIII,...
    vt_fIV,vt_fV,vt_fVI);
figure(2)
plot_cs100(n,nn,up_fI,up_fII,up_fIII,up_fIV,up_fV,up_fVI);
colorbar