%% ************************************************************************
% Solve SWEC on the Cubed-Sphere.
% *** options :
% test = 0 : test 2 of Williamson & al.,
%        1 : test 5 of Williamson & al..
%        2 : test 5 of Williamson with smooth mountain,
%        3 : stationnary Galewsky (exp),
%        4 : Galewsky with perturbation (exp),
%        -1 : test with Earth topography
% scheme : numerical spatial scheme used. 
% video : 'yes' ou 'no', do a video or not,
% nper  :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet) or adaptative filtering.

%% ************************************************************************
clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test scheme


comment='.';
test=1;
video = 'no';
nper=1;
sauvegarde = 0;
filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';
snapshot='yes';

n=63; % for snapshot, n must be odd !
mod101

[vt_I,vortvt_I] = sol_exacte_vort(x_fI,y_fI,z_fI);
[vt_II,vortvt_II] = sol_exacte_vort(x_fII,y_fII,z_fII);
[vt_III,vortvt_III] = sol_exacte_vort(x_fIII,y_fIII,z_fIII);
[vt_IV,vortvt_IV] = sol_exacte_vort(x_fIV,y_fIV,z_fIV);
[vt_V,vortvt_V] = sol_exacte_vort(x_fV,y_fV,z_fV);
[vt_VI,vortvt_VI] = sol_exacte_vort(x_fVI,y_fVI,z_fVI);

[vort_I,vort_II,vort_III,vort_IV,vort_V,vort_VI]=...
    vort101(vt_I,vt_II,vt_III,vt_IV,vt_V,vt_VI,n,nn);

err_I=vort_I-vortvt_I;
err_II=vort_II-vortvt_II;
err_III=vort_III-vortvt_III;
err_IV=vort_IV-vortvt_IV;
err_V=vort_V-vortvt_V;
err_VI=vort_VI-vortvt_VI;

figure(1)
plot_cs100(n,nn,err_I,err_II,err_III,err_IV,err_V,err_VI);
colorbar
title('err.')

figure(2)
plot_cs100(n,nn,vort_I,vort_II,vort_III,vort_IV,vort_V,vort_VI);
colorbar
title('calc.')

figure(3)
plot_cs100(n,nn,vortvt_I,vortvt_II,vortvt_III,vortvt_IV,vortvt_V,vortvt_VI);
colorbar
title('exact')

fig_placier