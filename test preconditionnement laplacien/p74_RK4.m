%% ************************************************************************
% Resolution de SWEC sur la Cubed-Sphere.
% *** options :
% test = 0 : test 2 of Williamson & al.,
%        1 : test 5 of Williamson & al..
%        2 : test 5 of Williamson with smooth mountain,
%        -1 : test maison avec véritables reliefs,
% scheme : numerical spatial scheme used. 
% video : 'yes' ou 'no', do a video or not,
% nper  :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet) or implicit (visbal) filtering.
% delta_ftr : is between 0 and 1. If the filter of u is note Fu, then, the
%             filtering action is :
%                        (1-delta_ftr)*u + delta_ftr*Fu
%% ************************************************************************
clc; clear all; close all;
format long

global n nn
global scheme
% -----------------------------------
global G11_fI G12_fI G22_fI
global G11_fII G12_fII G22_fII
global G11_fIII G12_fIII G22_fIII
global G11_fIV G12_fIV G22_fIV
global G11_fV G12_fV G22_fV
global G11_fVI G12_fVI G22_fVI

scheme='compact4';
snapshot='yes';

n=51; % for snapshot, n must be in the form 2^m-1 !
ndaymax=3;
mod74

figure(1)
plot_cs11(n,nn,G11_fI,G11_fII,G11_fIII,G11_fIV,G11_fV,G11_fVI)
title('g^{\xi, \xi}')

figure(2)
plot_cs11(n,nn,G12_fI,G12_fII,G12_fIII,G12_fIV,G12_fV,G12_fVI)
title('g^{\xi, \eta}')

figure(3)
plot_cs11(n,nn,G22_fI,G22_fII,G22_fIII,G22_fIV,G22_fV,G22_fVI)
title('g^{\eta, \eta}')

figure(4)
plot_cs11(n,nn,G12_fI./G11_fI,G12_fII./G11_fII,G12_fIII./G11_fIII,G12_fIV./G11_fIV,G12_fV./G11_fV,G12_fVI./G11_fVI)
title('ordre de grandeur g^{\xi, \eta} / g^{\xi, \xi}')

str='infty';
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,G12G11]=nrm74(G12_fI./G11_fI,G12_fII./G11_fII,G12_fIII./G11_fIII,G12_fIV./G11_fIV,G12_fV./G11_fV,G12_fVI./G11_fVI,n,nn,str);
G12G11

figure(5)
plot_cs11(n,nn,G12_fI./G22_fI,G12_fII./G22_fII,G12_fIII./G22_fIII,G12_fIV./G22_fIV,G12_fV./G22_fV,G12_fVI./G22_fVI)
title('ordre de grandeur g^{\xi, \eta} / g^{\eta, \eta}')

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,G12G22]=nrm74(G12_fI./G22_fI,G12_fII./G22_fII,G12_fIII./G22_fIII,G12_fIV./G22_fIV,G12_fV./G22_fV,G12_fVI./G22_fVI,n,nn,str);
G12G22

fig_placier

