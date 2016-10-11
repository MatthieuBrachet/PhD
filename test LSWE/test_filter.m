%% test new filter

clc; clear all; close all; format short;

global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global hp gp radius omega
global teta0 teta1

%% ************************************************************************
% Resolution de LSWEC sur la Cubed-Sphere.
%
% *** options :
% test : 0, ..., 5 choix du test à lancer.
%          test = 0 : solution stationnaire de type Galewski
%          test = 5 : solution en exp(-sigma*t) (forcage/ammortissement).
% video : 'yes' ou 'no', faire une video ou non.
% sauvegarde = 0 (ne rien sauvegarder), 1 (sauvegarder toutes les valeurs
%          finales).
% opt_ftr : filtre explicite de S. Redonnet (=2, 4, 6, 8, 10, ordre du
%          filtre) (=0, pas de filtrage)
% type_ftr : caracteristic or classic.
%
%% ************************************************************************

test=0;
opt_ftr=6;
n=40;
teta0=-3*pi/16;
teta1=3*pi/16;
mod72

cgrav=sqrt(gp*hp);
ccor=radius*omega;
c=max(cgrav,ccor);
cfl=0.5;
ddt=radius*dxi*cfl/c;
ndaymax=2;
Tmax=ndaymax*3600*24;
itermax=1;

tstart=cputime;
%% *** initialisation des données
t=0;
[ ht_fI,    vt_fI] = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII] = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV] = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

%% Filtrage
[htf1_fI, htf1_fII, htf1_fIII, htf1_fIV, htf1_fV, htf1_fVI, vtf_fI, vtf_fII, vtf_fIII, vtf_fIV, vtf_fV, vtf_fVI ]= ftrcar72(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);
[htf2_fI, htf2_fII, htf2_fIII, htf2_fIV, htf2_fV, htf2_fVI]=ftr72(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);


figure(1)
plot_cs15(n,nn,ht_fI-htf1_fI,ht_fII-htf1_fII,ht_fIII-htf1_fIII,ht_fIV-htf1_fIV,ht_fV-htf1_fV,ht_fVI-htf1_fVI);
title('new method')


figure(2)
plot_cs15(n,nn,ht_fI-htf2_fI,ht_fII-htf2_fII,ht_fIII-htf2_fIII,ht_fIV-htf2_fIV,ht_fV-htf2_fV,ht_fVI-htf2_fVI);
title('old method')


