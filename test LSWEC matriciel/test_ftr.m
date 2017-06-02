%% ************************************************************************
% Resolution de LSWEC sur la Cubed-Sphere.
%
% *** options :
% test :   test = 0 : solution stationnaire de type Galewsky
%          test = 1 : solution en exp(-sigma*t) (forcage/ammortissement).
%          test = 2 : N. Paldor test case.
%          test = 3 : personnal test case.
% video : 'yes' ou 'no', faire une video ou non.
% sauvegarde = 0 (ne rien sauvegarder), 1 (sauvegarder toutes les valeurs
%          finales).
% opt_ftr : filtre explicite de S. Redonnet (=2, 4, 6, 8, 10, ordre du
%          filtre) (=0, pas de filtrage)
% type_ftr : classic (on the xi and eta variables) or caracteristic 
%            variables.
%
%% ************************************************************************
clc; clear all; close all; format long;
global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test scheme
global hp gp radius omega
global teta0 teta1

test=0;
video = 'no';
sauvegarde = 0;
opt_ftr='redonnet10';
type_ftr='nonsymetric';
scheme='compact4';
n=63;
mod101

teta0=-pi/3;
teta1=pi/3;

cgrav=sqrt(gp*hp);
ccor=radius*omega;
c=max(cgrav,ccor);

cfl=0.9;
ddt=radius*dxi*cfl/c;
ndaymax=5;
Tmax=ndaymax*3600*24;
itermax=10000;

%% *** initialisation des données
t=0;
[ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

[ht_fI1, ht_fII1, ht_fIII1, ht_fIV1, ht_fV1, ht_fVI1]=ftr72(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);

[ht_fI2, ht_fII2, ht_fIII2, ht_fIV2, ht_fV2, ht_fVI2]=ftr101(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, n, nn);


funfI=ht_fI1-ht_fI;
funfII=ht_fII1-ht_fII;
funfIII=ht_fIII1-ht_fIII;
funfIV=ht_fIV1-ht_fIV;
funfV=ht_fV1-ht_fV;
funfVI=ht_fVI1-ht_fVI;

figure(1)
plot_cs11(n,nn,ht_fI1-ht_fI,ht_fII1-ht_fII,ht_fIII1-ht_fIII,ht_fIV1-ht_fIV,ht_fV1-ht_fV,ht_fVI1-ht_fVI)

figure(2)
plot_cs11(n,nn,ht_fI2-ht_fI,ht_fII2-ht_fII,ht_fIII2-ht_fIII,ht_fIV2-ht_fIV,ht_fV2-ht_fV,ht_fVI2-ht_fVI)

figure(3)
plot_cs11(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI)

figure(4)
surf(funfI)

fig_placier