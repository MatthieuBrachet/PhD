clc; clear all; close all; format short;
format long;

global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test
global hp gp u0 radius omega
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
%
%% ************************************************************************


test=0;
opt_ftr=10;
n=30;
teta0=-3*pi/16;
teta1=3*pi/16;
mod72

cgrav=sqrt(gp*hp);
ccor=radius*omega;
c=max(cgrav,ccor);

cfl=0.5;
ddt=radius*dxi*cfl/c;
%ndaymax=2;
% JPC
ndaymax=3;
Tmax=ndaymax*3600*24;
itermax=5000;

comment='no filtering, (hnew-h) au lieu de (hnew-h)/h';

tstart=cputime;
mm=-5; MM=0;
%% *** initialisation des données
t=0;
[ ht_fI,    vt_fI] = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII] = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV] = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV] = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI] = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

[adve_fI]=adv_exacte(x_fI,y_fI,z_fI);
[adve_fII]=adv_exacte(x_fII,y_fII,z_fII);
[adve_fIII]=adv_exacte(x_fIII,y_fIII,z_fIII);
[adve_fIV]=adv_exacte(x_fIV,y_fIV,z_fIV);
[adve_fV]=adv_exacte(x_fV,y_fV,z_fV);
[adve_fVI]=adv_exacte(x_fVI,y_fVI,z_fVI);

[adv_fI,adv_fII,adv_fIII, adv_fIV, adv_fV, adv_fVI]=adv72(vt_fI,vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

e_fI=max(max(max(abs(adv_fI-adve_fI))));
e_fII=max(max(max(abs(adv_fII-adve_fII))));
e_fIII=max(max(max(abs(adv_fIII-adve_fIII))));
e_fIV=max(max(max(abs(adv_fIV-adve_fIV))));
e_fV=max(max(max(abs(adv_fV-adve_fV))));
e_fVI=max(max(max(abs(adv_fVI-adve_fVI))));

E=max([e_fI, e_fII, e_fIII, e_fIV, e_fV, e_fVI])



