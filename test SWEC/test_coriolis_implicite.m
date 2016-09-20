clc; clear all; close all

global n dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test scheme
global gp h0 u0 radius omega
global teta0 teta1


test=1;
video = 'no';
sauvegarde = 0;
opt_ftr=10;
scheme='compact4';
snapshot='no';

n=40;
teta0=pi/7;
teta1=pi/2-teta0;
mod74

ccor=radius*omega;
cgrav=sqrt(h0*gp);
cvit=u0;
c=max([cgrav,ccor,cvit]);

cfl=0.1;
ddt=radius*dxi*cfl/c;
ndaymax=2;
Tmax=ndaymax*3600*24;
itermax=10000;

comment='Start Galewsky benchmark.';

tstart=cputime;
%% *** initialisation des données *****************************************
t=0;
[ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);

alpha=200;
[f_fI] = thcoriolis_mat75(vt_fI,alpha,x_fI,y_fI,z_fI);
[u_fI] = coriolis_impli75(f_fI,alpha,x_fI,y_fI,z_fI);

[f_fII] = thcoriolis_mat75(vt_fII,alpha,x_fII,y_fII,z_fII);
[u_fII] = coriolis_impli75(f_fII,alpha,x_fII,y_fII,z_fII);

[f_fIII] = thcoriolis_mat75(vt_fIII,alpha,x_fIII,y_fIII,z_fIII);
[u_fIII] = coriolis_impli75(f_fIII,alpha,x_fIII,y_fIII,z_fIII);

[f_fIV] = thcoriolis_mat75(vt_fIV,alpha,x_fIV,y_fIV,z_fIV);
[u_fIV] = coriolis_impli75(f_fIV,alpha,x_fIV,y_fIV,z_fIV);

[f_fV] = thcoriolis_mat75(vt_fV,alpha,x_fV,y_fV,z_fV);
[u_fV] = coriolis_impli75(f_fV,alpha,x_fV,y_fV,z_fV);

[f_fVI] = thcoriolis_mat75(vt_fVI,alpha,x_fVI,y_fVI,z_fVI);
[u_fVI] = coriolis_impli75(f_fVI,alpha,x_fVI,y_fVI,z_fVI);

e_fI=max(max(max(abs(vt_fI-u_fI))))
e_fII=max(max(max(abs(vt_fII-u_fII))))
e_fIII=max(max(max(abs(vt_fIII-u_fIII))))
e_fIV=max(max(max(abs(vt_fIV-u_fIV))))
e_fV=max(max(max(abs(vt_fV-u_fV))))
e_fVI=max(max(max(abs(vt_fVI-u_fVI))))
