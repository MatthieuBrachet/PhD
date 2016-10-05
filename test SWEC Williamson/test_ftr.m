%% test filtre
clc; clear all; close all;
format short
global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr alfa_ftr test scheme
global gp h0 u0 radius omega
global alpha

test=1;
video = 'no';
nper=1;
sauvegarde = 1;
opt_ftr='redonnet0';
alfa_ftr=0.45;
delta_ftr=0;
scheme='compact4';
snapshot='no';

n=63; % for snapshot, n must be in the form 2^m-1 !
mod74

ccor=radius*omega;
cgrav=sqrt(h0*gp);
cvit=u0;
c=max([cgrav,ccor,cvit]);

cfl=0.9;
ddt=radius*dxi*cfl/c;
ndaymax=5;
Tmax=ndaymax*3600*24;
itermax=5000;

tstart=cputime;
ref=floor(10000*now);
jour=date;
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


t=0;
[ ht_fI,    vt_fI]   = sol_exacte(x_fI,   y_fI,   z_fI,   t);
[ ht_fII,   vt_fII]  = sol_exacte(x_fII,  y_fII,  z_fII,  t);
[ ht_fIII,  vt_fIII] = sol_exacte(x_fIII, y_fIII, z_fIII, t);
[ ht_fIV,   vt_fIV]  = sol_exacte(x_fIV,  y_fIV,  z_fIV,  t);
[ ht_fV,    vt_fV]   = sol_exacte(x_fV,   y_fV,   z_fV,   t);
[ ht_fVI,   vt_fVI]  = sol_exacte(x_fVI,  y_fVI,  z_fVI,  t);


[htf_fI,htf_fII,htf_fIII,htf_fIV,htf_fV,htf_fVI]=...
            ftr72(ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI,n,nn);

e_fI=abs(htf_fI-ht_fI)./abs(ht_fI);
e_fII=abs(htf_fII-ht_fII)./abs(ht_fII);
e_fIII=abs(htf_fIII-ht_fIII)./abs(ht_fIII);
e_fIV=abs(htf_fIV-ht_fIV)./abs(ht_fIV);
e_fV=abs(htf_fV-ht_fV)./abs(ht_fV);
e_fVI=abs(htf_fVI-ht_fVI)./abs(ht_fVI);
        
  
figure(1)
plot_cs7(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title('original data')

figure(2)
plot_cs7(n,nn,htf_fI,htf_fII,htf_fIII,htf_fIV,htf_fV,htf_fVI);
title('filtred data')

figure(3)
plot_cs7(n,nn,e_fI, e_fII, e_fIII, e_fIV, e_fV, e_fVI);
title('relative error')

fig_placier








