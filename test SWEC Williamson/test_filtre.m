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

global n nn dxi
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr test scheme
global gp h0 u0 radius omega
global alpha

test=0;
opt_ftr=10;
scheme='compact4';

n=15; % for snapshot, n must be in the form 2^m-1 !
mod74

ccor=radius*omega;
cgrav=sqrt(h0*gp);
cvit=u0;
c=max([cgrav,ccor,cvit]);

cfl=0.7;
ddt=radius*dxi*cfl/c;
ndaymax=5;
Tmax=ndaymax*3600*24;
itermax=1;
comment='correction in filtering.';

tstart=cputime;
ref=floor(10000*now);
jour=date;
%% *** test data **********************************************************

if test == 0
    alpha=pi/2;
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

%% filtrage
[ hf_fI, hf_fII, hf_fIII, hf_fIV, hf_fV, hf_fVI, vf_fI, vf_fII, vf_fIII, vf_fIV, vf_fV, vf_fVI ]...
    = ftrcar74(ht_fI, ht_fII, ht_fIII, ht_fIV, ht_fV, ht_fVI, vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI);

err_fI=abs(hf_fI-ht_fI)./abs(ht_fI);
err_fII=abs(hf_fII-ht_fII)./abs(ht_fII);
err_fIII=abs(hf_fIII-ht_fIII)./abs(ht_fIII);
err_fIV=abs(hf_fIV-ht_fIV)./abs(ht_fIV);
err_fV=abs(hf_fV-ht_fV)./abs(ht_fV);
err_fVI=abs(hf_fVI-ht_fVI)./abs(ht_fVI);

figure(1)
plot_cs11(n,nn,ht_fI,ht_fII,ht_fIII,ht_fIV,ht_fV,ht_fVI);
title('initial function')

figure(2)
plot_cs11(n,nn,hf_fI, hf_fII, hf_fIII, hf_fIV, hf_fV, hf_fVI);
title('filtered function')

figure(3)
plot_cs11(n,nn,err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI);
title('error')

fig_placier