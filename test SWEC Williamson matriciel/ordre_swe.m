%% ************************************************************************
% Solve SWEC on the Cubed-Sphere.
% *** options :
% test = 0  : test 2 of Williamson & al.,
%        1  : test 5 of Williamson & al..
%        2  : test 5 of Williamson with smooth mountain,
%        3  : stationnary Galewsky (exp),
%        4  : Galewsky with perturbation (exp),
%        5  : Rossby-Haurwitz waves,
%        -1 : test with Earth topography
%        -2 : bump
% scheme : numerical spatial scheme used. 
% video  : 'yes' ou 'no', do a video or not,
% nper   :  periodicity of frames in the video.
% sauvegarde = 0 (do not save data), 1 (save all data).
% opt_ftr : explicit (redonnet).
%% ************************************************************************
clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr filtre test scheme
global gp h0 u0 radius omega

comment=' .';
test=5;
video = 'no';
sauvegarde = 1;
filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';
snapshot='yes';

n1=31;
n2=63;

radius=6.37122d+06;
omega=7.292d-05;
gp=9.80616;
ddx1=2*pi*radius./(4.*(n1+1));
ddx2=2*pi*radius./(4.*(n2+1));

ccor=radius*omega;
cgrav=sqrt(h0*gp);
cvit=u0;
c=max([cgrav,ccor,cvit]);
cfl=0.9;
ddtime=radius*cfl/c*min(ddx1,ddx2);

nday=5;
Tmax=nday*ddtime%*3600*24;

disp('part 1...')
[ht1_fI, ht1_fII, ht1_fIII, ht1_fIV, ht1_fV, ht1_fVI,...
    vt1_fI, vt1_fII, vt1_fIII, vt1_fIV, vt1_fV, vt1_fVI] = iter_swe( n1, ddtime , Tmax);

disp('part 2...')
[ht2_fI, ht2_fII, ht2_fIII, ht2_fIV, ht2_fV, ht2_fVI,...
    vt2_fI, vt2_fII, vt2_fIII, vt2_fIV, vt2_fV, vt2_fVI] = iter_swe( n2, ddtime , Tmax);

err_fI=ht1_fI-ht2_fI(1:2:end,1:2:end);
err_fII=ht1_fII-ht2_fII(1:2:end,1:2:end);
err_fIII=ht1_fIII-ht2_fIII(1:2:end,1:2:end);
err_fIV=ht1_fIV-ht2_fIV(1:2:end,1:2:end);
err_fV=ht1_fV-ht2_fV(1:2:end,1:2:end);
err_fVI=ht1_fVI-ht2_fVI(1:2:end,1:2:end);

nrmI=max(max(abs(err_fI)));
nrmII=max(max(abs(err_fII)));
nrmIII=max(max(abs(err_fIII)));
nrmIV=max(max(abs(err_fIV)));
nrmV=max(max(abs(err_fV)));
nrmVI=max(max(abs(err_fVI)));
err=max([nrmI nrmII nrmIII nrmIV nrmV nrmVI]);

err