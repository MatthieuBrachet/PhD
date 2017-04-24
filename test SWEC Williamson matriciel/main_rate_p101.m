clc; clear all; close all;
global n
global opt_ftr test scheme filtre
global h0 u0
% *************************************************************************
% *** options :
% test = 0  : test 2 of Williamson & al.,
%        1  : test 5 of Williamson & al..
%        2  : test 5 of Williamson with smooth mountain,
%        3  : stationnary Galewsky (exp),
%        4  : Galewsky with perturbation (exp),
%        5  : Rossby-Haurwitz waves,
%        6  : Polar rotating low-high,
%        -1 : test with Earth topography
%        -2 : bump
% *************************************************************************
test=5;
filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';

ddt=150;
ndaymax=14;
n_init=31;

%% iterations
n=n_init;
[t1,ht_fI1,ht_fII1,ht_fIII1,ht_fIV1,ht_fV1,ht_fVI1,vt_fI1,vt_fII1,vt_fIII1,vt_fIV1,vt_fV1,vt_fVI1] = main_p101(n,ddt,ndaymax);
n=2*n_init+1;
[t2,ht_fI2,ht_fII2,ht_fIII2,ht_fIV2,ht_fV2,ht_fVI2,vt_fI2,vt_fII2,vt_fIII2,vt_fIV2,vt_fV2,vt_fVI2] = main_p101(n,ddt,ndaymax);
n=2*(2*n_init+1)+1;
[t3,ht_fI3,ht_fII3,ht_fIII3,ht_fIV3,ht_fV3,ht_fVI3,vt_fI3,vt_fII3,vt_fIII3,vt_fIV3,vt_fV3,vt_fVI3] = main_p101(n,ddt,ndaymax);

%% erreur
err1h_fI=ht_fI1-ht_fI2(1:2:end,1:2:end);
err1h_fII=ht_fII1-ht_fII2(1:2:end,1:2:end);
err1h_fIII=ht_fIII1-ht_fIII2(1:2:end,1:2:end);
err1h_fIV=ht_fIV1-ht_fIV2(1:2:end,1:2:end);
err1h_fV=ht_fV1-ht_fV2(1:2:end,1:2:end);
err1h_fVI=ht_fVI1-ht_fVI2(1:2:end,1:2:end);
errh1=max(max(abs([err1h_fI, err1h_fII, err1h_fIII, err1h_fIV, err1h_fV, err1h_fVI])))./h0;

err1v_fI=vt_fI1(:,:,3)-vt_fI2(1:2:end,1:2:end,3);
err1v_fII=vt_fII1(:,:,3)-vt_fII2(1:2:end,1:2:end,3);
err1v_fIII=vt_fIII1(:,:,3)-vt_fIII2(1:2:end,1:2:end,3);
err1v_fIV=vt_fIV1(:,:,3)-vt_fIV2(1:2:end,1:2:end,3);
err1v_fV=vt_fV1(:,:,3)-vt_fV2(1:2:end,1:2:end,3);
err1v_fVI=vt_fVI1(:,:,3)-vt_fVI2(1:2:end,1:2:end,3);
errv1=max(max(max(abs([err1v_fI, err1v_fII, err1v_fIII, err1v_fIV, err1v_fV, err1v_fVI]))))./u0;

err2h_fI=ht_fI2-ht_fI3(1:2:end,1:2:end);
err2h_fII=ht_fII2-ht_fII3(1:2:end,1:2:end);
err2h_fIII=ht_fIII2-ht_fIII3(1:2:end,1:2:end);
err2h_fIV=ht_fIV2-ht_fIV3(1:2:end,1:2:end);
err2h_fV=ht_fV2-ht_fV3(1:2:end,1:2:end);
err2h_fVI=ht_fVI2-ht_fVI3(1:2:end,1:2:end);
errh2=max(max(abs([err2h_fI, err2h_fII, err2h_fIII, err2h_fIV, err2h_fV, err2h_fVI])))./h0;

err2v_fI=vt_fI2(:,:,3)-vt_fI3(1:2:end,1:2:end,3);
err2v_fII=vt_fII2(:,:,3)-vt_fII3(1:2:end,1:2:end,3);
err2v_fIII=vt_fIII2(:,:,3)-vt_fIII3(1:2:end,1:2:end,3);
err2v_fIV=vt_fIV2(:,:,3)-vt_fIV3(1:2:end,1:2:end,3);
err2v_fV=vt_fV2(:,:,3)-vt_fV3(1:2:end,1:2:end,3);
err2v_fVI=vt_fVI2(:,:,3)-vt_fVI3(1:2:end,1:2:end,3);
errv2=max(max(max(abs([err2v_fI, err2v_fII, err2v_fIII, err2v_fIV, err2v_fV, err2v_fVI]))))./u0;

%% convergence
rate_h=log(errh1./errh2)/log(2)
rate_v=log(errv1./errv2)/log(2)

hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(n_init,err1h_fI./h0,err1h_fII./h0,err1h_fIII./h0,err1h_fIV./h0,err1h_fV./h0,err1h_fVI./h0)
title('Relative error coarse')
colorbar

hFig = figure(2);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(2*n_init+1,err2h_fI./h0,err2h_fII./h0,err2h_fIII./h0,err2h_fIV./h0,err2h_fV./h0,err2h_fVI./h0)
title('Relative error fine')
colorbar

hFig = figure(3);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(n_init,log2(abs(err1h_fI./err2h_fI(1:2:end,1:2:end))),log2(abs(err1h_fII./err2h_fII(1:2:end,1:2:end))),log2(abs(err1h_fIII./err2h_fIII(1:2:end,1:2:end)))...
    ,log2(abs(err1h_fIV./err2h_fIV(1:2:end,1:2:end))),log2(abs(err1h_fV./err2h_fV(1:2:end,1:2:end))),log2(abs(err1h_fVI./err2h_fVI(1:2:end,1:2:end))))
title('Rate of error')
colorbar

hFig = figure(4);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(n_init,ht_fI1,ht_fII1,ht_fIII1,ht_fIV1,ht_fV1,ht_fVI1)
title('Very coarse solution at final time')
colorbar

hFig = figure(5);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(2*n_init+1,ht_fI2,ht_fII2,ht_fIII2,ht_fIV2,ht_fV2,ht_fVI2)
title('Coarse solution at final time')
colorbar

hFig = figure(6);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(2*(2*n_init+1)+1,ht_fI3,ht_fII3,ht_fIII3,ht_fIV3,ht_fV3,ht_fVI3)
title('Fine solution at final time')
colorbar

fig_placier
