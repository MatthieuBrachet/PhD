clc; clear all; close all;
global n
global opt_ftr test scheme filtre
test=1;
filtre='classic';
opt_ftr='redonnet10';
scheme='compact4';

n_init=74;
%% iterations
ddt=100;
n=n_init;
ndaymax=15;
[t1,ht_fI1,ht_fII1,ht_fIII1,ht_fIV1,ht_fV1,ht_fVI1,vt_fI1,vt_fII1,vt_fIII1,vt_fIV1,vt_fV1,vt_fVI1] = main_p101(n,ddt,ndaymax);

load('sol_ref_test_1_N_149.mat');
nref=n;
htref_fI=ht_fI;
htref_fII=ht_fII;
htref_fIII=ht_fIII;
htref_fIV=ht_fIV;
htref_fV=ht_fV;
htref_fVI=ht_fVI;

vtref_fI=vt_fI;
vtref_fII=vt_fII;
vtref_fIII=vt_fIII;
vtref_fIV=vt_fIV;
vtref_fV=vt_fV;
vtref_fVI=vt_fVI;

%% calcul d'erreur relative
pas=(nref+1)/(n_init+1);

err1h_fI=ht_fI1-htref_fI(1:pas:end,1:pas:end);
err1h_fII=ht_fII1-htref_fII(1:pas:end,1:pas:end);
err1h_fIII=ht_fIII1-htref_fIII(1:pas:end,1:pas:end);
err1h_fIV=ht_fIV1-htref_fIV(1:pas:end,1:pas:end);
err1h_fV=ht_fV1-htref_fV(1:pas:end,1:pas:end);
err1h_fVI=ht_fVI1-htref_fVI(1:pas:end,1:pas:end);
errh1=max(max(abs([err1h_fI, err1h_fII, err1h_fIII, err1h_fIV, err1h_fV, err1h_fVI])))./max(max(abs([htref_fI,htref_fII,htref_fIII,htref_fIV,htref_fV,htref_fVI])));

err1v_fI=vt_fI1-vtref_fI(1:pas:end,1:pas:end,:);
err1v_fII=vt_fII1-vtref_fII(1:pas:end,1:pas:end,:);
err1v_fIII=vt_fIII1-vtref_fIII(1:pas:end,1:pas:end,:);
err1v_fIV=vt_fIV1-vtref_fIV(1:pas:end,1:pas:end,:);
err1v_fV=vt_fV1-vtref_fV(1:pas:end,1:pas:end,:);
err1v_fVI=vt_fVI1-vtref_fVI(1:pas:end,1:pas:end,:);
errv1=max(max(max(abs([err1v_fI, err1v_fII, err1v_fIII, err1v_fIV, err1v_fV, err1v_fVI]))))./max(max(max(abs([vtref_fI,vtref_fII,vtref_fIII,vtref_fIV,vtref_fV,vtref_fVI]))));

disp([n_init nref])
disp([errh1 errv1])

%% courbe
hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(n_init,err1h_fI./abs(htref_fI(1:pas:end,1:pas:end)),err1h_fII./abs(htref_fII(1:pas:end,1:pas:end)),err1h_fIII./abs(htref_fIII(1:pas:end,1:pas:end)),err1h_fIV./abs(htref_fIV(1:pas:end,1:pas:end)),err1h_fV./abs(htref_fV(1:pas:end,1:pas:end)),err1h_fVI./abs(htref_fVI(1:pas:end,1:pas:end)))
title('Relative error coarse')
colorbar

hFig = figure(2);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(n_init,ht_fI1,ht_fII1,ht_fIII1,ht_fIV1,ht_fV1,ht_fVI1)
title(['Solution at time ' num2str(t1/(24*3600)) ' days'])
colorbar

hFig = figure(3);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [50 50 1000 500])
plot_err(nref,htref_fI,htref_fII,htref_fIII,htref_fIV,htref_fV,htref_fVI)
title(['Reference solution at time ' num2str(t/(24*3600)) ' days'])
colorbar

fig_placier