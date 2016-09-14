%% test rotational
clc; clear all; close all
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr
global teta0 teta1


teta0=-3*pi/16;
teta1=3*pi/16;
opt_ftr=0;
n=80;
mod72;

[vt_fI,vortvt_fI] = sol_exacte_vort(x_fI,y_fI,z_fI);
[vt_fII,vortvt_fII] = sol_exacte_vort(x_fII,y_fII,z_fII);
[vt_fIII,vortvt_fIII] = sol_exacte_vort(x_fIII,y_fIII,z_fIII);
[vt_fIV,vortvt_fIV] = sol_exacte_vort(x_fIV,y_fIV,z_fIV);
[vt_fV,vortvt_fV] = sol_exacte_vort(x_fV,y_fV,z_fV);
[vt_fVI,vortvt_fVI] = sol_exacte_vort(x_fVI,y_fVI,z_fVI);

[vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=...
    vort74(vt_fI, vt_fII, vt_fIII, vt_fIV, vt_fV, vt_fVI,n,nn);

err_I=abs(vortvt_fI-vort_fI);
err_II=abs(vortvt_fII-vort_fII);
err_III=abs(vortvt_fIII-vort_fIII);
err_IV=abs(vortvt_fIV-vort_fIV);
err_V=abs(vortvt_fV-vort_fV);
err_VI=abs(vortvt_fVI-vort_fVI);

e_I=max(max(err_I))
e_II=max(max(err_II))
e_III=max(max(err_III))
e_IV=max(max(err_IV))
e_V=max(max(err_V))
e_VI=max(max(err_VI))

err=max([e_I e_II e_III e_IV e_V e_VI])

figure(1)
plot_cs11(n,nn,err_I, err_II, err_III, err_IV, err_V, err_VI)
title('error')

figure(2)
plot_cs11(n,nn,vortvt_fI,vortvt_fII,vortvt_fIII,vortvt_fIV,vortvt_fV,vortvt_fVI)
title('exact')

figure(3)
plot_cs11(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
title('approx')

