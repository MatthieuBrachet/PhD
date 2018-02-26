clc; clear all; close all; format shorte
global n nn
global scheme
global opt_ftr

scheme='compact4';
opt_ftr='redonnet10';

n=31;
mod101
disp('mod101 : ok')

[mfunfI]   = fun(x_fI,y_fI,z_fI);
[mfunfII]  = fun(x_fII,y_fII,z_fII);
[mfunfIII] = fun(x_fIII,y_fIII,z_fIII);
[mfunfIV]  = fun(x_fIV,y_fIV,z_fIV);
[mfunfV]   = fun(x_fV,y_fV,z_fV);
[mfunfVI]  = fun(x_fVI,y_fVI,z_fVI);

[funfI,funfII,funfIII,funfIV,funfV,funfVI]=ftr_xi101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI, n, nn);
figure(1)
plot_cs102(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI)
title('initial map')
colorbar

[funfI,funfII,funfIII,funfIV,funfV,funfVI]=ftr_eta101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI, n, nn);
figure(2)
plot_cs102(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI)
title('initial map')
colorbar