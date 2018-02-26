clc; clear all; close all; format shorte
global n nn
global scheme
global opt_ftr

scheme='compact4';
opt_ftr='redonnet10';

n=15;
mod101
disp('mod101 : ok')

[mfunfI]   = fun(x_fI,y_fI,z_fI);
[mfunfII]  = fun(x_fII,y_fII,z_fII);
[mfunfIII] = fun(x_fIII,y_fIII,z_fIII);
[mfunfIV]  = fun(x_fIV,y_fIV,z_fIV);
[mfunfV]   = fun(x_fV,y_fV,z_fV);
[mfunfVI]  = fun(x_fVI,y_fVI,z_fVI);

[funfI,funfII,funfIII,funfIV,funfV,funfVI]=ftr_xi101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI, n, nn);
[funfI,funfII,funfIII,funfIV,funfV,funfVI]=ftr_eta101(funfI,funfII,funfIII,funfIV,funfV,funfVI, n, nn);
figure(1)
plot_cs102(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI)
title('initial map')
colorbar

[funfI2,funfII2,funfIII2,funfIV2,funfV2,funfVI2]=ftr_eta101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI, n, nn);
[funfI2,funfII2,funfIII2,funfIV2,funfV2,funfVI2]=ftr_xi101(funfI2,funfII2,funfIII2,funfIV2,funfV2,funfVI2, n, nn);
figure(2)
plot_cs102(n,nn,funfI2,funfII2,funfIII2,funfIV2,funfV2,funfVI2)
title('initial map')
colorbar

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'infty');
MMM=1;
errI   = abs(funfI-funfI2);
errII  = abs(funfII-funfII2);
errIII = abs(funfIII-funfIII2);
errIV  = abs(funfIV-funfIV2);
errV   = abs(funfV-funfV2);
errVI  = abs(funfVI-funfVI2);

figure(3)
plot_cs11(n,nn,errI,errII,errIII,errIV,errV,errVI)
colorbar