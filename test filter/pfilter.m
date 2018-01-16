clc; clear all; close all; format shorte
global n nn
global scheme
global opt_ftr

scheme='compact4';
opt_ftr='bogey6';

n=31;
mod101
disp('mod101 : ok')

[mfunfI]   = fun(x_fI,y_fI,z_fI);
[mfunfII]  = fun(x_fII,y_fII,z_fII);
[mfunfIII] = fun(x_fIII,y_fIII,z_fIII);
[mfunfIV]  = fun(x_fIV,y_fIV,z_fIV);
[mfunfV]   = fun(x_fV,y_fV,z_fV);
[mfunfVI]  = fun(x_fVI,y_fVI,z_fVI);

[funfI,funfII,funfIII,funfIV,funfV,funfVI]=ftr_mixte101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI, n, nn);
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'infty');

errI   = abs(funfI-mfunfI)./MMM;
errII  = abs(funfII-mfunfII)./MMM;
errIII = abs(funfIII-mfunfIII)./MMM;
errIV  = abs(funfIV-mfunfIV)./MMM;
errV   = abs(funfV-mfunfV)./MMM;
errVI  = abs(funfVI-mfunfVI)./MMM;


[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'1');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'1');
e1=nrmg./MMM

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'2');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'2');
e2=nrmg./MMM

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'infty');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'infty');
einf=nrmg./MMM

figure(1)
plot_cs102(n,nn,mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI)
title('initial map')
colorbar

figure(2)
plot_cs102(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI)
title('filtered map')
colorbar

figure(3)
plot_cs102(n,nn,errI,errII,errIII,errIV,errV,errVI)
title('error')
colorbar

fig_placier