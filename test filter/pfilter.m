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

[funfI,funfII,funfIII,funfIV,funfV,funfVI]=ftr_mixte101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI, n, nn);

errI   = abs(funfI-mfunfI);
errII  = abs(funfII-mfunfII);
errIII = abs(funfIII-mfunfIII);
errIV  = abs(funfIV-mfunfIV);
errV   = abs(funfV-mfunfV);
errVI  = abs(funfVI-mfunfVI);


[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'1');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'1');
e1=nrmg./MMM

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'2');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'2');
e2=nrmg./MMM

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(errI,errII,errIII,errIV,errV,errVI,n,nn,'infty');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MMM]=nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,'infty');
einf=nrmg./MMM
% 
% hFig = figure(1);
% set(gcf,'PaperPositionMode','auto')
% set(hFig, 'Position', [50 50 1000 500])
% plot_cs102(n,nn,mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI)
% %title('initial map')
% colorbar
% 
% hFig = figure(2);
% set(gcf,'PaperPositionMode','auto')
% set(hFig, 'Position', [50 50 1000 500])
% plot_cs102(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI)
% %title('filtered map')
% colorbar
% 
% hFig = figure(3);
% set(gcf,'PaperPositionMode','auto')
% set(hFig, 'Position', [50 50 1000 500])
% plot_cs102(n,nn,errI./MMM,errII./MMM,errIII./MMM,errIV./MMM,errV./MMM,errVI./MMM)
% %title('error')
% colorbar
