clc; clear all; close all;
format shorte
global radius
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global scheme nrm

scheme='compact4';
nrm='int';

n=255;
mod101
disp('mod101 : ok')

[mfunfI] = fun(x_fI,y_fI,z_fI);
[mfunfII] = fun(x_fII,y_fII,z_fII);
[mfunfIII] = fun(x_fIII,y_fIII,z_fIII);
[mfunfIV] = fun(x_fIV,y_fIV,z_fIV);
[mfunfV] = fun(x_fV,y_fV,z_fV);
[mfunfVI] = fun(x_fVI,y_fVI,z_fVI);
Esurf=4*pi*radius^2;

[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=gr101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);
[vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=vort101(grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,n,nn);

% figure(1)
% plot_cs100(n,nn,mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI)
% title('h')
% colorbar
% 
% 
% figure(2)
% plot_cs100(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
% title('vort. numerique')
% colorbar


[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI,n,nn,'1');
nrmg
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI,n,nn,'2');
nrmg
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI,n,nn,'infty');
nrmg