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

[mfunfI,vorte_fI] = func(x_fI,y_fI,z_fI);
[mfunfII,vorte_fII] = func(x_fII,y_fII,z_fII);
[mfunfIII,vorte_fIII] = func(x_fIII,y_fIII,z_fIII);
[mfunfIV,vorte_fIV] = func(x_fIV,y_fIV,z_fIV);
[mfunfV,vorte_fV] = func(x_fV,y_fV,z_fV);
[mfunfVI,vorte_fVI] = func(x_fVI,y_fVI,z_fVI);

[vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI]=vort101(mfunfI,mfunfII,mfunfIII,mfunfIV,mfunfV,mfunfVI,n,nn);

%% calcul de l'erreur
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MM]=nrm101(vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI,n,nn,'infty');
err_fI=abs(vorte_fI-vort_fI);
err_fII=abs(vorte_fII-vort_fII);
err_fIII=abs(vorte_fIII-vort_fIII);
err_fIV=abs(vorte_fIV-vort_fIV);
err_fV=abs(vorte_fV-vort_fV);
err_fVI=abs(vorte_fVI-vort_fVI);

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,'1');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MM]=nrm101(vorte_fI,vorte_fII,vorte_fIII,vorte_fIV,vorte_fV,vorte_fVI,n,nn,'1');
nrmg./MM

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,'2');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MM]=nrm101(vorte_fI,vorte_fII,vorte_fIII,vorte_fIV,vorte_fV,vorte_fVI,n,nn,'2');
nrmg./MM

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,'infty');
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,MM]=nrm101(vorte_fI,vorte_fII,vorte_fIII,vorte_fIV,vorte_fV,vorte_fVI,n,nn,'infty');
nrmg./MM

% figure(1)
% plot_cs100(n,nn,vort_fI,vort_fII,vort_fIII,vort_fIV,vort_fV,vort_fVI)
% title('vort. numerique')
% colorbar
% 
% figure(2)
% plot_cs100(n,nn,vorte_fI,vorte_fII,vorte_fIII,vorte_fIV,vorte_fV,vorte_fVI)
% title('vort. exacte')
% colorbar
% 
% figure(3)
% plot_cs100(n,nn,err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI)
% title('vort. app.')
% colorbar
% 
% fig_placier