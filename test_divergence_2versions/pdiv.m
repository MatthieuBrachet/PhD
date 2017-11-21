clc; clear all; close all;
format shorte
global radius
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global scheme nrm

scheme='compact4';
nrm='int';

n=63;
mod101
disp('mod101 : ok')

[ funfI, div_fIe ] = func( x_fI,y_fI,z_fI );
[ funfII, div_fIIe ] = func( x_fII,y_fII,z_fII );
[ funfIII, div_fIIIe ] = func( x_fIII,y_fIII,z_fIII );
[ funfIV, div_fIVe ] = func( x_fIV,y_fIV,z_fIV );
[ funfV, div_fVe ] = func( x_fV,y_fV,z_fV );
[ funfVI, div_fVIe ] = func( x_fVI,y_fVI,z_fVI );

[div_I,div_II,div_III,div_IV,div_V,div_VI]=div101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

err_fI=div_I-div_fIe;
err_fII=div_II-div_fIIe;
err_fIII=div_III-div_fIIIe;
err_fIV=div_IV-div_fIVe;
err_fV=div_V-div_fVe;
err_fVI=div_VI-div_fVIe;

str='1';
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,M]=nrm101(div_fIe,div_fIIe,div_fIIIe,div_fIVe,div_fVe,div_fVIe,n,nn,str);
nrmg./M

str='2';
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,M]=nrm101(div_fIe,div_fIIe,div_fIIIe,div_fIVe,div_fVe,div_fVIe,n,nn,str);
nrmg./M

str='infty';
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,M]=nrm101(div_fIe,div_fIIe,div_fIIIe,div_fIVe,div_fVe,div_fVIe,n,nn,str);
nrmg./M

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=nrm101(div_I,div_II,div_III,div_IV,div_V,div_VI,n,nn,nrm);
surf_earth=4*pi*radius.^2;
nrmg./(M*surf_earth)

% hFig=figure(1)
% set(gcf,'PaperPositionMode','auto')
% set(hFig, 'Position', [50 50 1000 500])
% plot_cs100(n,nn,abs(err_fI)./M,abs(err_fII)./M,abs(err_fIII)./M,abs(err_fIV)./M,abs(err_fV)./M,abs(err_fVI)./M)
% colorbar