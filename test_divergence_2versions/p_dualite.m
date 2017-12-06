clc; clear all; close all;
format shorte
global radius
global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global scheme nrm

scheme='compact6';
nrm='int';

n=255;
mod101
disp('mod101 : ok')

[h_fI,u_fI] = f_dual( x_fI,y_fI,z_fI );
[h_fII,u_fII] = f_dual( x_fII,y_fII,z_fII );
[h_fIII,u_fIII] = f_dual( x_fIII,y_fIII,z_fIII );
[h_fIV,u_fIV] = f_dual( x_fIV,y_fIV,z_fIV );
[h_fV,u_fV] = f_dual( x_fV,y_fV,z_fV );
[h_fVI,u_fVI] = f_dual( x_fVI,y_fVI,z_fVI );

[div_I,div_II,div_III,div_IV,div_V,div_VI]=div103(u_fI,u_fII,u_fIII,u_fIV,u_fV,u_fVI,n,nn);
[gr_I,gr_II,gr_III,gr_IV,gr_V,gr_VI]=gr101(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn);

funI=div_I.*h_fI;
funII=div_II.*h_fII;
funIII=div_III.*h_fIII;
funIV=div_IV.*h_fIV;
funV=div_V.*h_fV;
funVI=div_VI.*h_fVI;

% hFig=figure(1);
% set(gcf,'PaperPositionMode','auto');
% set(hFig, 'Position', [50 50 1000 500]);
% plot_cs100(n,nn,funI,funII,funIII,funIV,funV,funVI);
% colorbar

fun2I=scal(gr_I,u_fI);
fun2II=scal(gr_II,u_fII);
fun2III=scal(gr_III,u_fIII);
fun2IV=scal(gr_IV,u_fIV);
fun2V=scal(gr_V,u_fV);
fun2VI=scal(gr_VI,u_fVI);

% hFig=figure(2);
% set(gcf,'PaperPositionMode','auto');
% set(hFig, 'Position', [50 50 1000 500]);
% plot_cs100(n,nn,fun2I,fun2II,fun2III,fun2IV,fun2V,fun2VI);
% colorbar

[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg1]=nrm101(funI,funII,funIII,funIV,funV,funVI,n,nn,nrm);
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg2]=nrm101(fun2I,fun2II,fun2III,fun2IV,fun2V,fun2VI,n,nn,nrm);

nrmg1/(4*pi*radius.^2)
nrmg2/(4*pi*radius^2)
(nrmg1+nrmg2)/(4*pi*radius.^2)

fig_placier