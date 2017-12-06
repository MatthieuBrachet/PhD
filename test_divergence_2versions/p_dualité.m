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

[ u_fI, div_fIe ] = func( x_fI,y_fI,z_fI );
[ u_fII, div_fIIe ] = func( x_fII,y_fII,z_fII );
[ u_fIII, div_fIIIe ] = func( x_fIII,y_fIII,z_fIII );
[ u_fIV, div_fIVe ] = func( x_fIV,y_fIV,z_fIV );
[ u_fV, div_fVe ] = func( x_fV,y_fV,z_fV );
[ u_fVI, div_fVIe ] = func( x_fVI,y_fVI,z_fVI );

[h_fI,gr] = fun2(x_fI,y_fI,z_fI);
[h_fII,gr] = fun2(x_fII,y_fII,z_fII);
[h_fIII,gr] = fun2(x_fIII,y_fIII,z_fIII);
[h_fIV,gr] = fun2(x_fIV,y_fIV,z_fIV);
[h_fV,gr] = fun2(x_fV,y_fV,z_fV);
[h_fVI,gr] = fun2(x_fVI,y_fVI,z_fVI);

[div_I,div_II,div_III,div_IV,div_V,div_VI]=div101(u_fI,u_fII,u_fIII,u_fIV,u_fV,u_fVI,n,nn);
[gr_I,gr_II,gr_III,gr_IV,gr_V,gr_VI]=gr101(h_fI,h_fII,h_fIII,h_fIV,h_fV,h_fVI,n,nn);

