clc; clear all; close all;
format long

global n nn
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI
global opt_ftr opt_ftr1 opt_detec test scheme

comment='.';
test=1;
video = 'no';
nper=2;
sauvegarde = 1;
filtre='classic';
opt_ftr='inf';
opt_detec='redonnet10';
opt_ftr1='redonnet4';

scheme='compact4';
snapshot='yes';

n=127; % for snapshot, n must be in the form 2^m-1 !
ndaymax=5;
mod74

[alt_fI]=topography(x_fI,y_fI,z_fI);
[alt_fII]=topography(x_fII,y_fII,z_fII);
[alt_fIII]=topography(x_fIII,y_fIII,z_fIII);
[alt_fIV]=topography(x_fIV,y_fIV,z_fIV);
[alt_fV]=topography(x_fV,y_fV,z_fV);
[alt_fVI]=topography(x_fVI,y_fVI,z_fVI);

figure(1)
plot_cs11(n,nn,alt_fI,alt_fII,alt_fIII,alt_fIV,alt_fV,alt_fVI);

figure(2)
plot_cs7(n,nn,alt_fI,alt_fII,alt_fIII,alt_fIV,alt_fV,alt_fVI);