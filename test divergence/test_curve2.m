clear all;clc; close all;
global n nn;
global radius
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global test
%%%%%%%%%%%%%%%%%%%

n=20;
test=4;
mod67;

[u_fI] = test_curve(x_fI,y_fI,z_fI);
[u_fII] = test_curve(x_fII,y_fII,z_fII);
[u_fIII] = test_curve(x_fIII,y_fIII,z_fIII);
[u_fIV] = test_curve(x_fIV,y_fIV,z_fIV);
[u_fV] = test_curve(x_fV,y_fV,z_fV);
[u_fVI] = test_curve(x_fVI,y_fVI,z_fVI);

figure(1)
plot_cs5(n,nn,0*u_fI,0*u_fII,0*u_fIII,0*u_fIV,0*u_fV,0*u_fVI); hold on

[funfI,dfunI]=fun6(x_fI,y_fI,z_fI); 
[funfII,dfunII]=fun6(x_fII,y_fII,z_fII);
[funfIII,dfunIII]=fun6(x_fIII,y_fIII,z_fIII);
[funfIV,dfunIV]=fun6(x_fIV,y_fIV,z_fIV);
[funfV,dfunV]=fun6(x_fV,y_fV,z_fV);
[funfVI,dfunVI]=fun6(x_fVI,y_fVI,z_fVI);

plot_quiver(nn,funfI,funfII, funfIII, funfIV, funfV, funfVI)