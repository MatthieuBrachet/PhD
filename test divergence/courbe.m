clear all;clc; close all;
global n nn;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global test;
%%%%%%%%%%%%%%%%%%%
test=4;
n=40;
mod67;

[mfunfI,dfunI]=fun6(x_fI,y_fI,z_fI); 
[mfunfII,dfunII]=fun6(x_fII,y_fII,z_fII);
[mfunfIII,dfunIII]=fun6(x_fIII,y_fIII,z_fIII);
[mfunfIV,dfunIV]=fun6(x_fIV,y_fIV,z_fIV);
[mfunfV,dfunV]=fun6(x_fV,y_fV,z_fV);
[mfunfVI,dfunVI]=fun6(x_fVI,y_fVI,z_fVI);

ppp=1
figure(1)
plot_cs5(n,nn,mfunfI(:,:,ppp),mfunfII(:,:,ppp),mfunfIII(:,:,ppp),mfunfIV(:,:,ppp),mfunfV(:,:,ppp),mfunfVI(:,:,ppp))

test=2;
[mfunfI,dfunI]=fun6(x_fI,y_fI,z_fI); 
[mfunfII,dfunII]=fun6(x_fII,y_fII,z_fII);
[mfunfIII,dfunIII]=fun6(x_fIII,y_fIII,z_fIII);
[mfunfIV,dfunIV]=fun6(x_fIV,y_fIV,z_fIV);
[mfunfV,dfunV]=fun6(x_fV,y_fV,z_fV);
[mfunfVI,dfunVI]=fun6(x_fVI,y_fVI,z_fVI);

ppp=3
figure(2)
plot_cs5(n,nn,mfunfI(:,:,ppp),mfunfII(:,:,ppp),mfunfIII(:,:,ppp),mfunfIV(:,:,ppp),mfunfV(:,:,ppp),mfunfVI(:,:,ppp))