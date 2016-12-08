% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
clc; close all; clear all;

global n nn;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global scheme

scheme='compact4';


n=63;
mod98; 

[funfI,funx,funy,funz]=fun3(x_fI,y_fI,z_fI);
[funfII,funx,funy,funz]=fun3(x_fII,y_fII,z_fII);
[funfIII,funx,funy,funz]=fun3(x_fIII,y_fIII,z_fIII);
[funfIV,funx,funy,funz]=fun3(x_fIV,y_fIV,z_fIV);
[funfV,funx,funy,funz]=fun3(x_fV,y_fV,z_fV);
[funfVI,funx,funy,funz]=fun3(x_fVI,y_fVI,z_fVI);

[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,va1,va2]=...
      gr100(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
  
  