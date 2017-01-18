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


n=31;
mod98; 

[funfI,funx,funy,funz]=fun3(x_fI,y_fI,z_fI);
[funfII,funx,funy,funz]=fun3(x_fII,y_fII,z_fII);
[funfIII,funx,funy,funz]=fun3(x_fIII,y_fIII,z_fIII);
[funfIV,funx,funy,funz]=fun3(x_fIV,y_fIV,z_fIV);
[funfV,funx,funy,funz]=fun3(x_fV,y_fV,z_fV);
[funfVI,funx,funy,funz]=fun3(x_fVI,y_fVI,z_fVI);

[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,vad_fI]=...
      gr100(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
  
  
[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI,vad_fIt]=...
      gr100t(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);
  
  
figure(1)
surf(vad_fI-vad_fIt)
title('erreur globale')
% 
% figure(2)
% surf(vad_fI(1:nn-1,1:nn-1)-vad_fIt(1:nn-1,1:nn-1))
% title('erreur face A')
% 
% figure(3)
% surf(vad_fI(1:nn-1,nn:2*(nn-1))-vad_fIt(1:nn-1,nn:2*(nn-1)))
% title('erreur face B')
% 
% figure(4)
% surf(vad_fI(1:nn-1,2*nn-1:3*(nn-1))-vad_fIt(1:nn-1,2*nn-1:3*(nn-1)))
% title('erreur face C')
% 
% figure(5)
% surf(vad_fI(1:nn-1,3*(nn-1)+1:end)-vad_fIt(1:nn-1,3*(nn-1)+1:end))
% title('erreur face D')
%  

fig_placier;


