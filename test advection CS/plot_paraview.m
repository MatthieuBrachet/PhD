clc; clear all; close all;

% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
% authors : Matthieu Brachet
%           Jean-Pierre Croisille
% ----------------------------------
% RK4 + Filtrage à l'ordre 10
clear all; clc; close all;
%% construction des variables globales
global n nn;
global radius;
global dxi;
global x_fI y_fI z_fI;
global x_fII y_fII z_fII;
global x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV;
global x_fV y_fV z_fV;
global x_fVI y_fVI z_fVI;
global ite aaa bbb itestop
global coef u0
% test de Williamson
global alphad
global tetac lambdac
% test de Nair et Machenhauer
global gamma rho0
global teta_p lambda_p
% test de Nair et Jablonowski
global teta0 lambda0
%% *** OPTIONS ************************************************************
%
% si coef = 0, test 1 de Williamson (solid body rotation on the sphere)
%    coef = 1, test de Nair et Machenhauer  (deformational flow test - 
%                                                    stationnary vortex)
%    coef = 2, test de Nair, Jablonowski (moving vortices on the sphere)
coef = 2;
% si film = 1 : faire le film,
%    film = 0 : ne pas faire.
film = 0;
% si qquiv = 1 : tracer le champ de vecteurs
%    qquiv = 0 : ne pas tracer
qquiv = 0;
% si save_graph = 1 : enregistrer les graphiques
%    save_graph = 0 : ne pas enregistrer
save_graph = 0;
%
%% *** Benchmarks data ****************************************************

%% 
 n=150;
 nn=n+2;
 nday=12;
 time=24*3600*nday;
 
  if coef == 0
 % test de Williamson
 alphad=pi/2;  
 lambdac=3*pi/2;                                                           % longitude BUMP
 tetac=0;                                                                  % latitude BUMP
 lambda_p=pi;                                                              % position du pole nord, i.e. position du vortex nord
 teta_p=pi/2 - alphad;

 elseif coef == 1
 % test de Nair et Machenhauer
 lambda_p=3*pi/4;                                                            % position du pole nord, i.e. position du vortex nord
 teta_p=-pi/4;
 rho0=3;
 gamma=5;
 
 elseif coef == 2
 % test de Nair et Jablonowski
 alphad=3*pi/4; 
 
 lambda0 = 3*pi/2;
 teta0 = 0;
 lambda_p=pi;                                                              % position du pole nord à t=0, i.e. position du vortex nord à t=0
 teta_p=pi/2 - alphad;
 rho0=3;
 gamma=5;
 
  end
 
  mod_1b
  disp('mod ok')
funfIe=fun4_b(x_fI,y_fI,z_fI,time);
funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
funfVe=fun4_b(x_fV,y_fV,z_fV,time);
funfVIe=fun4_b(x_fVI,y_fVI,z_fVI,time);
  
data_save(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe);
 
 
 
 
 
 
 
 
 