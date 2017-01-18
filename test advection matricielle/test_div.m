% MODULE PROBLEM FOR THE CUBED SPHERE
% ADVECTION EQUATION.
% ----------------------------------
% authors : Matthieu Brachet
%           Jean-Pierre Croisille
% ----------------------------------
clear all; clc; close all;
vvv=[1 1 1];
%% construction des variables globales
global n nn;
global radius u0 dxi;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI;
global ite aaa bbb itestop
global coef opt_ftr scheme
% test de Williamson
global alphad tetac lambdac
% test de Nair et Machenhauer
global gamma rho0 teta_p lambda_p
% test de Nair et Jablonowski
global teta0 lambda0
% test de Nair et Lauritzen
global lambdac1 tetac1 lambdac2 tetac2

time1=cputime;
%% *** OPTIONS ************************************************************
% si coef = 0, test 1 de Williamson (solid body rotation on the sphere)
%    coef = 1, test de Nair et Machenhauer  (deformational flow test - 
%                                                    stationnary vortex)
%    coef = 2, test de Nair, Jablonowski (moving vortices on the sphere)
%    coef = 3, test de Nair, Lauritzen (slotted cylinder) ( = Zaleska)
coef = 0;
% si film = 1 : faire le film,
%    film = 0 : ne pas faire.
film = 0;
% si save_graph = 1 : enregistrer les graphiques et les données dans TEST_SAVE.txt
%    save_graph = 0 : ne pas enregistrer
save_graph = 0;
% option de filtre : opt_ftr = ordre souhaité pour le filtre
% opt = 0 (sans filtre), 2, 4, 6, 8, 10
opt_ftr ='redonnet10';
% snapshot = 0 : pas de snapshot
%          = 1 : snapshot ( n must be odd. )
snapshot = 0;
% coupe = 0 : pas de coupe le long de l'équateur de la face 2
%         1 : coupe.
coupe = 0;
% sauvegarde = 1 : sauvegarde toutes les données,
%            = 0 : ne les sauvegarde pas, (utiliser load('namefile') pour
%            recharger les données).
sauvegarde = 0;
% choix du schéma aux différences finies
scheme='compact8'; % compact ou explicite
%% *** Benchmarks data ****************************************************
 n=31;
 nn=n+2;
 cfl=0.9;
 ndaymax=30;
 err=2;
 mm=0;
 MM=1000;
%% ************************************************************************
 if coef == 0
     %% test de Williamson
     alphad=3*pi/4;  
     lambdac=0;                                                           % longitude BUMP
     tetac=0;                                                                  % latitude BUMP
     lambda_p=pi;                                                              % position du pole nord, i.e. position du vortex nord
     teta_p=pi/2 - alphad;
 elseif coef == 1
     %% test de Nair et Machenhauer
     lambda_p=pi/4;                                                            % position du pole nord, i.e. position du vortex nord
     teta_p=-pi/4;
     rho0=3;
     gamma=5;
 elseif coef == 2
     %% test de Nair et Jablonowski
     alphad=pi/3; 
     lambda0 = 0;
     teta0 = 0;
     lambda_p=pi;                                                              % position du pole nord à t=0, i.e. position du vortex nord à t=0
     teta_p=pi/2 - alphad;
     rho0=3;
     gamma=5;
 elseif coef == 3
     %% test de Nair et Lauritzen
     alphad=3*pi/4;                                                                 % latitude BUMP
     lambda_p=pi;                                                              % position du pole nord, i.e. position du vortex nord
     teta_p=pi/2 - alphad;
     lambdac1=-pi/2;
     tetac1=0;
     lambdac2=pi/2;
     tetac2=0;
 end

%% données du problème
itestop=10000;
tstart=cputime;
mod101

ndays=6;
time=ndays*60*60*24;
[vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
        vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
        vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time);

vit_I(1:nn,1:nn,1)=vitx_I;
vit_I(1:nn,1:nn,2)=vity_I;
vit_I(1:nn,1:nn,3)=vitz_I;

vit_II(1:nn,1:nn,1)=vitx_II;
vit_II(1:nn,1:nn,2)=vity_II;
vit_II(1:nn,1:nn,3)=vitz_II;

vit_III(1:nn,1:nn,1)=vitx_III;
vit_III(1:nn,1:nn,2)=vity_III;
vit_III(1:nn,1:nn,3)=vitz_III;

vit_IV(1:nn,1:nn,1)=vitx_IV;
vit_IV(1:nn,1:nn,2)=vity_IV;
vit_IV(1:nn,1:nn,3)=vitz_IV;

vit_V(1:nn,1:nn,1)=vitx_V;
vit_V(1:nn,1:nn,2)=vity_V;
vit_V(1:nn,1:nn,3)=vitz_V;

vit_VI(1:nn,1:nn,1)=vitx_VI;
vit_VI(1:nn,1:nn,2)=vity_VI;
vit_VI(1:nn,1:nn,3)=vitz_VI;

[div_I,div_II,div_III,div_IV,div_V,div_VI]=...
    vort101(vit_I,vit_II,vit_III,vit_IV,vit_V,vit_VI,n,nn);

figure(1)
plot_cs11(n,nn,div_I,div_II,div_III,div_IV,div_V,div_VI);