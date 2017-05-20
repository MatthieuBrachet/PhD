% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
% authors : Matthieu Brachet
%           Jean-Pierre Croisille
% ----------------------------------
clear all; clc; close all;
%% construction des variables globales
global n nn;
global radius u0 dxi;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI;
global ite itestop
global coef opt_ftr scheme
global alphad tetac lambdac
global gamma rho0 teta_p lambda_p
global teta0 lambda0
global lambdac1 tetac1 lambdac2 tetac2
%% *** OPTIONS ************************************************************
coef = 2;
opt_ftr = 'redonnet10';
scheme='compact4';
coupe = 1;
sauvegarde = 0;
save_graph = 0;
%% *** Benchmarks data ****************************************************
 Na=29;
 Nb=59;
 ndaymax=12;
%% ************************************************************************
 if coef == 0
 %% test de Williamson
 alphad=0;  
 lambdac=0;                                                           % longitude BUMP
 tetac=0;                                                                  % latitude BUMP
 lambda_p=pi;                                                              % position du pole nord, i.e. position du vortex nord
 teta_p=pi/2 - alphad;
 elseif coef == 1
 %% test de Nair et Machenhauer
 lambda_p=0;                                                            % position du pole nord, i.e. position du vortex nord
 teta_p=0;
 rho0=3;
 gamma=5;
 elseif coef == 2
 %% test de Nair et Jablonowski
 alphad=pi/4; 
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
itestop=1;
tstart=cputime;


[xa,fa] = iterations_coupe(Na,ndaymax);
[xb,fb] = iterations_coupe(Nb,ndaymax);


%% graphiques
time=12*24*3600;
if sauvegarde == 1
    mkdir(['./filtre-' date ])
    save(['./filtre-' date '/ref_' num2str(ref) '_erreurdata_test_' num2str(coef) '.mat'],'funfI','funfII','funfIII','funfIV','funfV','funfVI')
end

if coupe == 1
    n=200;
    nn=n+2;
    mod101
    funfIe=fun4_b(x_fI,y_fI,z_fI,time);
    funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
    funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
    funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
    [ xe,fe ] = coupe_eq(funfIe,funfIIe,funfIIIe,funfIVe);
    figure(10)
    plot(xa,fa,'ko',xb,fb,'kx',xe,fe,'k-')
    grid on;
    legend(['approximate solution with N=' num2str(Na+1)],['approximate solution with N=' num2str(Nb+1)],'exact solution')
    xlabel('path Front')
    if save_graph==1
        mkdir(['./filtre-' date '/']);
        print('-dpng', ['./filtre-' date '/ref_' num2str(ref) '_coupefaceI_equateur_test_' num2str(coef) '.png'])
        savefig(['./filtre-' date '/ref_' num2str(ref) '_coupefaceI_test_' num2str(coef)]);
    end
end