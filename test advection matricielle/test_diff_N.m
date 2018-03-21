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
global itestop
global coef opt_ftr scheme
global teta_p lambda_p rho0 gamma
global time
%% *** OPTIONS ************************************************************
scheme='compact4';
save_graph=0;
coupe=1;

%% *** Benchmarks data ****************************************************
 ndaymax=12;

 coef=1;
 lambda_p=3*pi/4;
 teta_p=0;
 rho0=3;
 gamma=5;

%% données du problème
itestop=1;
tstart=cputime;


opt_ftr = 'redonnet10';
[xa,fa] = iterations_coupe(24,ndaymax);

[xb,fb] = iterations_coupe(49,ndaymax);



%% graphiques
n=200;
nn=n+2;
mod101
funfIe=fun4_b(x_fI,y_fI,z_fI,time);
funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
[ xex,fex ] = coupe_eq(funfIe,funfIIe,funfIIIe,funfIVe);


figure(10)
plot(xex,fex,'k-',xa,fa,'rx',xb,fb,'bo','Linewidth',2)
grid minor;
legend('Exacte','N=25','N=50','Location','NorthWest')
xticks([3*pi/2 7*pi/4 2*pi])
xticklabels({'3\pi/2','7\pi/4','2\pi'})
axis([3*pi/2-.05 2*pi+.05 0.4 1.6])