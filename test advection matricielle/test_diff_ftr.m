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
global alphad tetac lambdac
global teta_p lambda_p
global time
%% *** OPTIONS ************************************************************
scheme='compact4';
save_graph=0;
coupe=1;

%% *** Benchmarks data ****************************************************
 na=79;
 ndaymax=12;

 coef=0;
 alphad=pi/4;  
 lambdac=3*pi/2;
 tetac=0;
 lambda_p=pi;
 teta_p=pi/2 - alphad;

%% données du problème
itestop=1;
tstart=cputime;


opt_ftr = 'redonnet10';
[xa,fa] = iterations_coupe(na,ndaymax);

opt_ftr = 'redonnet8';
[xb,fb] = iterations_coupe(na,ndaymax);

opt_ftr = 'redonnet6';
[xc,fc] = iterations_coupe(na,ndaymax);

opt_ftr = 'redonnet4';
[xd,fd] = iterations_coupe(na,ndaymax);

opt_ftr = 'redonnet2';
[xe,fe] = iterations_coupe(na,ndaymax);


%% graphiques
if coupe == 1
    n=200;
    nn=n+2;
    mod101
    funfIe=fun4_b(x_fI,y_fI,z_fI,time);
    funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
    funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
    funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
    [ xex,fex ] = coupe_eq(funfIe,funfIIe,funfIIIe,funfIVe);
    
    
    figure(10)
    plot(xa,fa,xb,fb,xc,fc,xd,fd,xe,fe,xex,fex,'Linewidth',2)
    grid minor;
    legend('filtre 10','filtre 8','filtre 6','filtre 4','filtre 2','Solution exacte','Location','NorthWest')
    xticks([0 pi/4 pi/2 3*pi/4 pi 5*pi/4 3*pi/2 7*pi/4 2*pi])
    xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'})
    axis([pi-.05 2*pi+.05 -50 1050])
    
    if save_graph==1
        mkdir(['./filtre-' date '/']);
        print('-dpng', ['./filtre-' date '/coupefaceI_equateur_test_' num2str(coef) '_N20.png'])
        savefig(['./filtre-' date '/coupefaceI_test_' num2str(coef)]);
    end
end