% MODULE PROBLEM FOR THE CUBED SPHERE
% ADVECTION EQUATION.
% ----------------------------------
% authors : Matthieu Brachet
%           Jean-Pierre Croisille
% ----------------------------------
clear all; clc; close all; format shorte
vvv=[1 1 1];
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

time1=cputime;
%% *** OPTIONS ************************************************************
% si coef = 0, test 1 de Williamson (solid body rotation on the sphere)
%    coef = 1, test de Nair et Machenhauer  (deformational flow test - 
%                                                    stationnary vortex)
%    coef = 2, test de Nair, Jablonowski (moving vortices on the sphere)
%    coef = 3, test de Nair, Lauritzen (slotted cylinder) ( = Zaleska)
coef = 0;
save_graph = 1;
opt_ftr ='redonnet10';
coupe = 0;
sauvegarde = 1;
scheme='compact4';
%% *** Benchmarks data ****************************************************
 n=31;
 nn=n+2;
 cfl=.7;
 ndaymax=12;
 err=2;
 mm=0;
 MM=1000;
%% ************************************************************************
 if coef == 0
     %% test 1 de Williamson
     alphad=pi/4;  
     lambdac=3*pi/2;                             % longitude BUMP
     tetac=0;                                    % latitude BUMP
     lambda_p=pi;                                % position du pole nord, i.e. position du vortex nord
     teta_p=pi/2 - alphad;
 elseif coef == 1
     %% test de Nair et Machenhauer
     lambda_p=pi/4;                              % position du pole nord, i.e. position du vortex nord
     teta_p=pi/4;
     rho0=3;
     gamma=5;
 elseif coef == 2
     %% test de Nair et Jablonowski
     alphad=pi/4; 
     lambda0 = 3*pi/2;
     teta0 = 0;
     lambda_p=pi;                                % position du pole nord à t=0, i.e. position du vortex nord à t=0
     teta_p=pi/2 - alphad;
     rho0=3;
     gamma=5;
 elseif coef == 3
     %% test de Nair et Lauritzen
     alphad=3*pi/4;                              % latitude BUMP
     lambda_p=pi;                                % position du pole nord, i.e. position du vortex nord
     teta_p=pi/2 - alphad;
     lambdac1=-pi/2;
     tetac1=0;
     lambdac2=pi/2;
     tetac2=0;
 end

%% données du problème
itestop=2000;
tstart=cputime;
mod101
tmax=24*3600*ndaymax;
ddt=cfl*radius*dxi/u0;
itemax=floor(tmax/ddt);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUT BOUCLE EN TEMPS %
% ------------------------------------------------------------------------
% 5 - CALCUL DERIVEES HERMITIENNES SUR LES 6 RESEAUX DE GRANDS CERCLES
% *** PHASE 1 CHARGEMENT DES DONNEES, CELA COMPREND 
%      A- LE TRANSFERT SIMPLES DES VALEURS SUR DEUX FACES 
%      B- L'INTERPOLATION SPLINE CUBIQUE SUR LES DEUX AUTRES FACES (DE TYPE
%                                                                   CROSS).
% *** PHASE 2: INTERPOLATION DE LA DERIVEE HERMITIENNE.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -----------------------------------------------------------------------
% ETAPE 1: CALCUL DE LA VITESSE DE CONVECTION SOLIDE SUR CHAQUE FACE.
% ------------------------------------------------------------------------
vitx_I=zeros(nn,nn);vity_I=zeros(nn,nn);vitz_I=zeros(nn,nn);
vitx_II=zeros(nn,nn);vity_II=zeros(nn,nn);vitz_II=zeros(nn,nn);
vitx_III=zeros(nn,nn);vity_III=zeros(nn,nn);vitz_III=zeros(nn,nn);
vitx_IV=zeros(nn,nn);vity_IV=zeros(nn,nn);vitz_IV=zeros(nn,nn);
vitx_V=zeros(nn,nn);vity_V=zeros(nn,nn);vitz_V=zeros(nn,nn);
vitx_VI=zeros(nn,nn);vity_VI=zeros(nn,nn);vitz_VI=zeros(nn,nn);

% appel une fois fictive de fun4 pour chargement de u0=module vitesse  et
% alphad=angle du plan de propagation et axe polaire.
% fun4 = solution exacte
%[vwk0]=fun4_b(0,0,0,0);
lambda=zeros(nn,nn);
teta=zeros(nn,nn);
radius1=zeros(nn,nn);
elambda_x=zeros(nn,nn);
elambda_y=zeros(nn,nn);
elambda_z=zeros(nn,nn);
eteta_x=zeros(nn,nn);
eteta_y=zeros(nn,nn);
eteta_z=zeros(nn,nn);

%% données pour RK4
tinit=0;
% initial condition
[funfI]=fun4_b(x_fI,y_fI,z_fI,tinit);
[funfII]=fun4_b(x_fII,y_fII,z_fII,tinit);
[funfIII]=fun4_b(x_fIII,y_fIII,z_fIII,tinit);
[funfIV]=fun4_b(x_fIV,y_fIV,z_fIV,tinit);
[funfV]=fun4_b(x_fV,y_fV,z_fV,tinit);
[funfVI]=fun4_b(x_fVI,y_fVI,z_fVI,tinit);

%% masse de reference
str='int';
[nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
    nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,str);
mass_ref=nrmg;


time=tinit;



%% Boucles RK 4 avec filtrage
xdays(1)=0;
ite =1;
erinfty(1)=0; er2(1)=0; er1(1)=0;
ref=floor(10000*now);
while ite<itemax & erinfty(end)<10
    clc; disp(num2str([ite itemax er1(end) er2(end) erinfty(end)]));

    % Filtrage
    [funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_mixte101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

    funfI=funftI;funfII=funftII;funfIII=funftIII;
    funfIV=funftIV;funfV=funftV;funfVI=funftVI;

    % iterations
    %% CALCUL KK0
    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

    [vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
        vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
        vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time);

    kk0_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
    kk0_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
    kk0_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
    kk0_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
    kk0_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
    kk0_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 

    %% CALCUL KK1
    fun1_fI=funfI+0.5*ddt*kk0_I;
    fun1_fII=funfII+0.5*ddt*kk0_II;
    fun1_fIII=funfIII+0.5*ddt*kk0_III;
    fun1_fIV=funfIV+0.5*ddt*kk0_IV;
    fun1_fV=funfV+0.5*ddt*kk0_V;
    fun1_fVI=funfVI+0.5*ddt*kk0_VI;

    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr101(fun1_fI,fun1_fII,fun1_fIII,fun1_fIV,fun1_fV,fun1_fVI,n,nn);

    [vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
        vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
        vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time+ddt/2);

    kk1_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
    kk1_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
    kk1_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
    kk1_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
    kk1_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
    kk1_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 

    %% CALCUL KK2
    fun2_fI=funfI+0.5*ddt*kk1_I;
    fun2_fII=funfII+0.5*ddt*kk1_II;
    fun2_fIII=funfIII+0.5*ddt*kk1_III;
    fun2_fIV=funfIV+0.5*ddt*kk1_IV;
    fun2_fV=funfV+0.5*ddt*kk1_V;
    fun2_fVI=funfVI+0.5*ddt*kk1_VI;

    [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
        gr101(fun2_fI,fun2_fII,fun2_fIII,fun2_fIV,fun2_fV,fun2_fVI,n,nn);

    [vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
        vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
        vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time+ddt/2);

    kk2_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
    kk2_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
    kk2_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
    kk2_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
    kk2_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
    kk2_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 

    %% CALCUL KK3
    fun3_fI=funfI+ddt*kk2_I;
    fun3_fII=funfII+ddt*kk2_II;
    fun3_fIII=funfIII+ddt*kk2_III;
    fun3_fIV=funfIV+ddt*kk2_IV;
    fun3_fV=funfV+ddt*kk2_V;
    fun3_fVI=funfVI+ddt*kk2_VI;

     [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
         gr101(fun3_fI,fun3_fII,fun3_fIII,fun3_fIV,fun3_fV,fun3_fVI,n,nn);

     [vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
        vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
        vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time+ddt);

    kk3_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
    kk3_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
    kk3_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
    kk3_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
    kk3_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
    kk3_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ; 

    %% ASSEMBLAGE RK4
    funfInew(1:nn,1:nn)=funfI(1:nn,1:nn)+ddt*...
        ((1/6)*kk0_I(1:nn,1:nn)+(1/3)*kk1_I(1:nn,1:nn)+(1/3)*kk2_I(1:nn,1:nn)+(1/6)*kk3_I(1:nn,1:nn));
    funfIInew(1:nn,1:nn)=funfII(1:nn,1:nn)+ddt*...
        ((1/6)*kk0_II(1:nn,1:nn)+(1/3)*kk1_II(1:nn,1:nn)+(1/3)*kk2_II(1:nn,1:nn)+(1/6)*kk3_II(1:nn,1:nn));
    funfIIInew(1:nn,1:nn)=funfIII(1:nn,1:nn)+ddt*...
        ((1/6)*kk0_III(1:nn,1:nn)+(1/3)*kk1_III(1:nn,1:nn)+(1/3)*kk2_III(1:nn,1:nn)+(1/6)*kk3_III(1:nn,1:nn));
    funfIVnew(1:nn,1:nn)=funfIV(1:nn,1:nn)+ddt*...
        ((1/6)*kk0_IV(1:nn,1:nn)+(1/3)*kk1_IV(1:nn,1:nn)+(1/3)*kk2_IV(1:nn,1:nn)+(1/6)*kk3_IV(1:nn,1:nn));
    funfVnew(1:nn,1:nn)=funfV(1:nn,1:nn)+ddt*...
        ((1/6)*kk0_V(1:nn,1:nn)+(1/3)*kk1_V(1:nn,1:nn)+(1/3)*kk2_V(1:nn,1:nn)+(1/6)*kk3_V(1:nn,1:nn));
    funfVInew(1:nn,1:nn)=funfVI(1:nn,1:nn)+ddt*...
        ((1/6)*kk0_VI(1:nn,1:nn)+(1/3)*kk1_VI(1:nn,1:nn)+(1/3)*kk2_VI(1:nn,1:nn)+(1/6)*kk3_VI(1:nn,1:nn));

    [funfInew,funfIInew,funfIIInew,funfIVnew,funfVnew,funfVInew]=...
        ds101(funfInew,funfIInew,funfIIInew,funfIVnew,funfVnew,funfVInew,n,nn);
    
    %% mise à jour du temps
    time=time+ddt;
    % mise a jour des données
    funfI=funfInew;funfII=funfIInew;funfIII=funfIIInew;
    funfIV=funfIVnew;funfV=funfVnew;funfVI=funfVInew;
    % HISTORIQUE
    xdays(ite)=(time-tinit)/(24*3600);

    % solution exacte
    funfIe=fun4_b(x_fI,y_fI,z_fI,time);
    funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
    funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
    funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
    funfVe=fun4_b(x_fV,y_fV,z_fV,time);
    funfVIe=fun4_b(x_fVI,y_fVI,z_fVI,time);

    % erreur mesurée
    err_fI=funfI-funfIe;
    err_fII=funfII-funfIIe;
    err_fIII=funfIII-funfIIIe;
    err_fIV=funfIV-funfIVe;
    err_fV=funfV-funfVe;
    err_fVI=funfVI-funfVIe;

    %% calcul d'erreur
    % en norme 1

    str='1';
    [~,~,~,~,~,~,nrmger]=...
        nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [~,~,~,~,~,~,nrmge]=...
        nrm101(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
    er1(ite)=nrmger/nrmge;

    % en norme 2

    str='2';
    [~,~,~,~,~,~,nrmger]=...
        nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [~,~,~,~,~,~,nrmge]=...
        nrm101(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
    er2(ite)=nrmger/nrmge;

    % en norme infinie

    str='infty';
    [nrmerI,nrmerII,nrmerIII,nrmerIV,nrmerV,nrmerVI,nrmger]=...
        nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
    [nrmeI,nrmeII,nrmeIII,nrmeIV,nrmeV,nrmeVI,nrmge]=...
        nrm101(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
    erinfty(ite)=nrmger/nrmge;

    %erreur 'brut'

    ermax(ite)=max(max([funfI-funfIe,funfII-funfIIe,funfIII-funfIIIe,funfIV-funfIVe,funfV-funfVe,funfVI-funfVIe]))/...
        max(max([abs(funfIe),abs(funfIIe),abs(funfIIIe),abs(funfIVe),abs(funfVe),abs(funfVIe)]));
    ermin(ite)=min(min([funfI-funfIe,funfII-funfIIe,funfIII-funfIIIe,funfIV-funfIVe,funfV-funfVe,funfVI-funfVIe]))/...
        max(max([abs(funfIe),abs(funfIIe),abs(funfIIIe),abs(funfIVe),abs(funfVe),abs(funfVIe)]));

    %% conservation
    str='int';
    [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,mass]=...
        nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,str);
    cons_mass(ite)=mass./mass_ref;
    
    ite=ite+1;
end
time2=cputime-time1;
disp(['temps de fonctionnement : ', num2str(time2)])

%% solution au temps final
 funfIe=fun4_b(x_fI,y_fI,z_fI,time);
 funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
 funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
 funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
 funfVe=fun4_b(x_fV,y_fV,z_fV,time);
 funfVIe=fun4_b(x_fVI,y_fVI,z_fVI,time);
 
%% erreur au temps final (erreur par point)
 err_fI=funfI-funfIe;
 err_fII=funfII-funfIIe;
 err_fIII=funfIII-funfIIIe;
 err_fIV=funfIV-funfIVe;
 err_fV=funfV-funfVe;
 err_fVI=funfVI-funfVIe;

%% erreur globale par face
 errI_38=max(max(abs(err_fI)));
 errII_38=max(max(abs(err_fII)));
 errIII_38=max(max(abs(err_fIII)));
 errIV_38=max(max(abs(err_fIV)));
 errV_38=max(max(abs(err_fV)));
 errVI_38=max(max(abs(err_fVI)));
 
 str='infty';
  [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
    nrm101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,str);
 [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
    nrm101(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);

time_res=cputime-tstart;

%% graphiques
if sauvegarde == 1
    mkdir(['./results-' date ])
    % save(['./results-' date '/ref_' num2str(ref) '_erreurdata_test_' num2str(coef) '.mat'])
end


figure(1);
plot_cs11(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe);
view(vvv)
title('Solution exacte')
if save_graph==1
    print('-dpng', ['./results-' date '/ref_' num2str(ref) '_solexacte_test_' num2str(coef) '.png'])
    savefig(['./results-' date '/ref_' num2str(ref) '_solexacte_test_' num2str(coef)]);
end

figure(2);
plot_cs11(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI);
title('Solution approchée')
view(vvv)
if save_graph==1
    print('-dpng', ['./results-' date '/ref_' num2str(ref) '_solapprochee_test_' num2str(coef) '.png'])
    savefig(['./results-' date '/ref_' num2str(ref) '_solapprochee_test_' num2str(coef)]);
end

figure(3);
plot_cs11(n,nn,err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI);colorbar;
title('Erreur')
view(vvv)
if save_graph==1
    print('-dpng', ['./results-' date '/ref_' num2str(ref) '_erreur_test_' num2str(coef) '.png'])
    savefig(['./results-' date '/ref_' num2str(ref) '_erreur_test_' num2str(coef)]);
end
  
figure(4);
plot(xdays,er1,'k-');hold on;grid;
plot(xdays,er2,'k--');hold on;
plot(xdays,erinfty,'k.-');hold off
%title('relative error - RK4')
legend('norme 1','norme 2','norme \infty','Location','NorthWest')
if save_graph==1
    print('-dpng', ['./results-' date '/ref_' num2str(ref) '_normerreur_test_' num2str(coef) '.png'])
    savefig(['./results-' date '/ref_' num2str(ref) '_normerreur_test_' num2str(coef)]);
end

figure(5)
plot(xdays,cons_mass-1)
title('Conservation de la masse')
if save_graph==1
    print('-dpng', ['./results-' date '/ref_' num2str(ref) '_conservation_test_' num2str(coef) '.png'])
    savefig(['./results-' date '/ref_' num2str(ref) '_conservation_test_' num2str(coef)]);
end

if save_graph==1
    data = fopen('AAA_RESULTS_SAVE.txt','a');
    fprintf(data,'%s\n',['date : ', date]);
    fprintf(data,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n',['test : ', num2str(coef)]);
    fprintf(data,'%s\n','---------- numerical data ---------');
    fprintf(data,'%s\n',['number of points  : ', num2str(n)] );
    fprintf(data,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(data,'%s\n',['cfl               : ', num2str(cfl)] );
    fprintf(data,'%s\n',['ordre du filtre   : ', num2str(opt_ftr)] );
    fprintf(data,'%s\n','-------- mathematical data --------');
    fprintf(data,'%s\n',['ndaymax           : ', num2str(ndaymax)] );
    fprintf(data,'%s\n',['angle in degree   : ', num2str(alphad*180/pi)] );
    fprintf(data,'%s\n',['(lambdap,thetap)  : (', num2str(lambda_p*180/pi),',',num2str(teta_p*180/pi),')']);
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n',['max(er_1)         : ', num2str(max(er1))] );
    fprintf(data,'%s\n',['max(er_2)         : ', num2str(max(er2))] );
    fprintf(data,'%s\n',['max(er_infty)     : ', num2str(max(erinfty))] );
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','  ');
    fprintf(data,'%s\n','  ');
    fclose(data);
end

if coupe == 1
    [ x,f ] = coupe_eq(funfI,funfII,funfIII,funfIV);
    n=500;
    nn=n+2;
    mod101
    funfIe=fun4_b(x_fI,y_fI,z_fI,time);
    funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
    funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
    funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
    [ xe,fe ] = coupe_eq(funfIe,funfIIe,funfIIIe,funfIVe);
    figure(10)
    plot(x,f,'o',xe,fe,'-')
    grid on;
    legend('Solution approchée','Solution exacte')
    xlabel('path N')
    if save_graph==1
        print('-dpng', ['./results-' date '/ref_' num2str(ref) '_coupefaceI_equateur_test_' num2str(coef) '.png'])
        savefig(['./results-' date '/ref_' num2str(ref) '_coupefaceI_test_' num2str(coef)]);
    end
end

fig_placier;
[er1(end) er2(end) erinfty(end)]