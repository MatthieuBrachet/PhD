function [xa,fa] = iterations_coupe(Na,ndaymax)



global n nn;
global radius u0 dxi;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII;
global x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI;
global ite

%% ************************************************************************
%
%                              BLOC A :
%
%  ************************************************************************
n=Na;
nn=n+2;
mod101
cfl=.9;
tmax=24*3600*ndaymax;
ddt=cfl*radius*dxi/u0;
itemax=floor(tmax/ddt);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUT BOUCLE EN TEMPS %
% ------------------------------------------------------------------------
% 5 - CALCUL DERIVEES HERMITIENNES SUR LES 6 RESEAUX DE GRANDS CERCLES
% *** PHASE 1 CHARGEMENT DES DONNEES, CELA COMPREND 
%      A- LE TRANSFERT SIMPLES DES VALEURS SUR DEUX FACES 
%      B- L'INTERPOLATION SPLINE CUBIQUE SUR LES DEUX AUTRES FACES (DE TYPE CROSS).
% *** PHASE 2: INTERPOLATION DE LA DERIVEE HERMITIENNE.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - -----------------------------------------------------------------------
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
[vwk0]=fun4_b(0,0,0,0);
lambda=zeros(nn,nn);
teta=zeros(nn,nn);
radius1=zeros(nn,nn);
elambda_x=zeros(nn,nn);
elambda_y=zeros(nn,nn);
elambda_z=zeros(nn,nn);
eteta_x=zeros(nn,nn);
eteta_y=zeros(nn,nn);
eteta_z=zeros(nn,nn);

% données pour RK4
tinit=0.;
% initial condition
[funfI]=fun4_b(x_fI,y_fI,z_fI,tinit);
[funfII]=fun4_b(x_fII,y_fII,z_fII,tinit);
[funfIII]=fun4_b(x_fIII,y_fIII,z_fIII,tinit);
[funfIV]=fun4_b(x_fIV,y_fIV,z_fIV,tinit);
[funfV]=fun4_b(x_fV,y_fV,z_fV,tinit);
[funfVI]=fun4_b(x_fVI,y_fVI,z_fVI,tinit);
funfInew=zeros(nn,nn);funfIInew=zeros(nn,nn);funfIIInew=zeros(nn,nn);
funfIVnew=zeros(nn,nn);funfVnew=zeros(nn,nn);funfVInew=zeros(nn,nn);
kk0_I=zeros(nn,nn);kk1_I=zeros(nn,nn);kk2_I=zeros(nn,nn);kk3_I=zeros(nn,nn);
kk0_II=zeros(nn,nn);kk1_II=zeros(nn,nn);kk2_II=zeros(nn,nn);kk3_II=zeros(nn,nn);
kk0_III=zeros(nn,nn);kk1_III=zeros(nn,nn);kk2_III=zeros(nn,nn);kk3_III=zeros(nn,nn);
kk0_IV=zeros(nn,nn);kk1_IV=zeros(nn,nn);kk2_IV=zeros(nn,nn);kk3_IV=zeros(nn,nn);
kk0_V=zeros(nn,nn);kk1_V=zeros(nn,nn);kk2_V=zeros(nn,nn);kk3_V=zeros(nn,nn);
kk0_VI=zeros(nn,nn);kk1_VI=zeros(nn,nn);kk2_VI=zeros(nn,nn);kk3_VI=zeros(nn,nn);

time=tinit;

% Boucles RK 4 avec filtrage
xdays(1)=0;
for ite=1:itemax
clc; [ite itemax]
% ------------------------------------------------------------------------
%                 CALCUL DES ITERATIONS
%  -----------------------------------------------------------------------
% filtrage avant de commencer le calcul...
% N.B. le filtrage est effectué sur les grands cercles complets mais une
% seule fois pas iteration.

[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_mixte101(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

funfI=funftI;funfII=funftII;funfIII=funftIII;
funfIV=funftIV;funfV=funftV;funfVI=funftVI;

% iterations
% CALCUL KK0
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

% CALCUL KK1
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

% CALCUL KK2
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

% CALCUL KK3
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

% ASSEMBLAGE RK4
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
% mise à jour du temps
time=time+ddt;
% mise a jour des données
funfI=funfInew;funfII=funfIInew;funfIII=funfIIInew;
funfIV=funfIVnew;funfV=funfVnew;funfVI=funfVInew;
% HISTORIQUE
xdays(ite)=(time-tinit)/(24*3600);
end
ref=floor(10000*now);
[ xa,fa ] = coupe_eq(funfI,funfII,funfIII,funfIV);


end

