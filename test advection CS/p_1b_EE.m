% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
% authors : Matthieu Brachet
%           Jean-Pierre Croisille
%
% ----------------------------------
%
% Euler Explicite + Filtrage à l'ordre 10
%
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
%% *** OPTIONS
%**************************************************************************
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
%
%% ************************************************************************
%
% test de transport avec simplification de J.-P. Croisille
%
%% Benchmarks data
 n=20;
 nn=n+2;
 cfl=0.9;
 itestop=100000;

 if coef == 0
 % test de Williamson
 alphad=pi/2;  
 lambdac=3*pi/2;                                                           % longitude BUMP
 tetac=pi/4;                                                               % latitude BUMP
 lambda_p=pi;                                                              % i.e. position du vortex nord
 teta_p=pi/2 - alphad;

 elseif coef == 1
 % test de Nair et Machenhauer
 lambda_p=pi/4;                                                            % i.e. position du vortex nord
 teta_p=-pi/4;
 rho0=3;
 gamma=5;
 
 elseif coef == 2
 % test de Nair et Jablonowski
 alphad=0; 
 
 lambda0 = 3*pi/2;
 teta0 = 0;
 lambda_p=pi;                                                              % position du pole nord, i.e. position du vortex nord
 teta_p=pi/2 - alphad;
 rho0=3;
 gamma=5;
 
 end
%% données du problème
mod_1b
ndaymax=12;
tmax=24*3600*ndaymax;
ddt=cfl*radius*dxi/u0;
itemax=floor(tmax/ddt);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUT BOUCLE EN TEMPS %
% ------------------------------------------------------------------------
% 5 - CALCUL DERIVEES HERMITIENNES SUR LES 6 RESEAUX DE GRANDS CERCLES
% *** PHASE 1 CHARGEMENT DES DONNEES, CELA COMPREND 
%      A- LE TRANSFERT SIMPLES DES VALEURS SUR DEUX FACES 
%      B- L'INTERPOLATION SPLINE CUBIQUE SUR LES DEUX AUTRES FACES (DE TYPE CROSS).
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

%% données pour RK4
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

if film==1
    %% film data
    nFrames = 20; % frame per second
    aviobj = avifile(['CSapprox_test' num2str(coef) '.avi'], 'fps',10);
end


%% Boucles RK 4 avec filtrage

for ite=1:itemax
clc; [ite itemax]

%% -----------------------------------------------------------------------
%%                 CALCUL DES ITERATIONS
%%  -----------------------------------------------------------------------

%% filtrage avant de commencer le calcul...
% N.B. le filtrage est effectué sur les grands cercles complets mais une
% seule fois pas iteration.

[funftI,funftII,funftIII,funftIV,funftV,funftVI]=...
    ftr_1b(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

funfI=funftI;funfII=funftII;funfIII=funftIII;
funfIV=funftIV;funfV=funftV;funfVI=funftVI;


%%   CALCUL RK4

if ite==itestop,
    figure(3); plot(aaa,'-r');
    figure(4); plot(bbb,'-b');
    break
end

%% iterations *************************************************************

%% CALCUL KK0
[grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr_1b(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn);

[vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
    vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
    vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI] = vitesse_1b(time);

kk0_I   = -(vitx_I.*grad_I(1:nn,1:nn,1) +  vity_I.*grad_I(1:nn,1:nn,2) + vitz_I.*grad_I(1:nn,1:nn,3)) ;
kk0_II  = -(vitx_II.*grad_II(1:nn,1:nn,1)+  vity_II.*grad_II(1:nn,1:nn,2) + vitz_II.*grad_II(1:nn,1:nn,3)) ;     
kk0_III = -(vitx_III.*grad_III(1:nn,1:nn,1)+  vity_III.*grad_III(1:nn,1:nn,2) + vitz_III.*grad_III(1:nn,1:nn,3)) ; 
kk0_IV  = -(vitx_IV.*grad_IV(1:nn,1:nn,1)+  vity_IV.*grad_IV(1:nn,1:nn,2) + vitz_IV.*grad_IV(1:nn,1:nn,3)) ; 
kk0_V   = -(vitx_V.*grad_V(1:nn,1:nn,1)+  vity_V.*grad_V(1:nn,1:nn,2) + vitz_V.*grad_V(1:nn,1:nn,3)) ; 
kk0_VI  = -(vitx_VI.*grad_VI(1:nn,1:nn,1)+  vity_VI.*grad_VI(1:nn,1:nn,2) + vitz_VI.*grad_VI(1:nn,1:nn,3)) ;

%% ASSEMBLAGE EE
funfInew(1:nn,1:nn)=funfI(1:nn,1:nn)+ddt*kk0_I(1:nn,1:nn);
funfIInew(1:nn,1:nn)=funfII(1:nn,1:nn)+ddt*kk0_II(1:nn,1:nn);
funfIIInew(1:nn,1:nn)=funfIII(1:nn,1:nn)+ddt*kk0_III(1:nn,1:nn);
funfIVnew(1:nn,1:nn)=funfIV(1:nn,1:nn)+ddt*kk0_IV(1:nn,1:nn);
funfVnew(1:nn,1:nn)=funfV(1:nn,1:nn)+ddt*kk0_V(1:nn,1:nn);
funfVInew(1:nn,1:nn)=funfVI(1:nn,1:nn)+ddt*kk0_VI(1:nn,1:nn);

[funfInew,funfIInew,funfIIInew,funfIVnew,funfVnew,funfVInew]=...
    ds_1b(funfInew,funfIInew,funfIIInew,funfIVnew,funfVnew,funfVInew,n,nn);

%% TIME UPDATE

time=time+ddt;

%% mise a jour des données
funfI=funfInew;funfII=funfIInew;funfIII=funfIIInew;
funfIV=funfIVnew;funfV=funfVnew;funfVI=funfVInew;

%% HISTORIC
% gestion en temps réel
xdays(ite)=(time-tinit)/(24*3600);

% solution exacte
funfIe=fun4_b(x_fI,y_fI,z_fI,time);funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time); funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
funfVe=fun4_b(x_fV,y_fV,z_fV,time); funfVIe=fun4_b(x_fVI,y_fVI,z_fVI,time);

% erreur mesurée
err_fI=funfI-funfIe;err_fII=funfII-funfIIe;
err_fIII=funfIII-funfIIIe;err_fIV=funfIV-funfIVe;
err_fV=funfV-funfVe;err_fVI=funfVI-funfVIe;

%% ERROR
% en norme 1

str='1';
[nrmerI,nrmerII,nrmerIII,nrmerIV,nrmerV,nrmerVI,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[nrmeI,nrmeII,nrmeIII,nrmeIV,nrmeV,nrmeVI,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
er1(ite)=nrmger/nrmge;

% en norme 2
str='2';
[nrmerI,nrmerII,nrmerIII,nrmerIV,nrmerV,nrmerVI,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[nrmeI,nrmeII,nrmeIII,nrmeIV,nrmeV,nrmeVI,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
er2(ite)=nrmger/nrmge;

% en norme infinie
str='infty';
[nrmerI,nrmerII,nrmerIII,nrmerIV,nrmerV,nrmerVI,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[nrmeI,nrmeII,nrmeIII,nrmeIV,nrmeV,nrmeVI,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
erinfty(ite)=nrmger/nrmge;

%erreur 'brut'
ermax(ite)=max(max([funfI-funfIe,funfII-funfIIe,funfIII-funfIIIe,funfIV-funfIVe,funfV-funfVe,funfVI-funfVIe]))/...
    max(max([abs(funfIe),abs(funfIIe),abs(funfIIIe),abs(funfIVe),abs(funfVe),abs(funfVIe)]));
ermin(ite)=min(min([funfI-funfIe,funfII-funfIIe,funfIII-funfIIIe,funfIV-funfIVe,funfV-funfVe,funfVI-funfVIe]))/...
    max(max([abs(funfIe),abs(funfIIe),abs(funfIIIe),abs(funfIVe),abs(funfVe),abs(funfVIe)]));

%% OPTIONS
if film==1
    figure(100); % Makes sure you use your desired frame.
    plot_cs5(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI);
    colorbar;
 
    aviobj = addframe(aviobj, getframe(gca));
end
if qquiv == 1
    figure(40)
    plot_quiver(nn,vitx_I,vitx_II, vitx_III, vitx_IV, vitx_V, vitx_VI, ...
        vity_I,vity_II, vity_III, vity_IV, vity_V, vity_VI, ...
        vitz_I,vitz_II, vitz_III, vitz_IV, vitz_V, vitz_VI, time);
end

end
if film == 1
    close(gcf)
    aviobj = close(aviobj);
end

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
 disp('erreur algebrique infty : ')
 err_38=max([errI_38,errII_38,errIII_38,errIV_38,errV_38,errVI_38])
 
 str='infty';
  [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
    nrm_1b(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn,str);
 [nrmI,nrmII,nrmIII,nrmIV,nrmV,nrmVI,nrmg]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);


%% graphiques


figure(35);
plot_cs5(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe);colorbar;
title('solution exacte')

 
figure(37);
plot_cs5(n,nn,funfI,funfII,funfIII,funfIV,funfV,funfVI);colorbar;
title('solution approchee - EE')

figure(39);
plot_cs5(n,nn,err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI);colorbar;
title('erreur algébrique - EE')

figure(1);
plot(xdays,er1,'k-');hold on;grid;
plot(xdays,er2,'k--');hold on;
plot(xdays,erinfty,'k.');
legend('norme 1','norme 2','norme infinie')
title('erreur globale - EE')
  
figure(2);
plot(xdays,ermax,'k-');hold on;grid;
plot(xdays,ermin,'k--');
legend('max','min')
title('erreur relative - EE')

format shortE
disp('erreur relative L2 : ')
max(er2)

data_save(n,nn,funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe)
