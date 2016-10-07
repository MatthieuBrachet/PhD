clc; clear all; close all;

% MODULE PROBLEM FOR THE CUBED SPHERE
% ----------------------------------
% authors : Matthieu Brachet
%           Jean-Pierre Croisille
% ----------------------------------
% RK4 + Filtrage
clear all; clc; close all;
%% construction des variables globales
global n nn radius u0 dxi;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI;
global ite itestop coef opt_ftr ftr
% test de Williamson
global alphad tetac lambdac
% test de Nair et Machenhauer
global gamma rho0 teta_p lambda_p
% test de Nair et Jablonowski
global teta0 lambda0
% test de Nair et Lauritzen
global lambdac1 tetac1 lambdac2 tetac2
%% *** OPTIONS ************************************************************
% si coef = 0, test 1 de Williamson (solid body rotation on the sphere)
%    coef = 1, test de Nair et Machenhauer  (deformational flow test - 
%                                                    stationnary vortex)
%    coef = 2, test de Nair, Jablonowski (moving vortices on the sphere)
%    coef = 3, test de Nair, Lauritzen (slotted cylinder) ( = Zaleska)
coef = 0;
% si save_graph = 1 : enregistrer les graphiques et les données dans TEST_SAVE.txt
%    save_graph = 0 : ne pas enregistrer
save_graph = 1;
% coupe = 0 : pas de coupe le long de l'équateur de la face 1
%         1 : coupe.
coupe = 1;
%% *** Benchmarks data ****************************************************
 n=80;
 nn=n+2;
 cfl=0.9;
 ndaymax=12;
%% *** filtres choisis ****************************************************
ftra=10;
ftrb=4;
ftrc=2;
%% ************************************************************************
opt_ftr=10;
 if coef == 0
 % test de Williamson
 alphad=0;  
 lambdac=pi/2;                                                           % longitude BUMP
 tetac=0;                                                                  % latitude BUMP
 lambda_p=pi;                                                              % position du pole nord, i.e. position du vortex nord
 teta_p=pi/2 - alphad;
 elseif coef == 1
 % test de Nair et Machenhauer
 lambda_p=0;                                                            % position du pole nord, i.e. position du vortex nord
 teta_p=3*pi/2;
 rho0=3;
 gamma=5;
 elseif coef == 2
 % test de Nair et Jablonowski
 alphad=0; 
 lambda0 = 3*pi/2;
 teta0 = 0;
 lambda_p=pi;                                                              % position du pole nord à t=0, i.e. position du vortex nord à t=0
 teta_p=pi/2 - alphad;
 rho0=3;
 gamma=5;
 elseif coef == 3
 % test de Nair, Lauritzen
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
mod_1b
tmax=24*3600*ndaymax;
ddt=cfl*radius*dxi/u0;
itemax=floor(tmax/ddt);
%% -----------------------------------------------------------------------
% DEBUT BOUCLE EN TEMPS %
% ------------------------------------------------------------------------


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
% initial condition filtre = 10
[funfI10]=fun4_b(x_fI,y_fI,z_fI,tinit); [funfII10]=fun4_b(x_fII,y_fII,z_fII,tinit);
[funfIII10]=fun4_b(x_fIII,y_fIII,z_fIII,tinit); [funfIV10]=fun4_b(x_fIV,y_fIV,z_fIV,tinit);
[funfV10]=fun4_b(x_fV,y_fV,z_fV,tinit); [funfVI10]=fun4_b(x_fVI,y_fVI,z_fVI,tinit);

% initial condition filtre = 8
[funfI8]=fun4_b(x_fI,y_fI,z_fI,tinit); [funfII8]=fun4_b(x_fII,y_fII,z_fII,tinit);
[funfIII8]=fun4_b(x_fIII,y_fIII,z_fIII,tinit); [funfIV8]=fun4_b(x_fIV,y_fIV,z_fIV,tinit);
[funfV8]=fun4_b(x_fV,y_fV,z_fV,tinit); [funfVI8]=fun4_b(x_fVI,y_fVI,z_fVI,tinit);

% initial condition filtre = 6
[funfI6]=fun4_b(x_fI,y_fI,z_fI,tinit); [funfII6]=fun4_b(x_fII,y_fII,z_fII,tinit);
[funfIII6]=fun4_b(x_fIII,y_fIII,z_fIII,tinit); [funfIV6]=fun4_b(x_fIV,y_fIV,z_fIV,tinit);
[funfV6]=fun4_b(x_fV,y_fV,z_fV,tinit); [funfVI6]=fun4_b(x_fVI,y_fVI,z_fVI,tinit);

time=tinit;

%% Boucles RK 4 avec filtrage
xdays(1)=0;
for ite=1:itemax
clc; [ite itemax]

%% ------------------------------------------------------------------------
%%                 CALCUL DES ITERATIONS
%%  -----------------------------------------------------------------------

%% HISTORIQUE
% gestion en temps réel
xdays(ite)=(time-tinit)/(24*3600);

% solution exacte
funfIe=fun4_b(x_fI,y_fI,z_fI,time);
funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
funfVe=fun4_b(x_fV,y_fV,z_fV,time);
funfVIe=fun4_b(x_fVI,y_fVI,z_fVI,time);



%% filtrage avant de commencer le calcul...

%% *** ITERATIONS *********************************************************

%% filtre ordre a
opt_ftr=ftra;
[ ftr ] = filtre( na , opt_ftr );
[funfInew, funfIInew, funfIIInew, funfIVnew, funfVnew, funfVInew] = iteration(funfI10, funfII10,funfIII10,funfIV10,funfV10,funfVI10,ddt,time);
funfI10=funfInew;funfII10=funfIInew;funfIII10=funfIIInew;
funfIV10=funfIVnew;funfV10=funfVnew;funfVI10=funfVInew;
% erreur mesurée
err_fI=funfI10-funfIe;
err_fII=funfII10-funfIIe;
err_fIII=funfIII10-funfIIIe;
err_fIV=funfIV10-funfIVe;
err_fV=funfV10-funfVe;
err_fVI=funfVI10-funfVIe;
%% calcul d'erreur
% en norme 1
str='1';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
er1_10(ite)=nrmger/nrmge;
% en norme 2
str='2';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
er2_10(ite)=nrmger/nrmge;
% en norme infinie
str='infty';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
erinfty_10(ite)=nrmger/nrmge;

% *************************************************************************

%% filtre ordre b
opt_ftr=ftrb;
[ ftr ] = filtre( na , opt_ftr );
[funfInew, funfIInew, funfIIInew, funfIVnew, funfVnew, funfVInew] = iteration(funfI8, funfII8,funfIII8,funfIV8,funfV8,funfVI8,ddt,time);
funfI8=funfInew;funfII8=funfIInew;funfIII8=funfIIInew;
funfIV8=funfIVnew;funfV8=funfVnew;funfVI8=funfVInew;
% erreur mesurée
err_fI=funfI8-funfIe;
err_fII=funfII8-funfIIe;
err_fIII=funfIII8-funfIIIe;
err_fIV=funfIV8-funfIVe;
err_fV=funfV8-funfVe;
err_fVI=funfVI8-funfVIe;
%% calcul d'erreur
% en norme 1
str='1';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
er1_8(ite)=nrmger/nrmge;
% en norme 2
str='2';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
er2_8(ite)=nrmger/nrmge;
% en norme infinie
str='infty';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
erinfty_8(ite)=nrmger/nrmge;


%% ************************************************************************

%% filtre ordre c
opt_ftr=ftrc;
[ ftr ] = filtre( na , opt_ftr );
[funfInew, funfIInew, funfIIInew, funfIVnew, funfVnew, funfVInew] = iteration(funfI6, funfII6,funfIII6,funfIV6,funfV6,funfVI6,ddt,time);
funfI6=funfInew;funfII6=funfIInew;funfIII6=funfIIInew;
funfIV6=funfIVnew;funfV6=funfVnew;funfVI6=funfVInew;
% erreur mesurée
err_fI=funfI6-funfIe;
err_fII=funfII6-funfIIe;
err_fIII=funfIII6-funfIIIe;
err_fIV=funfIV6-funfIVe;
err_fV=funfV6-funfVe;
err_fVI=funfVI6-funfVIe;
%% calcul d'erreur
% en norme 1
str='1';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
er1_6(ite)=nrmger/nrmge;
% en norme 2
str='2';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
er2_6(ite)=nrmger/nrmge;
% en norme infinie
str='infty';
[~,~,~,~,~,~,nrmger]=...
    nrm_1b(err_fI,err_fII,err_fIII,err_fIV,err_fV,err_fVI,n,nn,str);
[~,~,~,~,~,~,nrmge]=...
    nrm_1b(funfIe,funfIIe,funfIIIe,funfIVe,funfVe,funfVIe,n,nn,str);
erinfty_6(ite)=nrmger/nrmge;

%% ************************************************************************


%% mise à jour du temps
time=time+ddt;
end

ref=floor(10000*now);
%% graphiques  
figure(1);
subplot(131)
plot(xdays,er1_10,'r-');hold on;grid on;
plot(xdays,er1_8,'g-');hold on;
plot(xdays,er1_6,'b.'); hold on
legend(['filtre', num2str(ftra)],['filtre ', num2str(ftrb)],['filtre ', num2str(ftrc)]);
title('erreur en norme 1')

subplot(132)
plot(xdays,er2_10,'r-');hold on;grid on;
plot(xdays,er2_8,'g-');hold on;
plot(xdays,er2_6,'b.'); hold on
legend(['filtre', num2str(ftra)],['filtre ', num2str(ftrb)],['filtre ', num2str(ftrc)]);
title('erreur en norme 2')

subplot(133)
plot(xdays,erinfty_10,'r-');hold on;grid on;
plot(xdays,erinfty_8,'g-');hold on;
plot(xdays,erinfty_6,'b.'); hold on
legend(['filtre', num2str(ftra)],['filtre ', num2str(ftrb)],['filtre ', num2str(ftrc)])
title('erreur en norme infinie')

if save_graph==1
    print('-dpng', ['./results_ftr/' date 'ref_' num2str(ref) '_normerreur_test_' num2str(coef) '.png'])
end

if coupe == 1
    [ ~,f_10 ] = coupe_eq(funfI10,funfII10,funfIII10,funfIV10);
    [ ~,f_8 ] = coupe_eq(funfI8,funfII8,funfIII8,funfIV8);
    [ x,f_6 ] = coupe_eq(funfI6,funfII6,funfIII6,funfIV6);
    n=100;
    nn=n+2;
    mod_1b
    funfIe=fun4_b(x_fI,y_fI,z_fI,time);
    funfIIe=fun4_b(x_fII,y_fII,z_fII,time);
    funfIIIe=fun4_b(x_fIII,y_fIII,z_fIII,time);
    funfIVe=fun4_b(x_fIV,y_fIV,z_fIV,time);
    [ xe,fe ] = coupe_eq(funfIe,funfIIe,funfIIIe,funfIVe);
    
    figure(10)
    plot(x,f_10,'b-o'); hold on;
    plot(x,f_8,'g-d');hold on
    plot(x,f_6,'r-x'); hold on;
    plot(xe,fe,'k-')
    grid on;
    legend(['filtre ', num2str(ftra)],['filtre ', num2str(ftrb)],['filtre ', num2str(ftrc)],'solution exacte')
    %xlabel('equateur - face II')
    title('coupe de la solution le long de l''equateur')
    if save_graph==1
        print('-dpng', ['./results_ftr/' date 'ref_' num2str(ref) '_coupefaceI_equateur_test_' num2str(coef) '.png'])
    end
end

if save_graph==1
    %ouvre un fichier ou le créé
    fid = fopen('./results_ftr/TEST_SAVE.txt','a');
    %écrit dans ce fichier, fid est sa reference pour matlab
    fprintf(fid,'%s\n',['date : ', date]);
    fprintf(fid,'%s\n',['date : ', num2str(ref)]);
    fprintf(fid,'%s\n','******************************');
    fprintf(fid,'%s\n',['test : ', num2str(coef)]);
    fprintf(fid,'%s\n','------- numerical data -------');
    fprintf(fid,'%s\n',['number of points : ', num2str(n)] );
    fprintf(fid,'%s\n',['time step        : ', num2str(ddt)] );
    fprintf(fid,'%s\n',['cfl              : ', num2str(cfl)] );
    fprintf(fid,'%s\n','----- mathematical data ------');
    fprintf(fid,'%s\n',['angle in degree  : ', num2str(alphad*180/pi)] );
    fprintf(fid,'%s\n','******************************');
    fprintf(fid,'%s\n',['*** filtre ordre ', num2str(ftra) ]);
    fprintf(fid,'%s\n',['max(er_1)        : ', num2str(max(er1_10))] );
    fprintf(fid,'%s\n',['max(er_2)        : ', num2str(max(er2_10))] );
    fprintf(fid,'%s\n',['max(er_infty)    : ', num2str(max(erinfty_10))] );
    fprintf(fid,'%s\n',['max              : ', num2str(max(f_10))] );
    fprintf(fid,'%s\n',['*** filtre ordre ', num2str(ftrb) ]);
    fprintf(fid,'%s\n',['max(er_1)        : ', num2str(max(er1_8))] );
    fprintf(fid,'%s\n',['max(er_2)        : ', num2str(max(er2_8))] );
    fprintf(fid,'%s\n',['max(er_infty)    : ', num2str(max(erinfty_8))] );
    fprintf(fid,'%s\n',['max              : ', num2str(max(f_8))] );
    fprintf(fid,'%s\n',['*** filtre ordre ', num2str(ftrc) ]);
    fprintf(fid,'%s\n',['max(er_1)        : ', num2str(max(er1_6))] );
    fprintf(fid,'%s\n',['max(er_2)        : ', num2str(max(er2_6))] );
    fprintf(fid,'%s\n',['max(er_infty)    : ', num2str(max(erinfty_6))] );
    fprintf(fid,'%s\n',['max              : ', num2str(max(f_6))] );
    fprintf(fid,'%s\n','******************************');
    fprintf(fid,'%s\n',['ndaymax          : ', num2str(ndaymax)] );
    fprintf(fid,'%s\n','******************************');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    fprintf(fid,'%s\n','  ');
    %n'oublie pas de fermer le fichier sinon tu ne peux pas le lire
    fclose(fid);
end

