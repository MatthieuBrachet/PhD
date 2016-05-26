% TEST PANEL II
clear all; clc; close all;
%% construction des variables globales
global n nn;
global x_fI y_fI z_fI x_fII y_fII z_fII x_fIII y_fIII z_fIII x_fIV y_fIV z_fIV x_fV y_fV z_fV x_fVI y_fVI z_fVI coef opt_ftr
% test de Williamson
global alphad tetac lambdac
% test de Nair et Machenhauer
global gamma rho0 teta_p lambda_p
% test de Nair et Jablonowski
global teta0 lambda0
% test de Nair et Lauritzen
global lambdac1 tetac1 lambdac2 tetac2

global pts_beta ptscr_beta;
%% *** OPTIONS ************************************************************
coef = 1;
opt_ftr =10;
%% *** Benchmarks data ****************************************************
 n=50;
 nn=n+2;
%% ************************************************************************
 if coef == 0
 %% test de Williamson
 alphad=3*pi/4;  
 lambdac=-pi/2;                                                           % longitude BUMP
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
 alphad=3*pi/4; 
 lambda0 = pi/2;
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
mod_1b

ndaymax=12;
time=24*3600*ndaymax;

%% problème
% initial condition
[funfI]=fun4_b(x_fI,y_fI,z_fI,time);
[funfII]=fun4_b(x_fII,y_fII,z_fII,time);
[funfIII]=fun4_b(x_fIII,y_fIII,z_fIII,time);
[funfIV]=fun4_b(x_fIV,y_fIV,z_fIV,time);
[funfV]=fun4_b(x_fV,y_fV,z_fV,time);
[funfVI]=fun4_b(x_fVI,y_fVI,z_fVI,time);

%% méthode 1
disp('----- case 1 :');
t1=cputime;

vb_fII=zeros(nn,4*(nn-1));
% Face II : TRANSFERT OF DATA OF FACE II
for iline1=1:nn, % boucle sur les lignes iso-xi du reseau II-beta
    vb_fII(iline1,1:nn-1)=funfII(iline1,1:nn-1);
end

% FACE V: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE V
for i=1:nn-1 % LOOP ON THE XI OF FACE V. i= order of encountering !!!
    betaspline(1:nn)=beta(nn-i+1,1:nn);
    funspline(1:nn)=funfV(nn-i+1,1:nn);
    ppspline=spline(betaspline,funspline);
    funbV1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end

% FACE III: TRANSFERT OF DATA OF FACE III
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funfIV(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
end

% FACE VI: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE Vi
for i=1:nn-1 % LOOP ON THE XI OF FACE VI. i= order of encountering !!!
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfVI(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbVI1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end

vb_fII1=vb_fII;

time_meth1=cputime-t1

%% méthode 2
disp('----- case 2 :')
t1=cputime;
vb_fII=zeros(nn,4*(nn-1));
% Face II 
for iline1=1:nn, 
    vb_fII(iline1,1:nn-1)=funfII(iline1,1:nn-1);
end

% Face V
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfV(nn-j+1,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfV(nn-j+1,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta);
for i=1:nn-1
    if i > nn/2
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(1:nn,i)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(1:nn,i)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
vb_fII(1:nn,nn-1+[1:nn-1])=FUN(nn:-1:1,1:nn-1);

% FACE III: TRANSFERT OF DATA OF FACE III
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funfIV(iline1,nn:-1:2); % symetrie sur adresses en i + inversion en j !
end

% Face VI
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfVI(j,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfVI(j,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta);

for i=1:nn-1
    if i > nn/2
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(1:nn,i)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(1:nn,i)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(1:nn,i)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
vb_fII(1:nn,3*nn-3+[1:nn-1])=FUN(nn:-1:1,1:nn-1);



time_meth2=cputime-t1


%% courbes

figure(1)
surf(abs(vb_fII1-vb_fII))
title('error')

max(max(abs(vb_fII1-vb_fII)))

figure(2)
subplot(121)
surf(vb_fII1)
title('old method')

subplot(122)
surf(vb_fII)
title('new method')





