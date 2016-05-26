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
 n=30;
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

va_fII=zeros(4*(nn-1),nn); % VALEURS =DONNEES RESEAU FACE II SELON ALPHA
% Face II : TRANSFERT OF DATA OF FACE II
for jline1=1:nn, % boucle sur les lignes iso-eta du reseaux IIa
    va_fII(1:nn-1,jline1)=funfII(1:nn-1,jline1);
end

% FACE III: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE II
for i=1:nn-1 % LOOP ON THE XI OF FACE II
    betaspline(1:nn)=beta(i,1:nn);
    funspline(1:nn)=funfIII(i,1:nn);
    ppspline=spline(betaspline,funspline);
    funbIII1(1:nn)=ppval(ppspline,betacr(i,1:nn));
    va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
end

% FACE IV: TRANSFERT OF DATA
for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funfIV(1:nn-1,nn-jline1+1); % symetrie sur adresses en j !    
end

% FACE I: SPLINE INTERPOLATION A L'AIDE DES ANGLES BETA
% INVERSION OF THE LOOP ON THE ISO-ETA LINES AND ON THE XI OF FACE I
for i=1:nn-1, % LOOP ON THE XI OF FACE I
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funfI(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbI1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
end

va_fII1=va_fII;

time_meth1=cputime-t1

%% méthode 2
disp('----- case 2 :')
t1=cputime;

% Face II : TRANSFERT OF DATA OF FACE I 
va_fI=zeros(4*(nn-1),nn); 
for jline1=1:nn, % boucle sur les lignes iso-eta du reseaux Ia
    va_fII(1:nn-1,jline1)=funfII(1:nn-1,jline1);
end

% Face III
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfIII(j,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfIII(j,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta(1:nn*nn));

for i=1:nn-1
    if i <= nn/2
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(i,1:nn)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(i,1:nn)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
va_fII(nn-1+[1:nn-1],1:nn)=FUN(1:nn-1,1:nn);




% Face IV 
for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funfIV(1:nn-1,nn-jline1+1);
end

% Face I
compteur=1;
for j=1:nn
    for i=1:nn
        if rem(j,2) == 1
            fun(compteur)=funfI(j,i);
        elseif rem(j,2) == 0
            fun(compteur)=funfI(j,nn-i+1);
        end
            compteur=compteur+1;
    end
end
ppspline=spline(pts_beta,fun);
funspl=ppval(ppspline,ptscr_beta);

for i=1:nn-1
    if i > nn/2
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[1:nn]);
        else
            FUN(i,1:nn)=funspl(i*nn-[1:nn]+1);
        end
    else
        if rem(i,2)== 1
            FUN(i,1:nn)=funspl((i-1)*nn+[nn:-1:1]);
        else
            FUN(i,1:nn)=funspl(i*nn-[nn:-1:1]+1);
        end
    end
end
va_fII(3*nn-3+[1:nn-1],1:nn)=FUN(1:nn-1,1:nn);
time_meth2=cputime-t1


%% courbes

figure(1)
surf(abs(va_fII1-va_fII))
title('error')

max(max(abs(va_fII1-va_fII)))

figure(2)
subplot(121)
surf(va_fII1)
title('old method')

subplot(122)
surf(va_fII)
title('new method')
