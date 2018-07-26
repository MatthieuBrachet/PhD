clc; clear all; close all;
global n nn;
global radius u0 dxi
global coef opt_ftr scheme
global alphad tetac lambdac
global teta_p lambda_p



N=16;
n=N-1;
nn=N+1;
cfl=1;

alphaA=35/180*pi;
alphaB=40/180*pi;
alphaC=45/180*pi;

scheme='compact4';

%% *** bloc A *************************************************************
opt_ftr ='redonnet10';
coef=0;

alphad=alphaA; 
lambdac=3*pi/2;
tetac=0;
lambda_p=pi;
teta_p=pi/2 - alphad;

mod101
ddt=cfl*radius*dxi/u0;

% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
time=0;
for i=1:Nb
    clc;
    disp(['chargement matrice A : ' num2str(i/Nb*100)])
    vect=Id(1:Nb,i);
    [ sm ] = second_mem( vect,time );
    Da(1:Nb,i)=sm;
end

% recherche des valeurs propres
Ea = eig(Da);


%% *** bloc B *************************************************************
opt_ftr ='redonnet10';
coef=0;

alphad=alphaB; 
lambdac=3*pi/2;
tetac=0;
lambda_p=pi;
teta_p=pi/2 - alphad;

mod101
ddt=cfl*radius*dxi/u0;

% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
time=0;
for i=1:Nb
    clc;
    disp(['chargement matrice B : ' num2str(i/Nb*100)])
    vect=Id(1:Nb,i);
    [ sm ] = second_mem( vect,time );
    Db(1:Nb,i)=sm;
end

% recherche des valeurs propres
Eb = eig(Db);



%% *** bloc C *************************************************************
opt_ftr ='redonnet10';
coef=0;

alphad=alphaC; 
lambdac=3*pi/2;
tetac=0;
lambda_p=pi;
teta_p=pi/2 - alphad;

mod101
ddt=cfl*radius*dxi/u0;

% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
time=0;
for i=1:Nb
    clc;
    disp(['chargement matrice C : ' num2str(i/Nb*100)])
    vect=Id(1:Nb,i);
    [ sm ] = second_mem( vect,time );
    Dc(1:Nb,i)=sm;
end

% recherche des valeurs propres
Ec = eig(Dc);


%% *** PLOT **************************************************************
hFig=figure(1);
plot(Ea,'b+'); hold on
plot(Eb,'r+');
plot(Ec,'g+');hold off
legend(['\alpha = ' num2str(alphaA*180/pi)],['\alpha = ' num2str(alphaB*180/pi)],['\alpha = ' num2str(alphaC*180/pi)])
xlabel('real(\lambda)')
ylabel('imag(\lambda)')
set(hFig, 'Position', [50 50 500 500])

xxx=-3:0.01:3;
[X,Y]=meshgrid(xxx,xxx);
L=(X+1i*Y);
AMPLI=1+L+.5*L.^2+(L.^3)/6+(L.^4)/24;
AMPLI=abs(AMPLI);
hFig2=figure(2);
contour(X,Y,(AMPLI<1),'r-','Linewidth',2); hold on
axis([-3 .5 -3 3])
set(hFig2, 'Position', [50 50 350 600])
plot(ddt*Ea,'b+'); hold on
plot(ddt*Eb,'r+');
plot(ddt*Ec,'g+');hold off
xlabel('ddt \times real(\lambda)')
ylabel('ddt \times imag(\lambda)')

figure(3);
plot(sort(abs(Ea)),'b+'); hold on
plot(sort(abs(Eb)),'r+');
plot(sort(abs(Ec)),'g+'); hold off
legend(['\alpha = ' num2str(alphaA*180/pi)],['\alpha = ' num2str(alphaB*180/pi)],['\alpha = ' num2str(alphaC*180/pi)])
ylabel('|\lambda|')

figure(4);
plot(sort(real(Ea)),'b+'); hold on
plot(sort(real(Eb)),'r+');
plot(sort(real(Ec)),'g+'); hold off
legend(['\alpha = ' num2str(alphaA*180/pi)],['\alpha = ' num2str(alphaB*180/pi)],['\alpha = ' num2str(alphaC*180/pi)])
ylabel('real(\lambda)')