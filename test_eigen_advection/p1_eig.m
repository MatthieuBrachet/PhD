clc; clear all; close all;
format longe
global n nn
global radius u0 dxi
global coef opt_ftr scheme
global alphad tetac lambdac
global teta_p lambda_p

%% options
opt_ftr ='redonnet10';
scheme='compact4';
coef=0;
image='yes';

%% coef numériques
N=16;
n=N-1;
nn=N+1;
mod101
cfl=.7;
ddt=cfl*radius*dxi/u0;

%% coef. probleme
alphad=pi/4;  
lambdac=3*pi/2;
tetac=0;
lambda_p=pi;
teta_p=pi/2 - alphad;

%% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
time=0;
for i=1:Nb
    clc;
    disp(['chargement matrice : ' num2str(i/Nb*100) ' %'])
    vect=Id(1:Nb,i);
    [ sm ] = second_mem( vect,time );
    D(1:Nb,i)=sm;
end

%% recherche des valeurs propres
E = eig(D);
lambda=max_imag( E );
disp('valeurs propres : ok')

r=max(abs(E));
mr=max(abs(real(E)));
mi=max(abs(imag(E)));
disp(['paramètre de la CS : ' num2str(N)])
disp(['rayon spectrale    : ' num2str(r)])
disp(['max. real part     : ' num2str(mr)])
disp(['max. imag. part    : ' num2str(mi)])
theta=0:.001:2*pi;
x=r.*cos(theta); y=r.*sin(theta);

%% courbes
if strcmp(image,'yes')==1
    hFig=figure(1);
    plot(E,'b+'); hold on
    xlabel('real(\lambda)')
    ylabel('imag(\lambda)')
    plot(x,y,'k--','Linewidth',2); hold off
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
    plot(ddt*E,'b+'); hold on
    xlabel('real(\lambda)')
    ylabel('imag(\lambda)')
    plot(ddt*x,ddt*y,'k--','Linewidth',2); hold off

    figure(3)
    plot(E,'b+')
    xlabel('real(\lambda)')
    ylabel('imag(\lambda)')
    grid on

    figure(4)
    spy(D)

    fig_placier
end
lambda