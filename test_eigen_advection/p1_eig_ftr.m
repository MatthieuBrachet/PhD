clc; clear all; close all;
format longe
global n nn
global opt_ftr scheme

%% options
opt_ftr ='redonnet10';
scheme='compact4';
image='yes';

%% coef numériques
N=16;
n=N-1;
nn=N+1;
mod101

%% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
for i=1:Nb
    clc;
    disp(['chargement matrice : ' num2str(i/Nb*100) ' %'])
    vect=Id(1:Nb,i);
    [ sm ] = mat_ftr( vect );
    D(1:Nb,i)=sm;
end
D=full(D);

%% recherche des valeurs propres
E = eig(D);
lambda=max_real( E );
disp('valeurs propres : ok')

r=max(abs(E));
mr=max(abs(real(E)));
mi=max(abs(imag(E)));
disp(['paramètre de la CS : ' num2str(N)])
disp(['rayon spectrale    : ' num2str(r)])
disp(['max. real part     : ' num2str(mr)])
disp(['max. imag. part    : ' num2str(mi)])
theta=0:.001:2*pi;
x=r*cos(theta); y=r*sin(theta);

%% courbes
if strcmp(image,'yes')==1
    hFig=figure(1);
    plot(E,'b+'); hold on
    xlabel('real(\lambda)')
    ylabel('imag(\lambda)')
    plot(x,y,'k--','Linewidth',2); hold off
    set(hFig, 'Position', [50 50 500 500])
    grid on

    figure(2)
    plot(E,'b+')
    xlabel('real(\lambda)')
    ylabel('imag(\lambda)')
    grid on

    figure(3)
    surf(log10(abs(D)))
    shading interp
    colorbar
    view(2)

    figure(4)
    spy(D)
    
    figure(5)
    plot(sort(abs(E)),'k+')
    ylabel('|\lambda|')
    grid on

    fig_placier
end
lambda