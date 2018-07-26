clc; clear all; close all;
format shorte
global n nn
global opt_ftr scheme

ftr1='redonnet2';
ftr2='redonnet4';
ftr3='redonnet6';
ftr4='redonnet8';
ftr5='redonnet10';
N=16;

scheme='compact4';
n=N-1; nn=N+1;

%% *** BLOC A ***
% options
opt_ftr =ftr1;
mod101
% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
for i=1:Nb
    clc;
    disp(['chargement matrice A : ' num2str(i/Nb*100) ' %'])
    vect=Id(1:Nb,i);
    [ sm ] = mat_ftr( vect );
    Da(1:Nb,i)=sm;
end
Da=full(Da);

%% *** BLOC B *************************************************************
% options
opt_ftr =ftr2;
mod101
% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
for i=1:Nb
    clc;
    disp(['chargement matrice B : ' num2str(i/Nb*100) ' %'])
    vect=Id(1:Nb,i);
    [ sm ] = mat_ftr( vect );
    Db(1:Nb,i)=sm;
end
Db=full(Db);

%% *** BLOC C *************************************************************
% options
opt_ftr =ftr3;
mod101
% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
for i=1:Nb
    clc;
    disp(['chargement matrice C : ' num2str(i/Nb*100) ' %'])
    vect=Id(1:Nb,i);
    [ sm ] = mat_ftr( vect );
    Dc(1:Nb,i)=sm;
end
Dc=full(Dc);

%% *** BLOC D *************************************************************
% options
opt_ftr =ftr4;
mod101
% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
for i=1:Nb
    clc;
    disp(['chargement matrice D : ' num2str(i/Nb*100) ' %'])
    vect=Id(1:Nb,i);
    [ sm ] = mat_ftr( vect );
    Dd(1:Nb,i)=sm;
end
Dd=full(Dd);

%% *** BLOC E *************************************************************
% options
opt_ftr =ftr5;
mod101
% construction de la matrice
Nb=6*N^2+2;
Id=speye(Nb,Nb);
for i=1:Nb
    clc;
    disp(['chargement matrice E : ' num2str(i/Nb*100) ' %'])
    vect=Id(1:Nb,i);
    [ sm ] = mat_ftr( vect );
    De(1:Nb,i)=sm;
end
De=full(De);


%% *** recherche des valeurs propres **************************************
Ea=eig(Da);
Eb=eig(Db);
Ec=eig(Dc);
Ed=eig(Dd);
Ee=eig(De);

%% *** figures ************************************************************
figure(1)
plot(Ea,'b+'); hold on
plot(Eb,'r+');
plot(Ec,'m+');
plot(Ed,'g+');
plot(Ee,'k+'); hold off
legend(ftr1,ftr2,ftr3,ftr4,ftr5)

figure(2)
plot(sort(abs(Ea),'descend'),'b-'); hold on
plot(sort(abs(Eb),'descend'),'r-');
plot(sort(abs(Ec),'descend'),'m-');
plot(sort(abs(Ed),'descend'),'g-');
plot(sort(abs(Ee),'descend'),'k-'); hold off
legend(ftr1,ftr2,ftr3,ftr4,ftr5,'Location','northeast')
title('|\lambda| in decreasing order')
grid on

figure(3)
plot(sort(real(Ea),'descend'),'b+'); hold on
plot(sort(real(Eb),'descend'),'r+');
plot(sort(real(Ec),'descend'),'m+');
plot(sort(real(Ed),'descend'),'g+');
plot(sort(real(Ee),'descend'),'k+'); hold off
legend(ftr1,ftr2,ftr3,ftr4,ftr5,'Location','northeast')
title('real part in decreasing order')
grid on

figure(4)
plot(sort(imag(Ea),'descend'),'b+'); hold on
plot(sort(imag(Eb),'descend'),'r+');
plot(sort(imag(Ec),'descend'),'m+');
plot(sort(imag(Ed),'descend'),'g+');
plot(sort(imag(Ee),'descend'),'k+'); hold off
legend(ftr1,ftr2,ftr3,ftr4,ftr5,'Location','northeast')
title('imaginary part in decreasing order')
grid on

fig_placier