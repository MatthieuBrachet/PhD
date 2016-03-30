clear
clc
close all

N = input('le nombre est  :')
nbcercle = 0;
tic;
X = zeros(N,1);
Y = zeros(N,1);
for i=1:N;
    x(i)=rand;
    y(i)=rand;
    if x(i)^2+y(i)^2<=1
        nbcercle= nbcercle + 1;
    end
end
  
  tic;
  X=rand(N,1);
  Y=rand(N,1);
  Nbcercle=length(find(X.^2+Y.^2<=1));
  tempsrequis=toc;
  disp (['le temps nécessaire pour la méthode 2 vaut : ' num2str(tempsrequis) 's'])


piMMC= 4*nbcercle/N;
disp(['pi a été estimé à ' num2str(piMMC)])

Xcercle=0:0.01:1;
Ycercle=sqrt(1-Xcercle.^2);
hold on
plot(Xcercle,Ycercle,'r','lineWidth',2)
plot(X,Y,'.k','markerSize',2)
title('Méthode de Monte Carlo : ');