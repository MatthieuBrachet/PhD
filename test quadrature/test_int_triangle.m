clc; clear all; close all;

NN=[3:100:3000];

for p=1:length(NN)
    n=NN(p);
    x=linspace(0,1,n);
    [X,Y]=meshgrid(x,x);
    
    TU=triu(ones(n,n));
    pp=6;
    qq=7;
    Ft=(X.^pp.*Y.^qq).*TU(:,n:-1:1);

    MAT=ones(n,n)-.5*eye(n,n);
    MAT(1,:)=.5;
    MAT(n,:)=.5;
    MAT(:,1)=.5;
    MAT(:,n)=.5;
    MAT(1,1)=1/6;
    MAT(1,n)=1/4;
    MAT(n,n)=1/8;
    MAT=MAT(:,n:-1:1);
    
    INT(p)=(1/(n-1)^2)*sum(sum(MAT.*Ft));
    int=factorial(pp)*factorial(qq)./(factorial(pp+qq+2));

    e(p)=(INT(p)-int)/int;
end

figure(1)
loglog(1./(NN-1),abs(e))
grid on

figure(2)
plot(1./(NN-1),e)
title('erreur')
grid on

figure(3)
plot(1./(NN-1),INT-INT(1))
title('deviation - integrale')
grid on

% figure(4)
% surf(X,Y,Ft)
% shading interp
% colorbar

figure(5); spy(MAT-1);
figure(6); spy(Ft)
figure(7); spy(TU(:,n:-1:1))

fig_placier

p=polyfit(log(1./(NN-1)),log(abs(e)),1);
disp(['ordre estimé : ' num2str(p(1))])

