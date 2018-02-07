clc; clear all; close all;

lambda=1.6883;
teta=linspace(.0001,pi,1000);
time='rk4';
c=4;
space='implicit';
f=10;
filtre='redonnet';
delta=0;

[ g ] = atsf(lambda, teta,time,c,space,f,filtre,delta);

dissip=abs(g);
z=g./dissip;
ang=imag(log(z));
ang=ang.*(ang<0)+(ang-2*pi).*(ang>0);
disper=-ang./(lambda.*teta);



figure(1)
plot(teta,dissip)
axis([0 pi 0 1.1])

figure(2)
plot(g)

figure(3)
plot(teta, dissip,teta,disper,'Linewidth',2)
legend('dissipation','dispersion','Location','southwest')
xlabel('\theta')
xticks([0 pi/4 pi/2 3*pi/4 pi])
xticklabels({'0','\pi/4','\pi/2','3\pi /4','\pi'})
grid on
axis([0 pi 0 max(1,max(disper)+.01)])
title(['\lambda = ', num2str(lambda)])