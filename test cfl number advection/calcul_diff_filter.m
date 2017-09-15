clc; clear all; close all
format long

%% schema testé
order=10;
genre='visbal';

%% courbe
teta=linspace(0,pi,100000);

delta1=.475;
ff1=ftr(teta,order,genre,delta1);
max(ff1)

delta2=.45;
ff2=ftr(teta,order,genre,delta2);
max(ff2)

delta3=.4;
ff3=ftr(teta,order,genre,delta3);
max(ff3)

delta4=.3;
ff4=ftr(teta,order,genre,delta4);
max(ff4)

delta5=0;
genre2='redonnet';
ff5=ftr(teta,order,genre2,delta5);
max(ff5)

figure(1)
plot(teta,ff1,teta,ff2,teta,ff3,teta,ff4,teta,ff5,'--','Linewidth',2)
legend(['alfa_f = ',num2str(delta1)], ['alfa_f = ',num2str(delta2)], ['alfa_f = ',num2str(delta3)], ['alfa_f = ',num2str(delta4)], ['explicite filter'],'Location','southwest')
title(['filter type : ', genre, ', with order : ', num2str(order)])
grid on
axis([0 3.14 0 1.05])