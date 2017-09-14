clc; clear all; close all
format long

%% schema testé
order=4;
genre='visbal';

%% courbe
teta=linspace(0,pi,10000);

delta1=0.5;
ff1=ftr(teta,order,genre,delta1);
max(ff1)

delta2=0.48;
ff2=ftr(teta,order,genre,delta2);
max(ff2)

delta3=0.25;
ff3=ftr(teta,order,genre,delta3);
max(ff3)

delta4=0.1;
ff4=ftr(teta,order,genre,delta4);
max(ff4)

delta5=0;
genre2='redonnet';
ff5=ftr(teta,order,genre2,delta5);
max(ff5)

figure(1)
plot(teta,ff1,teta,ff2,teta,ff3,teta,ff4,teta,ff5,'--')
legend(['alfa_f = ',num2str(delta1)], ['alfa_f = ',num2str(delta2)], ['alfa_f = ',num2str(delta3)], ['alfa_f = ',num2str(delta4)], ['explicite filter'] )
title(['filter type : ', genre, ', with order : ', num2str(order)])
grid on