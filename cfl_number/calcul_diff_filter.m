clc; clear all; close all

%% schema testé
order=10;
genre='visbal';

%% courbe
teta=0:0.001:pi;

delta1=0.49;
ff1=ftr(teta,order,genre,delta1);

delta2=0.4;
ff2=ftr(teta,order,genre,delta2);

delta3=0.3;
ff3=ftr(teta,order,genre,delta3);

delta4=0.1;
ff4=ftr(teta,order,genre,delta4);

delta5=0;
genre2='redonnet';
ff5=ftr(teta,order,genre2,delta5);


figure(1)
plot(teta,ff1,teta,ff2,teta,ff3,teta,ff4,teta,ff5,'--')
legend(['alfa_f = ',num2str(delta1)],['alfa_f = ',num2str(delta2)],['alfa_f = ',num2str(delta3)],['alfa_f = ',num2str(delta4)],['explicite filter'] ,3)
title(['filter type : ', genre, ', with order : ', num2str(order)])
grid on