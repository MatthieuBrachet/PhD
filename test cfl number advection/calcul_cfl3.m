clc; clear all; close all;

cfl=1;
teta=0:0.001:2*pi;
c=4;
space='implicit';
f=10;
filtre='visbal';
delta=0.4;

%% cercle unité
cerc1=exp(1i*teta);

figure(1)
hold on
plot(cerc1,'k-')



%% schema 1
time1='rk4';
g1 = atsf(cfl, teta,time1,c,space,f,filtre,delta);

plot(g1,'bo')

%% schema 2
time2='dirk34a';
g2 = atsf(cfl, teta,time2,c,space,f,filtre,delta);

plot(g2,'ro')

%% schema 3
time3='dirk33b';
g3 = atsf(cfl, teta,time3,c,space,f,filtre,delta);

plot(g3,'go')



%% 
legend('unit circle',time1,time2,time3,2)
axis([-1 1 -1 1])
hold off