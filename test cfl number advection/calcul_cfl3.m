clc; clear all; close all;

cfl=2*sqrt(2/3);
teta=0:0.001:4*pi;
c=6;
space='implicit';
f=0;
filtre='redonnet';
delta=0.4;

%% cercle unité
cerc1=exp(1i*teta);

figure(1)
hold on
plot(cerc1,'k-','Linewidth',2)



%% schema 1
time1='rk4';
g1 = atsf(cfl, teta,time1,c,space,f,filtre,delta);

plot(g1,'b-','Linewidth',2)

%% schema 2
time2='dirk34a';
g2 = atsf(cfl, teta,time2,c,space,f,filtre,delta);

plot(g2,'r-','Linewidth',2)

%% schema 3
time3='dirk33b';
g3 = atsf(cfl, teta,time3,c,space,f,filtre,delta);

plot(g3,'g-','Linewidth',2)



%% 
legend('unit circle',time1,time2,time3,'Location','northwest')
axis([-1 1 -1 1])
hold off