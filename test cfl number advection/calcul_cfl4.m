clc; clear all; close all;

lambda=0.7;
teta=0:0.01:pi;

%% schema entièrement compact
time='dirk12a';
c=4;
space='explicit';
f=10;
filtre='redonnet';
delta=0;
[ g1 ] = atsf(lambda, teta,time,c,space,f,filtre,delta);
h1=abs(g1);

%% schema compact/classique
ftr = 772/1024+2*(210/1024*cos(teta)-120/1024*cos(2*teta)+45/1024*cos(3*teta)-10/1024*cos(4*teta)+1/1024*cos(5*teta));
%si=4/3*sin(teta)-1/3*sin(2*teta);
si=3/2*sin(teta)-3/5*sin(2*teta)+1/10*sin(3*teta);
s = sin(teta)./(2/3+2*(1/6*cos(teta)));
thetai=-1i.*lambda.*si;
theta=-1i.*lambda.*s;
g = (1+theta/2)./(1-thetai/2);
g = g.*ftr;
g2=g.*ftr;
h2=abs(g2);

figure(1)
plot(teta,h1,teta,h2)
legend('entièrement compact','partiellement compact')