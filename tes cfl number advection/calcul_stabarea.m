clc; clear all; close all;

lambda=0.9;
teta=-10:0.001:10;

c=4;
space='implicit';

figure(1);
hold on
grid on
%% without filter
f=0;
filtre='explicit';
delta=0;
[ g ] = ark4cf(lambda, teta,c,space,f,filtre,delta);
rk4c4f0=log(g);
plot(rk4c4f0,'bx');

%% filter order 10
f=10;
filtre='explicit';
delta=0;
[ g ] = ark4cf(lambda, teta,c,space,f,filtre,delta);
rk4c4f10=log(g);
plot(rk4c4f10,'rx');

%% filter order 6
f=6;
filtre='explicit';
delta=0;
[ g ] = ark4cf(lambda, teta,c,space,f,filtre,delta);
rk4c4f6=log(g);
plot(rk4c4f6,'gx');

%% filter order 2
f=2;
filtre='explicit';
delta=0;
[ g ] = ark4cf(lambda, teta,c,space,f,filtre,delta);
rk4c4f2=log(g);
plot(rk4c4f2,'kx')



legend('rk4c4f0','rk4c4f10','rk4c4f4','rk4c4f2',2)

