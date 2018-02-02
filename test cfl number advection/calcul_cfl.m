clc; clear all; close all

%% scheme we want to test
time='rk4';
c=4;
space='implicit';
f=2;
delta=0;
filtre='redonnet';

%% research of cfl number max
epsi=10^-10;
itermax=1000;
iter=0;

%% data for cfl research
teta = linspace(0,pi,100000);
a=0;
b=5;

while iter < itermax && abs(b-a)>eps
    lambda=0.5*(a+b);
    g = atsf(lambda, teta,time,c,space,f,filtre,delta);
    h=floor(max(abs(g))/epsi)*epsi;
    iter=iter+1;
    if h>1
        b=lambda;
    else 
        a=lambda;
    end
end
cfl=lambda
iter

abs(cfl-2*sqrt(2/3))

g = atsf(cfl, teta,time,c,space,f,filtre,delta);
figure(1)
plot(teta,abs(g))
grid on
xlabel('wave number')
ylabel('amplification')
axis([0 pi 0 1.1])