clc; clear all; close all

%% scheme we want to test
time='rk1';
c=4;
space='implicit';
f=8;
delta=0.1;
filtre='redonnet';

%% research of cfl number max
epsi=10^-8;
itermax=1000;
iter=0;

%% data for cfl research
teta = 0:0.001:pi;
a=0;
b=5;

while abs(b-a)>epsi && iter < itermax
    lambda=0.5*(a+b);
    g = atsf(lambda, teta,time,c,space,f,filtre,delta);
    h=max(abs(g));
    iter=iter+1;
    if h>1
        b=lambda;
    else 
        a=lambda;
    end
end
cfl=lambda
iter

g = atsf(lambda, teta,time,c,space,f,filtre,delta);
figure(1)
plot(teta,abs(g))
grid on
xlabel('wave number')
ylabel('amplification')
axis([0 pi 0 1.1])