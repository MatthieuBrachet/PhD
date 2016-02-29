clc; clear all; close all

%% scheme we want to test
time='rk1';
c=4;
space='implicit';
f=10;
delta=0.4;
filtre='visbal';

%% research of cfl number max
epsi=10^-10;
itermax=1000;
iter=0;

%% data for cfl research
teta = 0:0.001:pi;
a=0;
b=5;

while abs(b-a)>epsi && iter < itermax
    lambda=0.5*(a+b);
    g = atsf(lambda, teta,time,c,space,f,filtre,delta);
    h=floor(max(abs(g))*1/epsi)*epsi;
    %h=max(abs(g));
    iter=iter+1;
    if h>1
        b=lambda;
    else 
        a=lambda;
    end
end
cfl=lambda
iter

g = atsf(cfl, teta,time,c,space,f,filtre,delta);
figure(1)
plot(teta,abs(g))
grid on
xlabel('wave number')
ylabel('amplification')
axis([0 pi 0 1.1])


lambda=linspace(0,5,1000);
for k=1:size(lambda,2)
    g = atsf(lambda(k), teta,time,c,space,f,filtre,delta);
    h(k)=max(abs(g));
end

figure(10)
plot(lambda,h)
grid on