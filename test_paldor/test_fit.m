clc; clear all; close all;

n=200;
x=linspace(0,10,n)';
alp=2;
y=cos(alp*x)+rand(size(x));

f = fit(x,y,'fourier2');
Cn=f.w

figure(1)
hold on
plot(x,y)
plot(f)
hold off