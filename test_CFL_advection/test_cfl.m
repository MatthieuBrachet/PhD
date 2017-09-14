clc; clear all; close all;
format long

teta=linspace(0,pi,1000);
lambda=1.689;

cp=lambda*1i* 3*sin(teta)./(2+(cos(teta)));
ftr = 772/1024+2*(210/1024*cos(teta)-120/1024*cos(2*teta)+45/1024*cos(3*teta)-10/1024*cos(4*teta)+1/1024*cos(5*teta));
ampli=(1+cp+.5*cp.^2+cp.^3/6+cp.^4/24).*ftr;

figure(1)
plot(teta,abs(ampli))

max(abs(ampli))