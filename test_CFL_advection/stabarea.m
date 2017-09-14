clc; clear all; close all;

x=-3:0.001:3;
[X,Y]=meshgrid(x,x);

L=(X+1i*Y);
AMPLI=1+L+.5*L.^2+(L.^3)/6+(L.^4)/24;
AMPLI=abs(AMPLI);

figure(1)
contourf(X,Y,AMPLI.*(AMPLI<1))


x=0:0.01:3;
L=1i*x;
AMPLI=1+L+.5*L.^2+(L.^3)/6+(L.^4)/24;
AMPLI=abs(AMPLI);

figure(2)
plot(x,AMPLI.*(AMPLI<1))