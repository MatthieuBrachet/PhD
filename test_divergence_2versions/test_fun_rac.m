clc; clear all; close all
global radius
radius=1;
x=-1:.001:1;
y1=1./abs(x);
y2=fun5(x);

figure(1)
plot(x,y1,x,y2)
axis([-1 1 0 50])