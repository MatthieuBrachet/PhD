clc; clear all; close all;

x=-3:0.01:3;
y=-3:0.01:3;

[X,Y]=meshgrid(x,y);


dx=0.5;
dy=0.5;
h=sin(pi/dx*X).*sin(pi/dy*Y);
h=h./(pi^2*X.*Y/(dx*dy));

surf(X,Y,log(abs(h)))