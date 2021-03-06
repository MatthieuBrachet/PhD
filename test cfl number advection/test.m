clc; clear all; close all;

x = linspace(-10,10,200);
y = cos(x);
plot(x,y)

xticks([-3*pi -2*pi -pi 0 pi 2*pi 3*pi])
xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%yticks([-1 -0.8 -0.2 0 0.2 0.8 1])