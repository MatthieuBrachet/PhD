%% test;
clc; clear all; close all;


a=0;
b=pi/2;
n=100;
h=(b-a)/n;
radius=6371220;
omega=7.2920e-5;

fa=radius*funfu(a).*(omega*sin(a)+tan(a)/radius.*funfu(a));
fb=radius*funfu(b).*(omega*sin(b)+tan(b)/radius.*funfu(b));

int=0.5*(fa+fb).*h;
for kk=1:n-1
    xx=a+kk*h;
    fxx=radius*funfu(xx).*(omega*sin(xx)+tan(xx)/radius.*funfu(xx));
    int=int+fxx.*h;
end
int=int
