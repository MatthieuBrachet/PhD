clc; clear all; close all;

syms f g s1 dx s2 dy H dt

A=[0 f -g*1i*s1/dx; -f 0 -g*1i*s2/dy; -H*1i*s1/dx -H*1i*s2/dy 0];
Id=eye(3,3);

%% valeurs propres de Id + dt*A

simplify(eig(Id+dt*A))