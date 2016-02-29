clear all; close all; clc;

teta1=0:0.1:pi;
teta2=0:0.1:pi;
[TETA1,TETA2]=meshgrid(teta1,teta2);

dx=0.01;
dy=0.01;
dt=0.0039;

[ stab ] = ampli(teta1, teta2,dx,dy,dt);

figure(1)
surf(TETA1,TETA2,stab)
