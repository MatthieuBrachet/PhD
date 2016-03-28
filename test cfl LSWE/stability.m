clear all; close all; clc;
format long
timescheme='rk4';

teta1=linspace(0,pi,50);
teta2=teta1;
[TETA1,TETA2]=meshgrid(teta1,teta2);

dx=0.01;
dy=0.01;
dt=0.0039;

[ stab ] = ampli(teta1, teta2,dx,dy,dt,timescheme);

figure(1)
surf(TETA1,TETA2,stab)

MAX=max(max(stab))
MIN=min(min(stab))