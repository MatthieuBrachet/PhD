function [v]=fun5b(x)
% 1/x
%----------------
global radius;
%----------------
n1=size(x,1);
n2=size(x,2);
v=zeros(n1,n2);
%---------------
% -Choix 1
% te0=2.028757838110434; % cf t40.m
% %
% k1=2/(radius*sin(te0));
% k2=2*te0/radius;
% v=k1*sin(k2*x);
% % - Choix 2
% polphi=[48 -64 24 0];
% ax=abs(x);
% ax1=max(ax,0.5);
% %%%!!! ATTENTION RADIUS=1 ICI. !!!!
% v=sign(x).*(polyval(polphi,ax).*(ax<=0.5)+ (1./ax1).*(ax>0.5)); 
% - Choix 3:
polphi=[64 0 -48 0 12 0];
ax=abs(x);
ax1=max(ax,0.5);
%%%!!! ATTENTION RADIUS=1 ICI. !!!!
v=polyval(polphi,x).*(ax<=0.5)+ sign(x).*(1./ax1).*(ax>0.5); 