function []=plot_cs3(n,nn)
% PLOT THE CUBED SPHERE GRID
% J-P. CROISILLE - JULY 20 2010
% clear all;
% n=3,nn=5; 
% n=7, nn=9; % number of points by face
% n=15, nn=17; % number of points by face
%  n=31, nn=33; % number of points by face
% n =63, nn=65; % number of points by face
 % n=127, nn=129; % number of points by face
% n=255, nn=257; % number of points by face
% n=511, nn=513; % number of points by face
%
radius=1;r=radius;
%
xi=linspace(-pi/4, pi/4, nn); 
eta=linspace(-pi/4, pi/4, nn); 
for i=1:nn,
  for j=1:nn,
    xx(i,j)=tan(xi(i));
    yy(i,j)=tan(eta(j));
  end
end
%
for i=1:nn,
  for j=1:nn,
    delta(i,j)=1+xx(i,j)^2+yy(i,j)^2;
  end
end
%
%% PLOT OF THE BOUNDARY OF THE ALFA/BETA DOMAIN FOR ONE FACE
%%
alfaL=-atan(cos(eta));
betaL=atan(tan(eta)/sqrt(2));
%
alfaR=atan(cos(eta));
betaR=atan(tan(eta)/sqrt(2));
%
alfaB=atan(tan(xi)/sqrt(2));
betaB=-atan(cos(xi));
%
alfaT=atan(tan(xi)/sqrt(2));
betaT=atan(cos(xi));
%
plot(alfaL,betaL,alfaR,betaR,alfaB,betaB,alfaT,betaT); grid;
%



axis equal;

