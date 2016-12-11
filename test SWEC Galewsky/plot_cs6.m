function []=plot_cs6(n,nn)
% PLOT THE CUBED SPHERE GRID
% J-P. CROISILLE - SEP 20, 2012
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
%% PLOT OF THE PROJECTION OF THE GRID ON FACE I ON THE CUBE.
%
yc=zeros(nn);
zc=zeros(nn);
%
yc=tan(xi);
zc=tan(eta);
surf(yc,zc);


