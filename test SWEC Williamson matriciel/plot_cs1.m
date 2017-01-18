function []=plot_cs1(n,nn)
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
% FACE F - I;
x=sqrt(radius^2./delta);
y=xx.*x;
z=yy.*x;
% figure(2);
surf(x,y,z); colormap(gray);hold on;
axis([-r r -r r -r r]);
text(1,0,0,'F');
%
% FACE E - II;
y=sqrt(radius^2./delta);
x=-xx.*y;
z=yy.*y;
% figure(2);
surf(x,y,z); colormap(gray);hold on;
axis([-r r -r r -r r]);
text(0,1,0,'E');
%
% FACE B - III;
x=-sqrt(radius^2./delta);
y=xx.*x;
z=-yy.*x;
% figure(2);
surf(x,y,z); hold on;
axis([-r r -r r -r r]);
text(-1,0,0,'B');
%
% FACE W - IV;
y=-sqrt(radius^2./delta);
x=-xx.*y;
z=-yy.*y;
% figure(2);
surf(x,y,z); hold on;
axis([-r r -r r -r r]);
text(0,-1,0,'W');
%
% FACE N - V;
z=sqrt(radius^2./delta);
y=xx.*z;
x=-yy.*z;
% figure(2);
surf(x,y,z); hold on;
axis([-r r -r r -r r]);
text(0,0,1,'N');
%
% FACE S - VI;
z=-sqrt(radius^2./delta);
y=-xx.*z;
x=-yy.*z;
% figure(2);
surf(x,y,z); hold on;
axis([-r r -r r -r r]);
text(0,0,-1,'S');
%
axis equal;
