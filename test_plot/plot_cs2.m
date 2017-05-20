function []=plot_cs2(n,nn)
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
%
% PLOT OF A POLYNOMIAL FUNCTION ON THE SPHERE
% FACE F - I;
x1=sqrt(radius^2./delta);
y1=xx.*x1;
z1=yy.*x1;
[fun,funx,funy,funz]=fun3(x1,y1,z1);
% figure(3);
surf(x1,y1,z1,fun);hold on;
axis([-r r -r r -r r]);
text(1,0,0,'F');
%
% FACE E - II;
y2=sqrt(radius^2./delta);
x2=-xx.*y2;
z2=yy.*y2;
[fun,funx,funy,funz]=fun3(x2,y2,z2);
% figure(3);
surf(x2,y2,z2,fun);hold on;
axis([-r r -r r -r r]);
text(0,1,0,'E');
%
% FACE B - III;
x3=-sqrt(radius^2./delta);
y3=xx.*x3;
z3=-yy.*x3;
[fun,funx,funy,funz]=fun3(x3,y3,z3);
% figure(3);
surf(x3,y3,z3,fun);hold on;
axis([-r r -r r -r r]);
text(-1,0,0,'B');
%
% FACE W - IV;
y4=-sqrt(radius^2./delta);
x4=-xx.*y4;
z4=-yy.*y4;
[fun,funx,funy,funz]=fun3(x4,y4,z4);
% figure(3);
surf(x4,y4,z4,fun);hold on;
axis([-r r -r r -r r]);
text(0,-1,0,'W');
%
% FACE N - V;
z5=sqrt(radius^2./delta);
y5=xx.*z5;
x5=-yy.*z5;
[fun,funx,funy,funz]=fun3(x5,y5,z5);
% figure(3);
surf(x5,y5,z5,fun);hold on;
axis([-r r -r r -r r]);
text(0,0,1,'N');
%
% FACE S - VI;
z6=-sqrt(radius^2./delta);
y6=-xx.*z6;
x6=-yy.*z6;
[fun,funx,funy,funz]=fun3(x6,y6,z6);
% figure(3);
surf(x6,y6,z6,fun);hold on;
axis([-r r -r r -r r]);
text(0,0,-1,'S');
%
axis equal;

