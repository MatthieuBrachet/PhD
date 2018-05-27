% PLOT THE CUBED SPHERE GRID

clear all;
global n nn;

N=128;
make_cs_grid(N);

radius=1;r=radius;
A=[];B=[];C=[];
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
A=[A,x];B=[B,y];C=[C,z];
%
% FACE E - II;
y=sqrt(radius^2./delta);
x=-xx.*y;
z=yy.*y;
A=[A,x];B=[B,y];C=[C,z];
%
% FACE B - III;
x=-sqrt(radius^2./delta);
y=xx.*x;
z=-yy.*x;
A=[A,x];B=[B,y];C=[C,z];
%
% FACE W - IV;
y=-sqrt(radius^2./delta);
x=-xx.*y;
z=-yy.*y;
A=[A,x];B=[B,y];C=[C,z];
%
% FACE N - V;
z=sqrt(radius^2./delta);
y=xx.*z;
x=-yy.*z;
A=[A,x];B=[B,y];C=[C,z];
%
% FACE S - VI;
z=-sqrt(radius^2./delta);
y=-xx.*z;
x=-yy.*z;
A=[A,x];B=[B,y];C=[C,z];
%
% MESH
nhs=5; mhs=2;
fun=real(sph(nhs,mhs,A,B,C));
mesh(A,B,C,fun);
colorbar;
axis([-r r -r r -r r]);
text(1,0,0,'F');
text(0,1,0,'E');
text(-1,0,0,'B');
text(0,-1,0,'W');
text(0,0,1,'N');
text(0,0,-1,'S');
%
axis equal;
