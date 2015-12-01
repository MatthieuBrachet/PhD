function [ gp ] = geopotential( x,y,z )
% geopotential on the sphere.
% ***************************
% gp (output) : matrix containing geopotential calculated on the coordinate
% (x,y,z).
%
% geopotentiel = gp = gravity * h
%
% stationnary solution of SWE+Coriolis.
% depend only of the longitud (no latitude).
global g
global h0
global n_simpson

[~, phi, ~]=cart2sph(x,y,z);
gp=zeros(size(phi));
syms x
n=floor(size(phi,1)/2);
nn=size(phi,1);
for i=1:n
    for j=1:size(phi,2)
        pp=phi(i,j);
         [ I ] = simpson(@u_geo,0,pp,n_simpson);
         clc; disp([num2str(i*j), '/', num2str(nn*nn)]);
%             syms x
%             I=double(int(u_geo(x),0,pp));
        gp(i,j)=g*h0-I;
    end
end
gp(n+1:nn,:)=gp(n:-1:1,:);

end