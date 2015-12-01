function [ vect ] = coriolis(V,x,y,z)
% calcul the coriolis force with the cross product :
%      vect = f.k x V
% k is the unit vector normal to the sphere at the coordinate (x,y,z).
% V is a nxmx3 tabular where V(:,:,i) is the i-th component of a vector V.
global OMEGA

[lambda, phi,~]=cart2sph(x,y,z);
k(:,:,1)=cos(phi).*sin(lambda);
k(:,:,2)=cos(phi).*sin(lambda);
k(:,:,3)=sin(phi);
vect=zeros(size(V));
cor=2.*OMEGA.*sin(phi);
for i=1:size(V,1)
    for j=1:size(V,2)
        a=[k(i,j,:)];
        b=[V(i,j,:)];
        v=cross(a,b);
        vect(i,j,:)=cor(i,j).*v;
    end
end
end