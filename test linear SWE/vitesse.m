function [ V ] = vitesse( x,y,z )
[lambda, phi, ~]=cart2sph(x,y,z);
u=u_test(phi);
V(:,:,1)=-u.*sin(lambda);
V(:,:,2)=u.*cos(lambda);
V(:,:,3)=u.*0;
end