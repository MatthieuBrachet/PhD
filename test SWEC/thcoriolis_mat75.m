function [f] = thcoriolis_mat75(u,alpha,x,y,z)

global radius omega
[n1,n2]=size(x);
[~, teta,~]=cart2sph(x,y,z);
for i=1:n1
    for j=1:n2
        vect(1:3)=u(i,j,1:3);
        kx=x(i,j)./radius;
        ky=y(i,j)./radius;
        kz=z(i,j)./radius;
        COR=[0 -kz ky;kz 0 -kx; -ky kx 0 ];
        cc(i,j,1:3)=COR*vect';
        f(i,j,1:3)=u(i,j,1:3)+alpha.*2.*omega.*sin(teta(i,j)).*cc(i,j,1:3);
    end
end

end
