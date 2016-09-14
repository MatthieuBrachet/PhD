function [cc] = coriolis_mat74(u,x,y,z)
global radius
[n1,n2]=size(x);
for i=1:n1
    for j=1:n2
        vect(1:3)=u(i,j,1:3);
        kx=x(i,j)./radius;
        ky=y(i,j)./radius;
        kz=z(i,j)./radius;
        COR=[0 -kz ky;kz 0 -kx; -ky kx 0 ];
        cc(i,j,1:3)=COR*vect';
    end
end

end

