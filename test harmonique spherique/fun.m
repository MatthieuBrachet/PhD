function [h]=fun(x,y,z,m,l)
[n1,n2]=size(x);
for i=1:n1
    for j=1:n2
        [lambda, teta, ~]=cart2sph(x(i,j),y(i,j),z(i,j));
        h(i,j)= harmsph( m,l,teta,lambda );
    end
end
end