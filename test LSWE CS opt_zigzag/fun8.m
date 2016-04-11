function [v]=fun8(x,y,z,t)
global nn
global radius

v=zeros(nn,nn,4);
n1=nn;
n2=nn;
global kv keta;
kv=1;keta=1;
kw1=1;kw2=1;kw3=1;
gcv(1:n1,1:n2,1)=kv*kw1*(pi/radius)*cos(kw1*pi*x/radius);
gcv(1:n1,1:n2,2)=kv*kw2*(pi/radius)*cos(kw2*pi*y/radius);
gcv(1:n1,1:n2,3)=kv*kw3*(pi/radius)*cos(kw3*pi*z/radius);
 %
xnorm=zeros(n1,n2,3);
nnorm=sqrt(x.*x+y.*y+z.*z);
xnorm(1:n1,1:n2,1)=x./nnorm;
xnorm(1:n1,1:n2,2)=y./nnorm;
xnorm(1:n1,1:n2,3)=z./nnorm;
% SPHERICAL GRADIENT 
v(1:n1,1:n2,1)=gcv(1:n1,1:n2,1)...
                   -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,1);
v(1:n1,1:n2,2)=gcv(1:n1,1:n2,2)...
                   -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,2);
v(1:n1,1:n2,3)=gcv(1:n1,1:n2,3)...
                   -dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,3);

v(1:n1,1:n2,4)=keta*(sin(kw1*pi*x/radius)+sin(kw2*pi*y/radius)+sin(kw3*pi*z/radius));