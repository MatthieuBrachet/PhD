function [v,dv]=fun6(x,y,z)
global radius
n1=size(x,1);
n2=size(x,2);

%% vector field
phi(1:n1,1:n2,1)=1;
phi(1:n1,1:n2,2)=2;
phi(1:n1,1:n2,3)=3;

nnorm=radius;
xnorm(1:n1,1:n2,1)=x./nnorm;
xnorm(1:n1,1:n2,2)=y./nnorm;
xnorm(1:n1,1:n2,3)=z./nnorm;

v1(1:n1,1:n2,1)=xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,3)...
    - xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,2);
v1(1:n1,1:n2,2)=xnorm(1:n1,1:n2,3).*phi(1:n1,1:n2,1)...
    - xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,3);
v1(1:n1,1:n2,3)=xnorm(1:n1,1:n2,1).*phi(1:n1,1:n2,2)...
    - xnorm(1:n1,1:n2,2).*phi(1:n1,1:n2,1);

a=1;
kw1=2;kw2=6;kw3=10;
cv=a*(sin(kw1*pi*x)+sin(kw2*pi*y)+sin(kw3*pi*z));

v(1:n1,1:n2,1)=cv(1:n1,1:n2).*v1(1:n1,1:n2,1);
v(1:n1,1:n2,2)=cv(1:n1,1:n2).*v1(1:n1,1:n2,2);
v(1:n1,1:n2,3)=cv(1:n1,1:n2).*v1(1:n1,1:n2,3);

%% divergence associated
gcv(1:n1,1:n2,1)=a*kw1*pi*cos(kw1*pi*x);
gcv(1:n1,1:n2,2)=a*kw2*pi*cos(kw2*pi*y);
gcv(1:n1,1:n2,3)=a*kw3*pi*cos(kw3*pi*z);

gscv(1:n1,1:n2,1)=gcv(1:n1,1:n2,1)-dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,1);
gscv(1:n1,1:n2,2)=gcv(1:n1,1:n2,2)-dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,2);
gscv(1:n1,1:n2,3)=gcv(1:n1,1:n2,3)-dot(xnorm(1:n1,1:n2,1:3),gcv(1:n1,1:n2,1:3),3).*xnorm(1:n1,1:n2,3);
           
dv(1:n1,1:n2)=dot(gscv(1:n1,1:n2,1:3),v1(1:n1,1:n2,1:3),3);