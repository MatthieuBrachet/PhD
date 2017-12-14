function [ h ] = fun(x,y,z)
global radius
p=3;
q=2;
r=4;
h=x.^p.*y.^q.*z.^r;
h=h./radius^(p+q+r);
