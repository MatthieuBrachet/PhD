function [ h ] = fun(x,y,z)
global radius
p=3;
q=2;
r=4;
h=(x/radius).^p.*(y/radius).^q.*(z/radius).^r;
%h=(x).^p.*(y).^q.*(z).^r;
