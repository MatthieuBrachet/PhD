function [ h ] = fun(x,y,z)
global radius
% p=3;
% q=2;
% r=3;
% h=(x/radius).^p.*(y/radius).^q.*(z/radius).^r;
z=z/radius;
h=1./(101-100.*z);