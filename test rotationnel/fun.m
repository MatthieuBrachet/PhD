function [ h ] = fun(x,y,z)
global radius
% p=3;
% q=2;
% r=3;
% h=(x/radius).^p.*(y/radius).^q.*(z/radius).^r;
[lambda, teta, rr]=cart2sph(x,y,z);
h=cos(teta).^5.*sin(30.*lambda);