function [y] = fun12(x)
global radius omega
[ up ] = u_test3(x);
y=(2.*omega.*sin(x)+up.*tan(x)./radius).*up;