function [ y ] = u_geo( x )
global radius OMEGA
y=radius*u_test(x)*(2*OMEGA*sin(x)+(tan(x)/radius).*u_test(x));
end

