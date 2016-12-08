function [y] = fun11b(x)
global radius omega u0
 y=radius*u0*fun10b(x).*(2*omega*sin(x)+u0/radius*tan(x).*fun10b(x));