function [y] = fun11(x)
global radius teta0 teta1 omega u0
 y=radius*u0*fun10(x, teta0, teta1).*(2*omega*sin(x)+u0/radius*tan(x).*fun10(x, teta0, teta1));
 
 