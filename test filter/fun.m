function [ h ] = fun(x,y,z)
global radius
test=2;
if test ==1
    h=cos(2*pi*x/radius).*cos(2*pi*y/radius).*cos(2*pi*z/radius)+3;
elseif test == 2
    h=(x/radius).^2.*(z./radius).^4 + y./radius+3;
elseif test == 3
    h=exp(x/radius)+exp(y/radius)+exp(z/radius);
end
