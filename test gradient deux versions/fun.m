function [func] = fun(x,y,z)
global radius
xx=x./radius; 
yy=y./radius;
zz=z./radius;
func=xx.^3+yy.*zz+cos(xx).*zz.^7;
end

