function [dfu] = dfun1(x)
dfu=(-12*x+(24/sqrt(2)).*x.^2).*(x>=0).*(x<=sqrt(2)/2);
end