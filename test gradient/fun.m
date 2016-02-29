function [h]=fun(x,y,z)
[l,t,~]=cart2sph(x,y,z);
h=exp(l)+exp(t);
end