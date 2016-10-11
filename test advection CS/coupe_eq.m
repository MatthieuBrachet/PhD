function [x,f] = coupe_eq(funfI,funfII,funfIII,funfIV)
global x_fI y_fI z_fI
[lambda, ~, ~]=cart2sph(x_fI, y_fI, z_fI);
p=floor(size(lambda,1)/2);
f=funfI(:,p);
x=lambda(:,p);
end