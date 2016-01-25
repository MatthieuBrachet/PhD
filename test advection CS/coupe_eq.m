function [x,f] = coupe_eq(funfI,funfII,funfIII,funfIV)
global x_fIV y_fIV z_fIV
[lambda, ~, ~]=cart2sph(x_fIV, y_fIV, z_fIV);
p=floor(size(lambda,1)/2);
f=funfIV(:,p);
x=lambda(:,p);
end