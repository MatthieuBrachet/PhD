function [x,f] = coupe_eq(funfI,funfII,funfIII,funfIV)
global x_fII y_fII z_fII
[lambda, ~, ~]=cart2sph(x_fII, y_fII, z_fII);
p=floor(size(lambda,1)/2);
f=funfII(:,p);
x=lambda(:,p);
end