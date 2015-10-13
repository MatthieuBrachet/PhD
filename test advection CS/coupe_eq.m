function [x,f] = coupe_eq(funfI,funfII,funfIII,funfIV)
global x_fI y_fI z_fI
global x_fII y_fII z_fII
global x_fIII y_fIII z_fIII
global x_fIV y_fIV z_fIV

% [lambda, ~, ~]=cart2sph(x_fI, y_fI, z_fI);
% p=floor(size(lambda,1)/2);
% f1=funfI(:,p);
% 
% [lambda, ~, ~]=cart2sph(x_fII, y_fII, z_fII);
% p=floor(size(lambda,1)/2);
% f2=funfII(:,p);
% 
% [lambda, ~, ~]=cart2sph(x_fIII, y_fIII, z_fIII);
% p=floor(size(lambda,1)/2);
% f3=funfIII(:,p);
% 
% [lambda, ~, ~]=cart2sph(x_fIV, y_fIV, z_fIV);
% p=floor(size(lambda,1)/2);
% f4=funfIV(:,p);
% 
% f=[f1; f2; f3; f4];
% n=length(f);
% x=linspace(-pi/4,7*pi/4,n);

[lambda, ~, ~]=cart2sph(x_fII, y_fII, z_fII);
p=floor(size(lambda,1)/2);
f=funfII(:,p);
x=lambda(:,p);
end