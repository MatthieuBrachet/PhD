function [ h ] = dharmsph_lambda( m,l,teta,lambda )
% ********************************************
% partial derivative of Spherical Harmonic
% in lambda.
%
% Author :
%     - Matthieu Brachet
% ********************************************
teta=pi/2-teta;
phi=lambda;
nom=(2*l+1)*factorial(l-abs(m));
denom=4*pi*factorial(l+abs(m));
Kml=sqrt(nom./denom);
if (m>0)
    PP=legendre(m,l,cos(teta));
    h=m.*sqrt(2).*Kml.*sin(m*phi).*PP;
elseif (m==0)
    h=0;
elseif (m<0)
    PP=legendre(-m,l,cos(teta));
    h=-m.*sqrt(2).*Kml.*cos(-m*phi).*PP;
end
end