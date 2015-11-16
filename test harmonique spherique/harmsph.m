function [ h ] = harmsph( m,l,teta,lambda )
% ********************************************
% Spherical Harmonic
%
% Author :
%     - Matthieu Brachet
% ********************************************
nom=(2*l+1)*factorial(l-abs(m));
denom=4*pi*factorial(l+abs(m));
Kml=sqrt(nom./denom);
if (m>0)
    PP=legendre(m,l,sin(teta));
    h=sqrt(2).*Kml.*cos(m*lambda).*PP;
elseif (m==0)
    PP=legendre(0,l,sin(teta));
    h=sqrt((2*l+1)/(4*pi))*PP;
elseif (m<0)
    PP=legendre(-m,l,sin(teta));
    h=sqrt(2).*Kml.*sin(-m*lambda).*PP;
end
end