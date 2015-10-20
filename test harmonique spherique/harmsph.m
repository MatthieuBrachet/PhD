function [ h ] = harmsph( m,l,teta,lambda )
% ********************************************
% Spherical Harmonic
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
    h=sqrt(2).*Kml.*cos(m*phi).*PP;
elseif (m==0)
    PP=legendre(0,l,cos(teta));
    h=sqrt((2*l+1)/(4*pi))*PP;
elseif (m<0)
    PP=legendre(-m,l,cos(teta));
    h=sqrt(2).*Kml.*sin(-m*phi).*PP;
end
end