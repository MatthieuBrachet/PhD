function [ h ] = dharmsph_theta( m,l,teta,lambda )
% ********************************************
% partial derivative of Spherical Harmonic
% in theta.
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
    dPP=derlegendre(m,l,cos(teta));
    h=-sqrt(2).*Kml.*cos(m*phi).*sin(teta).*dPP;
elseif (m==0)
    dPP=derlegendre(0,l,cos(teta));
    h=-sin(teta).*sqrt((2*l+1)/(4*pi))*dPP;
elseif (m<0)
    dPP=derlegendre(-m,l,cos(teta));
    h=-sin(teta).*sqrt(2).*Kml.*sin(-m*phi).*dPP;
end
end