function [forh] = for_h(x,y,z,t)
% forcaque sur l'équation de conservation
global radius test
if test == 3
    sigma=10^-4;
    [lambda, teta,~]=cart2sph(x,y,z);
    forh = -sigma.*cos(teta).*sin(lambda).*exp(-sigma*t);
    
elseif test ==4
     sigma=10^-4;
    [lambda, teta,~]=cart2sph(x,y,z);
    %teta=teta+pi;
    pwa=1;
    pwb=5;
    forh = sigma.*(cos(teta).^pwb).*(cos(3*lambda).^pwa).*cos(sigma*t);   
else
    forh=zeros(size(x));
end

end