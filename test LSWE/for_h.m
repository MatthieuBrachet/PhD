function [forh] = for_h(x,y,z,t)
% forcaque sur l'équation de conservation
global radius test teta0 teta1
if test == 3
    sigma=10^-4;
    [lambda, teta,~]=cart2sph(x,y,z);
    forh = -sigma.*cos(teta).*sin(lambda).*exp(-sigma*t);
    
elseif test ==4
     sigma=10^-4;
    [lambda, teta,~]=cart2sph(x,y,z);
    pwa=1;
    pwb=5;
    forh = sigma.*(cos(teta).^pwb).*(cos(3*lambda).^pwa).*cos(sigma*t);   
    
elseif test ==5
     sigma=10^-4;
    [~, teta,~]=cart2sph(x,y,z);
%     forh = -sigma.*fun10(teta,-pi/4,pi/4).*exp(-sigma*t); 
      forh = -sigma.*fun10(teta,teta0,teta1).*exp(-sigma*t); 


else
    forh=zeros(size(x));
end

end