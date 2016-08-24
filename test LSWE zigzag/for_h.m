function [forh] = for_h(x,y,z,t)
% forcaque sur l'équation de conservation
    global radius test
    sigma=1/5000;
if test == 3
    [lambda, teta,~]=cart2sph(x,y,z);
    forh = -sigma.*cos(teta).*sin(lambda).*exp(-sigma*t);
else
    forh=zero(size(x));
end

end